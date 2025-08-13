#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <mutex>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>
#include <ctime>
#include <filesystem>

namespace po = boost::program_options;
using namespace std;

#define bits_per_nt 4
#define str_to_bit_factor 2
#define barcode_length 16
#define unassigned_barcode "NA"
#define unassigned_bitcode boost::dynamic_bitset<>(barcode_length*bits_per_nt)
#define max_dist 1
#define max_num_outfiles 64 //maximum number of supported outfiles for static array size, since each outfile creates one full read/write, the limit of 64 should be enough
std::mutex assignment_mutex; // Mutex for thread-safe access to shared data


boost::dynamic_bitset<> string_to_bool_vector(string barcode){
    boost::dynamic_bitset<> bitcode(bits_per_nt*barcode_length);
    for (int i = 0;i<barcode_length;i++){
        switch (barcode[i])
        {
        case 'A':
            bitcode[bits_per_nt*i] = 1;
            break;
        case 'C':
            bitcode[bits_per_nt*i+1] = 1;
            break;
        case 'G':
            bitcode[bits_per_nt*i+2] = 1;
            break;
        case 'T':
            bitcode[bits_per_nt*i+3] = 1;
            break;
        default:
            break;
        }
    }
    return bitcode;
}

string bool_vector_to_string(boost::dynamic_bitset<> bitcode){
    string barcode = "";
    if (bitcode==unassigned_bitcode){
        return unassigned_barcode;
    }
    for(int i=0;i<barcode_length;i++){
        if (bitcode[bits_per_nt*i]){
            barcode += "A";
        } else if (bitcode[bits_per_nt*i+1]) {
            barcode += "C";
        } else if (bitcode[bits_per_nt*i+2]) {
            barcode += "G";
        } else if (bitcode[bits_per_nt*i+3]) {
            barcode += "T";
        } else {
            return unassigned_barcode;
        }
    }
    return barcode;
}

string find_closest_match(string barcode,set<boost::dynamic_bitset<>> &accepted_bitcodes){
    boost::dynamic_bitset<> best_match = unassigned_bitcode;
    int best_distance = max_dist + 1;
    int distance;
    boost::dynamic_bitset<> bitcode = string_to_bool_vector(barcode);
    for (boost::dynamic_bitset<> ref_bitcode : accepted_bitcodes){
        distance = (bitcode ^ ref_bitcode).count();
        if (distance < best_distance){
            best_distance = distance;
            best_match = ref_bitcode;
        } else if (distance == best_distance) {
            //We only want to assign if we have exactly one acceptable matching barcode
            best_match = unassigned_bitcode;
        }
    }
    return bool_vector_to_string(best_match);
}

void process_lines_chunk(
    vector<string>::iterator start,
    vector<string>::iterator end,
    unordered_map<string, vector<string>> &assigned_cells,
    set<string> &accepted_barcodes,
    set<boost::dynamic_bitset<>> &accepted_bitcodes
) {
    string line;
    string read_id;
    string barcode;
    string assigned_to;
    for (vector<string>::iterator it=start; it != end;it++) {
        line = *it;
        read_id = line.substr(0,line.find(' '));
        barcode = read_id.substr(read_id.size()-barcode_length,barcode_length);
        
        if (accepted_barcodes.find(barcode) != accepted_barcodes.end()){
            //exact match, no need for search
            assigned_to = barcode;
        } else if ((barcode.find('N') != string::npos) | (read_id[0] != '@')) {
            //contains N, skip
            assigned_to = unassigned_barcode;
        } else {
            //potentially rescueable
            assigned_to = find_closest_match(barcode,accepted_bitcodes);
        }
       
        // Lock the mutex before modifying shared data
        std::lock_guard<std::mutex> lock(assignment_mutex);
        assigned_cells[assigned_to].push_back(read_id);
        
    }
}

void assign_all_reads(
    unordered_map<string,vector<string>> &assigned_cells,
    set<boost::dynamic_bitset<>> &accepted_bitcodes,
    set<string> &accepted_barcodes,
    string idfile_path,
    int num_threads
) {
    for (string barcode:accepted_barcodes){
        assigned_cells[barcode] = vector<string>();
    }
    assigned_cells[unassigned_barcode] = vector<string>();
    ifstream idfile;
    idfile.open(idfile_path);
    vector<boost::thread> active_threads;
    vector<string>::iterator start;
    vector<string>::iterator end;
    vector<string> contents;
    string line;
    while (getline(idfile,line)) {
        contents.push_back(line);
    }
    size_t chunk_size = contents.size()/num_threads;
    for (int i = 0; i < num_threads; ++i){
        start = contents.begin() + i* chunk_size;
        if (i == num_threads-1) {
            end = contents.end();
        } else {
            end = start + chunk_size;
        }
        active_threads.emplace_back(
            boost::thread(
                process_lines_chunk,
                start,
                end,
                std::ref(assigned_cells),
                std::ref(accepted_barcodes),
                std::ref(accepted_bitcodes)
            )
        );
    }
    for (boost::thread &thread:active_threads){
        thread.join();
    }
}

void parse_barcodes_file(set<string> &legal_barcodes,set<boost::dynamic_bitset<>> &legal_bitcodes,string barcodes_path){
    ifstream infile;    
    infile.open(barcodes_path);
    string line;
    string barcode;
    while (getline(infile,line)){
        barcode = line.substr(0,barcode_length);
        legal_barcodes.insert(barcode);
        legal_bitcodes.insert(string_to_bool_vector(barcode));
    }

}

void write_outfiles(
    unordered_map<string,vector<string>> assigned_cells,
    string outstem,
    string unassigned_file,
    int num_outfiles,
    size_t padding_size
){
    ofstream OutFile[max_num_outfiles];
    ofstream u_file;
    string padded_number;
    string outfile_path;
    for (int i=1;i<=num_outfiles;++i){
        padded_number = to_string(i);
        while (padded_number.size() < padding_size){
            padded_number = "0" + padded_number;
        }
        outfile_path = outstem + padded_number + ".tsv";
        OutFile[i-1].open(outfile_path);
    }
    cout << assigned_cells[unassigned_barcode].size() << "\n";
    // Write unassigned read ids to file
    cout << assigned_cells[unassigned_barcode].size() << "\n";
    u_file.open(unassigned_file);
    for (string line:assigned_cells[unassigned_barcode]){
        u_file << line << "\n";
    }

    int filestate = 0;
    for (auto cell:assigned_cells){
        if ((cell.second.empty())|(cell.first == unassigned_barcode)){
            continue;
        }
        for (string line:cell.second){
            OutFile[filestate] << line.substr(0,line.size()) << "\t" << cell.first << "\n" ; 
        }
        filestate ++;
        if (filestate >= num_outfiles){
            filestate = 0;
        }
    }
}

int main(int argc, char *argv[]){
    try {
        // Define and parse command-line options
        string id_file_stem;
        string barcodes_path;
        string output_stem;
        string unassigned_file;
        int num_threads;
        int num_outfiles;
        size_t padding_size;


        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "Produce help message")
            ("barcodes,b",po::value<string>(&barcodes_path),"Path to the file containing the selected barcodes")
            ("outstem,o",po::value<string>(&output_stem),"Output file stem, final files will be {stem}XXXXX.tsv")
            ("threads,t",po::value<int>(&num_threads),"Number of threads to use")
            ("unassigned,u",po::value<string>(&unassigned_file),"File to write all unassigned read ids to")
            ("num_outfiles,n",po::value<int>(&num_outfiles),"Number of outfiles to create for downstream processing")
            ("padding",po::value<size_t>(&padding_size),"Number of digits to pad the outfile numbers to")
            ("idfile,i",po::value<string>(&id_file_stem),"Read ID file stem");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }

        po::notify(vm);

        string numstring = getenv("SLURM_ARRAY_TASK_ID");
        while (numstring.size() < padding_size){
            numstring = "0"+numstring;
        }
        
        string readnames_path = id_file_stem + numstring + ".txt";
        output_stem = output_stem + numstring + ".";
        unassigned_file = unassigned_file + numstring + ".tsv";

        if (num_outfiles > max_num_outfiles){
            cerr << "Please use a smaller number of outfiles, at most " << max_num_outfiles << "outfiles are supported!\n";
            return 1;
        }

        //create dictionaries
        string ufile_parent = unassigned_file.substr(0,unassigned_file.find_last_of("/"));
        if (!filesystem::is_directory(ufile_parent)){
            filesystem::create_directories(ufile_parent);
        }

        //get list of all allowed barcodes
        set<string> legal_barcodes;
        set<boost::dynamic_bitset<>> legal_bitcodes;
        parse_barcodes_file(legal_barcodes,legal_bitcodes,barcodes_path);
        //memory object to keep track of read assignments --> memory intensive, but best way to handle it for now
        unordered_map<string,vector<string>> read_assignments;
        assign_all_reads(read_assignments,legal_bitcodes,legal_barcodes,readnames_path,num_threads);
        //write to output files
        write_outfiles(read_assignments,output_stem,unassigned_file,num_outfiles,padding_size);
    } catch (const po::error &e) {
        cerr << "Command-line error: " << e.what() << "\n";
        return 1;
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}