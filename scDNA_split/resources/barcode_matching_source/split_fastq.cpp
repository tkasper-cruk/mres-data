#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace po = boost::program_options;
using namespace std;

#define barcode_length 16

void parse_barcodes_file(unordered_map<string,string> &read_to_barcode,unordered_map<string,vector<string>> &readbuffer,string assignment_file){
    string line;
    string read_id;
    string barcode;
    ifstream assignment_stream;
    size_t split_pos;
    assignment_stream.open(assignment_file);
    while (getline(assignment_stream,line)) {
        split_pos = line.find('\t');
        read_id = line.substr(0,split_pos);
        barcode = line.substr(split_pos+1,barcode_length);
        read_to_barcode[read_id] = barcode;
        if (readbuffer.count(barcode) == 0){
            readbuffer[barcode]= vector<string>();
        }
    }
}

void process_read(
    array<string,4> read,
    size_t &current_buffer_size,
    unordered_map<string,string> &read_to_barcode,
    unordered_map<string,vector<string>> &readbuffer
    // ofstream &u_file
){
    string read_id = read[0].substr(0,read[0].find(' '));
    if (read_to_barcode.count(read_id)==0){ //unassigned
        // u_file << read[0] << "\n" << read[1] << "\n" << read[2] << "\n" << read[3] << "\n";
        return;
    }
    string barcode = read_to_barcode[read_id];
    for (int i=0;i<4;i++){
        readbuffer[barcode].push_back(read[i]);
    }
    current_buffer_size ++;
}

string get_outfile_name(
    string barcode,
    string slx,
    string sample,
    string flowcell,
    string lane,
    string split_dir,
    string read_direction
){
    return split_dir + slx + "." + sample + "-" + barcode + "." + flowcell + ".s_" + lane + ".r_" + read_direction + ".fq";
}

void dump_buffer(
    unordered_map<string,vector<string>> &readbuffer,
    string slx,
    string sample,
    string flowcell,
    string lane,
    string split_dir,
    string read_direction
) {
    ofstream writer;
    string filename;
    for (auto& cell:readbuffer){
        if (cell.second.empty()){
            continue;
        }
        
        filename = get_outfile_name(cell.first,slx,sample,flowcell,lane,split_dir,read_direction);
        writer.open(filename,ios_base::app);
        for (const string &line:cell.second){
            writer << line << "\n";
        }
        writer.close();
        cell.second.clear();
       // readbuffer[cell.first] = vector<string>();
    }
}



void split_fastq(
    unordered_map<string,string> &read_to_barcode,
    unordered_map<string,vector<string>> &readbuffer,
    string split_dir,
    string fastq_file,
    string slx,
    string sample,
    string flowcell,
    string lane,
    string read_direction,
    size_t buffer_size
){
    size_t buffered_reads = 0;
    // ofstream u_stream;
    //check that buffer is actually empty
    for (auto cell:readbuffer){
        buffered_reads += cell.second.size();
    }
    assert(buffered_reads == 0);
    // u_stream.open(u_filename);

    ifstream fqfile(fastq_file, ios_base::in | ios_base::binary);
    if (!fqfile.is_open()) {
        cerr << "Error: Unable to open file " << fastq_file << "\n";
        // u_stream.close();
        return;
    }
    boost::iostreams::filtering_stream<boost::iostreams::input> fq_in;
    fq_in.push(boost::iostreams::gzip_decompressor());
    fq_in.push(fqfile);

    string line;
    array<string,4> read;
    int state = 0;
    while (getline(fq_in,line)){
        if ((state == 0)& (line[0]=='@')) {//header
            read[0] = line;
            state ++;
        } else if (state == 0) {//wrong lines
            cout << line << "\n";
        } else {//complete read in memory
            read[state] = line;
            state ++;
        }   
        if (state >= 4){
            process_read(read,buffered_reads,read_to_barcode,readbuffer);
            state = 0;
            read.fill("");
            if (buffered_reads >= buffer_size){
                dump_buffer(readbuffer,slx,sample,flowcell,lane,split_dir,read_direction);
                buffered_reads = 0;
            }
        }
    }
    // u_stream.close();
    if (buffered_reads > 0){
        dump_buffer(readbuffer,slx,sample,flowcell,lane,split_dir,read_direction);
        buffered_reads = 0;
    }
}

int main(int argc, char *argv[]) {
    try {
        // Define and parse command-line options
        string r1_file;
        string r2_file;
        string output_dir;
        string assignment_stem;
        string slx;
        string sample;
        string flowcell;
        string lane;
        size_t max_buffer_reads;
        size_t padding_size;


        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "Produce help message")
            ("r1_file,1", po::value<string>(&r1_file), "Gzip-compressed R1 FastQ")
            ("r2_file,2", po::value<string>(&r2_file), "Gzip-compressed R2 FastQ")
            ("output_dir,o", po::value<string>(&output_dir), "File to write the split FastQ files to")
            ("assigments,a", po::value<string>(&assignment_stem), "List of read-id to cell assignments")
            ("slx", po::value<string>(&slx), "SLX ID")
            ("sample", po::value<string>(&sample), "Sample ID")
            ("flowcell", po::value<string>(&flowcell), "Flowcell ID")
            ("lane", po::value<string>(&lane), "Lane number")
            ("padding,p",po::value<size_t>(&padding_size),"Number of digits to pad numbers to")
            ("max_buffer_reads", po::value<size_t>(&max_buffer_reads), "Maximum number of reads to keep in memory");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }

        po::notify(vm);

        if (output_dir[output_dir.size()-1] != '/'){ //remove / from directory paths
            output_dir = output_dir+ "/";
        }
        if (!filesystem::is_directory(output_dir)){
            filesystem::create_directories(output_dir);
        }

        string numstring = getenv("SLURM_ARRAY_TASK_ID");
        while (numstring.size() < padding_size){
            numstring = "0"+numstring;
        }

        string read_assignment = assignment_stem + numstring + ".tsv";

        //parse assignment file
        unordered_map<string,string> read_to_barcode;
        unordered_map<string,vector<string>> buffered_reads;
        parse_barcodes_file(read_to_barcode,buffered_reads,read_assignment);

        //split fastq
        split_fastq(
            read_to_barcode,
            buffered_reads,
            output_dir,
            r1_file,
            slx,
            sample,
            flowcell,
            lane,
            "1",
            max_buffer_reads
        );
        
        split_fastq(
            read_to_barcode,
            buffered_reads,
            output_dir,
            r2_file,
            slx,
            sample,
            flowcell,
            lane,
            "2",
            max_buffer_reads
        );
    } catch (const po::error &e) {
        cerr << "Command-line error: " << e.what() << "\n";
        return 1;
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
