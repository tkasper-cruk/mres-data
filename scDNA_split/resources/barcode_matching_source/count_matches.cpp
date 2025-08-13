#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

#define barcode_length 16

void parse_barcodes_file(set<string> &legal_barcodes,string barcodes_path){
    ifstream infile(barcodes_path);
    string header;
    
    getline(infile,header);
    bool use_first = (header.find("DNA") < header.find("\t"));
    
    string line;
    string barcode;
    while (getline(infile,line)){
        if (use_first) {
            barcode = line.substr(0,barcode_length);
        } else {
            barcode = line.substr(line.find('\t')+1,barcode_length);
        }
        legal_barcodes.insert(barcode);
    }

}

void count_reads(unordered_map<string,int> &readcounts, set<string> legal_barcodes, string readpath) {
    //initialise counts
    for (string bcd : legal_barcodes) {
        readcounts[bcd] = 0;
    }

    ifstream readfile(readpath);
    string line;
    string barcode;
    while (getline(readfile,line)){
        barcode = line.substr(line.find(' ')-barcode_length,barcode_length);
        if (legal_barcodes.find(barcode)!=legal_barcodes.end()) {
            readcounts[barcode]++;
        }
    }
}

void write_outfile(unordered_map<string,int> readcounts,string outpath){
    ofstream outfile(outpath);
    outfile << "Barcode\tcount\n";
    for (auto cell:readcounts) {
        if (cell.second > 0){
            outfile << cell.first << "\t" << cell.second << "\n";
        }
    }
}

int main(int argc, char *argv[]){
    try {
        // Define and parse command-line options
        string readnames_stem;
        string barcodes_path;
        string output_stem;
        size_t padding_size;

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "Produce help message")
            ("readnames,i",po::value<string>(&readnames_stem),"Stem for the txt file containing the read ids (1 per line, uncompressed)")
            ("barcodes,b",po::value<string>(&barcodes_path),"Path to a file containing the line-matched DNA and RNA barcodes, with the DNA and RNA headers")
            ("outstem,o", po::value<string>(&output_stem), "Output file stem")
            ("padding,p",po::value<size_t>(&padding_size),"Number of digits to pad numbers to");

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
        
        string readnames_path = readnames_stem + numstring + ".txt";
        string output_path = output_stem + numstring + ".tsv";

        //get list of all allowed barcodes
        set<string> legal_barcodes;
        parse_barcodes_file(legal_barcodes,barcodes_path);
        //set up counter
        unordered_map<string,int> readcounts;
        count_reads(readcounts,legal_barcodes,readnames_path);
        //write counts to file
        write_outfile(readcounts,output_path);
    
    } catch (const po::error &e) {
        cerr << "Command-line error: " << e.what() << "\n";
        return 1;
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}