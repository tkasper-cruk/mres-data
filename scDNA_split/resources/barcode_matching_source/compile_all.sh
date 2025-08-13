#!/bin/bash
g++ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -O3 -msse2 -Wall -pedantic -ftree-vectorizer-verbose=1 -lboost_program_options -o ../count_matches count_matches.cpp
g++ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -O3 -msse2 -Wall -pedantic -ftree-vectorizer-verbose=1 -lz -lboost_program_options -pthread -lboost_thread -o ../assign_barcodes assign_barcodes.cpp
g++ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -O3 -msse2 -Wall -pedantic -ftree-vectorizer-verbose=1 -lboost_iostreams -lboost_program_options -o ../split_fastq split_fastq.cpp