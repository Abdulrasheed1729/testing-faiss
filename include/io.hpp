#pragma once

#include <string>

/* 
  A simple POD for storing a FASTA record
  	- header: the header line
  	- sequence: the sequence
 * */
struct FastaRecord
{
    std::string header;
    std::string sequence;
};

FastaRecord
read_fasta_record(const char* filepath);

// TODO: add the support for the fastq file


