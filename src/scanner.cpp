#include <cctype>
#include <fstream>
#include <iostream>

#include "../include/scanner.hpp"

namespace scanner {

FastaRecord
FastaScanner::next()
{
    FastaRecord record;

    std::string line;
    std::string sequence;

    while (std::getline(file, line)) {
        // TODO: (Abdulrasheed1729) fix this to support multi-sequence files
        if (line[0] == '>') {
            record.header = line.substr(1);
        } else {
            for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) {
                    sequence.push_back(c);
                }
            }
            record.sequence = sequence;
        }
    }

    return record;
}

bool
FastqScanner::hasNext()
{
    return file.peek() != EOF;
}

FastqRecord
FastqScanner::next()
{
    FastqRecord record;
    std::string line;

    // Line 1: header (starts with '@')
    if (!std::getline(file, line) || line.empty() || line[0] != '@') {
        throw std::runtime_error(
          "Invalid FASTQ: expected header line starting with '@'");
    }
    record.header = line.substr(1);

    // Line 2: sequence
    if (!std::getline(file, line)) {
        throw std::runtime_error("Invalid FASTQ: expected sequence line");
    }
    record.sequence = line;

    // Line 3: '+' separator
    if (!std::getline(file, line) || line.empty() || line[0] != '+') {
        throw std::runtime_error("Invalid FASTQ: expected '+' separator line");
    }

    // Line 4: quality scores
    if (!std::getline(file, line)) {
        throw std::runtime_error("Invalid FASTQ: expected quality line");
    }
    record.quality = line;

    if (record.quality.length() != record.sequence.length()) {
        throw std::runtime_error("Invalid FASTQ: quality length (" +
                                 std::to_string(record.quality.length()) +
                                 ") does not match sequence length (" +
                                 std::to_string(record.sequence.length()) +
                                 ") for record: " + record.header);
    }

    return record;
}
} // Namespace scanner
