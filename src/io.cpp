#include <fstream>
#include <string>

#include "../include/io.hpp"

FastaRecord
read_fasta_record(const char* filepath)
{
    std::ifstream file(filepath);
    FastaRecord record;

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    std::string line;

    while (std::getline(file, line)) {
        if (line[0] == '>') {
            record.header = line.substr(1);
        } else {
            for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) {
                    record.sequence += c;
                }
            }
        }
    }

    return record;
}

