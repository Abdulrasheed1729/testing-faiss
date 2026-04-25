#pragma once

#include <string>
#include <cstdint>
#include <vector>
#include <tuple>

using KmerVector = std::vector<uint8_t>;

struct WindowMeta
{
    std::string sequence_name;
    int start_pos = 0;
    // TODO(Abdulrasheed1729): add support for the reverse complement strand
    // char strand = '+'; // '+' for forward, '-' for reverse complement (unused
    // in canonical mode)
};


bool
create_faiss_database(const uint8_t* vectors,
                      const std::string& index_path,
                      size_t nb,
                      size_t d,
                      size_t nlist,
                      size_t nprobe);


std::tuple<std::vector<KmerVector>, std::vector<WindowMeta>>
process_fasta_file(const std::string fasta_path,
                   size_t window_size = 50,
                   size_t stride = 1);


