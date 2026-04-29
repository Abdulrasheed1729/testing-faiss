#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "scanner.hpp"

using KmerVector = std::vector<uint8_t>;

constexpr size_t KMER_K = 5;

struct WindowMetaData
{
    std::string sequence_name;
    int start_pos = 0;
    // NOTE: should this POD struct change to something of this form:
    // struct WindowMeta
    // {
    //     std::string sequence_name;
    //     std::vector<int> indexes;
    // };
    // ?
};

bool
create_faiss_database(const uint8_t* vectors,
                      const std::string& index_path,
                      size_t nb,
                      size_t d,
                      size_t nlist,
                      size_t nprobe);

std::tuple<std::vector<KmerVector>, std::vector<WindowMetaData>>
process_fasta_file(const std::string fasta_path,
                   size_t window_size = 50,
                   size_t stride = 1);

bool
build_index(const std::string& index_path,
            const std::string& data_path,
            // size_t d,
            size_t nlist,
            size_t nprobe);

std::vector<KmerVector>
load_index(const std::string& index_path);

void
query_index(const std::string& index_path,
            const scanner::FastqRecord& query_sequence,
            size_t nq,
            int k,
            size_t d);
// size_t nprobe);
