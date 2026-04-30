#include <assert.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <faiss/IndexBinaryFlat.h>
#include <faiss/IndexBinaryIVF.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/MetricType.h>
#include <faiss/index_io.h>

#include "../include/faiss_utils.hpp"
#include "../include/kmer_utils.hpp"

std::tuple<std::vector<KmerVector>, std::vector<WindowMetaData>>
process_fasta_file(const std::string fasta_path,
                   size_t window_size,
                   size_t stride)
{
    // load fasta file and return a string of sequence

    scanner::FastaScanner scanner(fasta_path);

    // NOTE: for now we are only assuming that we have a single sequence fasta
    // file
    scanner::FastaRecord record = scanner.next();
    std::string sequence = record.sequence;

    std::vector<KmerVector> sequences;
    std::vector<WindowMetaData> window_metas;

    // make a sliding window run of each overlapping strings and encode each
    // window as a vector of uint8_t

    for (size_t i = 0; i + window_size <= sequence.size(); i += stride) {
        std::string window = sequence.substr(i, window_size);
        // KmerVector kmer_vector(window.begin(), window.end());
        KmerVector kmer_vector = kmer_one_hot<KMER_K>(window);

        sequences.push_back(kmer_vector);
        WindowMetaData meta;
        // FIX: there should not be a copying operation herre
        meta.sequence_name = record.header;
        meta.start_pos = i;
        window_metas.push_back(meta);
    }

    return std::make_tuple(sequences, window_metas);
}

bool
build_index(const std::string& index_path,
            const std::string& data_path, // size_t d,
            size_t nlist,
            size_t nprobe)
{

    // assert(d % 8 == 0); // d must be greater than 0

    auto d_bits = kmer_vector_size<5>();

    auto [vectors, windows] = process_fasta_file(data_path);

    assert(!vectors.empty());

    auto quantiser = std::make_unique<faiss::IndexBinaryFlat>(
      static_cast<faiss::idx_t>(d_bits));

    faiss::IndexBinaryIVF index(
      quantiser.get(), static_cast<faiss::idx_t>(d_bits), nlist);

    std::vector<uint8_t> packed_vectors;
    std::vector<uint8_t> packed_buffer;
    auto nb = vectors.size();

    if (nb <= 0)
        return false;

    auto d_bytes = d_bits >> 3;
    packed_vectors.resize(nb * d_bytes);

    for (size_t i = 0; i < vectors.size(); ++i) {
        auto current_ptr = packed_vectors.data() + i * d_bytes;
        pack_kmer_one_hot(vectors[i], packed_buffer);
        std::memcpy(current_ptr, packed_buffer.data(), packed_buffer.size());
    }

    if (!index.is_trained && nb > nlist) {
        std::cout << "Training the Binary IVF index ...\n";
        index.train(static_cast<faiss::idx_t>(nb), packed_vectors.data());
    }

    index.add(static_cast<faiss::idx_t>(nb), packed_vectors.data());
    faiss::write_index_binary(&index, index_path.c_str());

    std::cout << "Index saved to " << index_path << "\n";

    std::cout << "Saving metadata ...\n";

    std::ofstream meta_file(index_path + ".meta.tsv");
    for (size_t i = 0; i < windows.size(); ++i) {
        meta_file << i << "\t" << windows[i].sequence_name << "\t"
                  << windows[i].start_pos << "\n";
    }

    std::cout << "Metadata saved to " << index_path << ".meta.tsv\n";

    return true;
}

template<typename T>
std::string
vec_to_string(const std::vector<T>& vec)
{
    std::ostringstream oss;
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(oss, ", "));

    return oss.str();
}

void
query_index(const std::string& index_path,
            const scanner::FastqRecord& query_sequence,
            size_t nq,
            int k,
            size_t d)
// size_t nprobe)
{
    // load index
    auto loaded_index =
      faiss::read_index_binary_up(index_path.c_str(), faiss::IO_FLAG_READ_ONLY);

    // load metadata
    std::unordered_map<size_t, WindowMetaData> meta_map;
    std::ifstream meta_file(index_path + ".meta.tsv");
    std::string line;
    while (std::getline(meta_file, line)) {
        size_t idx;
        std::string idx_str;
        std::string start_pos_str;
        std::string seq_name;
        WindowMetaData meta;
        std::istringstream iss(line);
        std::getline(iss, idx_str, '\t');
        std::getline(iss, seq_name, '\t');
        std::getline(iss, start_pos_str, '\t');
        idx = std::stoul(idx_str.c_str());
        meta.start_pos = std::stoi(start_pos_str);
        meta.sequence_name = seq_name;
        // iss >> idx >> meta.sequence_name >> meta.start_pos;
        meta_map[idx] = meta;
        // std::cout << meta.start_pos << "\n";
        // std::cout << meta.sequence_name << "\n";
    }

    KmerVector query_vec = kmer_one_hot<KMER_K>(query_sequence.sequence);
    std::vector<uint8_t> packed_query(nq * d / 8);
    pack_kmer_one_hot(query_vec, packed_query);

    std::cout << "Total Index: " << loaded_index->ntotal << "\n";
    std::cout << "Total d: " << loaded_index->d << "\n";

    std::vector<faiss::idx_t> indices(nq * k);
    std::vector<int32_t> distances(nq * k);
    // loaded_index->nprobe = nprobe; // set nprobe
    auto t_start = std::chrono::steady_clock::now();
    loaded_index->search(
      nq, packed_query.data(), k, distances.data(), indices.data());

    auto t_end = std::chrono::steady_clock::now();
    double elapsed =
      std::chrono::duration<double, std::milli>(t_end - t_start).count();
    std::cout << "Query time: " << elapsed << "ms.\n";
    for (int i = 0; i < k; ++i) {
        auto meta = meta_map[indices[i]];
        std::cout << meta.sequence_name << ":" << meta.start_pos
                  << " (distance: " << distances[i] << ") "
                  << query_sequence.header << "\n";
    }
}
