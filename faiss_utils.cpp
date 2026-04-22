#pragma once

#include <faiss/MetricType.h>
#include <string>
#include <vector>
#include <cstdint>

#include <faiss/IndexIVFFlat.h>
#include <faiss/IndexBinaryFlat.h>
#include <faiss/IndexBinaryIVF.h>
#include <faiss/index_io.h>


#include "scanner.cpp"
#include "kmer_utils.cpp"

using KmerVector = std::vector<uint8_t>;


struct WindowMeta {
    std::string sequence_name;
    int start_pos = 0;
    // TODO(Abdulrasheed1729): add support for the reverse complement strand
    // char strand = '+'; // '+' for forward, '-' for reverse complement (unused in canonical mode)
};


std::tuple<std::vector<KmerVector>, std::vector<WindowMeta>> process_fasta_file(
    const std::string fasta_path,
    size_t window_size = 50,
    size_t stride = 1)
{
   // load fasta file and return a string of sequence

  scanner::FastaScanner scanner(fasta_path);

  scanner::FastaRecord record = scanner.next();
  std::string sequence = record.sequence;

  std::vector<KmerVector> sequences;
  std::vector<WindowMeta> window_metas;

  // make a sliding window run of each overlapping strings and encode each window as a vector of uint8_t

  for (size_t i = 0; i + window_size <= sequence.size(); i += stride) {
      std::string window = sequence.substr(i, window_size);
      // KmerVector kmer_vector(window.begin(), window.end());
      KmerVector kmer_vector = kmer_one_hot<5>(window);

      sequences.push_back(kmer_vector);
      WindowMeta meta;
      meta.sequence_name = record.header;
      meta.start_pos = i;
      window_metas.push_back(meta);
  }

  return std::make_tuple(sequences, window_metas);
}


bool create_faiss_database(
    const uint8_t* vectors,
    const std::string& index_path,
    size_t nb,
    size_t d,
    size_t nlist,
    size_t nprobe
)
{
    // static_assert(, );
    auto quantiser = std::make_unique<faiss::IndexBinaryFlat>(static_cast<faiss::idx_t>(d));
    faiss::IndexBinaryIVF index(quantiser.get(), static_cast<faiss::idx_t>(d), nlist);
    if (!index.is_trained && nb > nlist)
    {
        std::cout << "Training the Binary IVF index ...\n";
        index.train(static_cast<faiss::idx_t>(nb), vectors);
    }
    index.add(static_cast<faiss::idx_t>(nb), vectors);
    faiss::write_index_binary(&index, index_path.c_str());

    std::cout << "Index saved to " << index_path << "\n";
    return true;
}
