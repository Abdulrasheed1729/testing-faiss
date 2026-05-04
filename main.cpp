#include <iostream>

#include "include/faiss_utils.hpp"
#include "include/kmer_utils.hpp"

int
main()
{
    const std::string data_path = "data/dengue_ref_sequences.fasta";
    const std::string ivf_index_path = "data/dengue_ref_sequences.ivf.index";
    const std::string flat_index_path = "data/dengue_ref_sequences.flat.index";
    const std::string dengue_left_seq = "data/left.fq";

    auto is_ivf_index_built =
      build_ivf_index(ivf_index_path, data_path, 256, 32);

    auto is_flat_index_built = build_flat_index(flat_index_path, data_path);

    if (!is_ivf_index_built) {
        std::cerr << "Error building reference database index :( \n";
    }

    if (!is_flat_index_built) {
        std::cerr << "Error building reference database index :( \n";
    }

    scanner::FastqScanner dengue_left_scanner =
      scanner::FastqScanner(dengue_left_seq);
    size_t nq = 1;
    int k = 10;
    auto d = kmer_vector_size<5>();

    auto record = dengue_left_scanner.next();

    auto ivf_indices = query_index(ivf_index_path, record, nq, k, d, false);

    auto flat_indices = query_index(flat_index_path, record, nq, k, d, true);

    compare_flat_to_ivf_index(ivf_indices, flat_indices, nq, k);

    return 0;
}
