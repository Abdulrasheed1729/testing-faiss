#include <iostream>

#include "include/faiss_utils.hpp"
#include "include/kmer_utils.hpp"
// #include "include/scanner.hpp"

int
main()
{
    const std::string data_path = "data/dengue_ref_sequences.fasta";
    const std::string index_path = "data/dengue_ref_sequences.index";
    const std::string dengue_left_seq = "data/left.fq";

    auto is_index_built = build_index(index_path, data_path, 128, 10);

    if (!is_index_built) {
        std::cerr << "Error building reference database index :( \n";
    }

    scanner::FastqScanner dengue_left_scanner =
      scanner::FastqScanner(dengue_left_seq);
    size_t nq = 1;
    int k = 10;
    auto d = kmer_vector_size<5>();

    // while (dengue_left_scanner.hasNext()) {
    auto record = dengue_left_scanner.next();

    // auto t_start = std::chrono::steady_clock::now();
    query_index(index_path, record, nq, k, d);
    // auto t_end = std::chrono::steady_clock::now();

    // double elapsed =
    // std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // std::cout << "Query time: " << elapsed << "ms.\n";
    // }

    return 0;
}
