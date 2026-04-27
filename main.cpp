#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

#include "include/faiss_utils.hpp"
#include "include/kmer_utils.hpp"
#include "include/scanner.hpp"

template<typename T>
std::string
vec_to_string(const std::vector<T>& vec)
{
    std::ostringstream oss;
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(oss, ", "));

    return oss.str();
}

int
main()
{
    scanner::FastaScanner fasta_scanner("data/dengue_ref_sequences.fasta");
    scanner::FastaRecord record = fasta_scanner.next();
    std::cout << "Header: " << record.header << std::endl;
    std::cout << "Sequence: " << record.sequence << std::endl;
    auto [vectors, windows] =
      process_fasta_file("data/dengue_ref_sequences.fasta");
    std::ofstream outfile("data/dengue_ref_sequences.kmer.bin");
    std::vector<uint8_t> packed_vectors;
    auto nb = vectors.size();
    auto d_bits = kmer_vector_size<5>();

    auto d_bytes = d_bits >> 3;
    packed_vectors.resize(nb * d_bytes);

    for (size_t i = 0; i < vectors.size(); ++i) {
        auto current_ptr = packed_vectors.data() + i * d_bytes;
        auto dd = pack_kmer_one_hot(vectors[i]);
        std::memcpy(current_ptr, dd.data(), dd.size());
    }

    auto db = create_faiss_database(packed_vectors.data(),
                                    "data/dengue_ref_sequences.kmer.index",
                                    vectors.size(),
                                    d_bits,
                                    128,
                                    8);
}
