/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iterator>
#include <sstream>

#include "faiss_utils.cpp"
#include "scanner.cpp"

// #include <faiss/IndexFlat.h>
// #include <faiss/IndexIVFPQ.h>
// #include <faiss/index_io.h>
#include <vector>

template<typename T>
std::string
vec_to_string(const std::vector<T>& vec)
{
    std::ostringstream oss;
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(oss, ", "));

    return oss.str();
}

// double elapsed()
// {
//     struct timeval tv;
//     gettimeofday(&tv, nullptr);
//     return tv.tv_sec + tv.tv_usec * 1e-6;
// }

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

    // scanner::FastqScanner fqscanner("data/left.fq");
    // std::ofstream outfile("data/left.kmer.bin", std::ios::binary);

    // while (fqscanner.hasNext())
    // {
    //     scanner::FastqRecord record = fqscanner.next();
    //     std::vector<bool> kmer_vec = kmer_one_hot(record.sequence);
    //     std::string kmer_str = vec_to_string(kmer_vec);
    //     std::cout << kmer_str << std::endl;
    //     outfile.write(kmer_str.c_str(), kmer_str.size());
    //     outfile.write("\n", 1);
    //     // std::cout << record.header << std::endl;
    //     std::cout << record.sequence << std::endl;
    //     std::cout << sizeof(kmer_vec[0]) << std::endl;
    //     std::cout << kmer_vec[0] << std::endl;
    //     std::cout << "Size of bool: " << sizeof(true) << std::endl;
    //     std::cout << typeid(kmer_vec[0]).name() << ", size=" <<
    //     sizeof(kmer_vec[0]) << "\n"; std::cout << typeid(true).name() << ",
    //     size=" << sizeof(true) << "\n";
    //     // std::cout << record.quality << std::endl;
    // }
}

// int main() {
//     double t0 = elapsed();

//     // dimension of the vectors to index
//     int d = 128;

//     // size of the database we plan to index
//     size_t nb = 200 * 1000;

//     // make a set of nt training vectors in the unit cube
//     // (could be the database)
//     size_t nt = 100 * 1000;

//     // make the index object and train it
//     faiss::IndexFlatL2 coarse_quantizer(d);

//     // a reasonable number of centroids to index nb vectors
//     int ncentroids = int(4 * sqrt(nb));

//     // the coarse quantizer should not be deallocated before the index
//     // 4 = nb of bytes per code (d must be a multiple of this)
//     // 8 = nb of bits per sub-code (almost always 8)
//     faiss::IndexIVFPQ index(&coarse_quantizer, d, ncentroids, 4, 8);

//     std::mt19937 rng;

//     { // training
//         printf("[%.3f s] Generating %ld vectors in %dD for training\n",
//                elapsed() - t0,
//                nt,
//                d);

//         std::vector<float> trainvecs(nt * d);
//         std::uniform_real_distribution<> distrib;
//         for (size_t i = 0; i < nt * d; i++) {
//             trainvecs[i] = distrib(rng);
//         }

//         printf("[%.3f s] Training the index\n", elapsed() - t0);
//         index.verbose = true;

//         index.train(nt, trainvecs.data());
//     }

//     { // I/O demo
//         const char* outfilename = "/tmp/index_trained.faissindex";
//         printf("[%.3f s] storing the pre-trained index to %s\n",
//                elapsed() - t0,
//                outfilename);

//         write_index(&index, outfilename);
//     }

//     size_t nq;
//     std::vector<float> queries;

//     { // populating the database
//         printf("[%.3f s] Building a dataset of %ld vectors to index\n",
//                elapsed() - t0,
//                nb);

//         std::vector<float> database(nb * d);
//         std::uniform_real_distribution<> distrib;
//         for (size_t i = 0; i < nb * d; i++) {
//             database[i] = distrib(rng);
//         }

//         printf("[%.3f s] Adding the vectors to the index\n", elapsed() - t0);

//         index.add(nb, database.data());

//         printf("[%.3f s] imbalance factor: %g\n",
//                elapsed() - t0,
//                index.invlists->imbalance_factor());

//         // remember a few elements from the database as queries
//         int i0 = 1234;
//         int i1 = 1243;

//         nq = i1 - i0;
//         queries.resize(nq * d);
//         for (int i = i0; i < i1; i++) {
//             for (int j = 0; j < d; j++) {
//                 queries[(i - i0) * d + j] = database[i * d + j];
//             }
//         }
//     }

//     { // searching the database
//         int k = 5;
//         printf("[%.3f s] Searching the %d nearest neighbors "
//                "of %ld vectors in the index\n",
//                elapsed() - t0,
//                k,
//                nq);

//         std::vector<faiss::idx_t> nns(k * nq);
//         std::vector<float> dis(k * nq);

//         index.search(nq, queries.data(), k, dis.data(), nns.data());

//         printf("[%.3f s] Query results (vector ids, then distances):\n",
//                elapsed() - t0);

//         for (int i = 0; i < nq; i++) {
//             printf("query %2d: ", i);
//             for (int j = 0; j < k; j++) {
//                 printf("%7ld ", nns[j + i * k]);
//             }
//             printf("\n     dis: ");
//             for (int j = 0; j < k; j++) {
//                 printf("%7g ", dis[j + i * k]);
//             }
//             printf("\n");
//         }

//         printf("note that the nearest neighbor is not at "
//                "distance 0 due to quantization errors\n");
//     }

//     return 0;
// }
