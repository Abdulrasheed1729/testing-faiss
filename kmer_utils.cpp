#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

template <size_t KMER>
static constexpr size_t kmer_vector_size()
{
    return (1u << (2u * KMER));
}

// Branchless lookup table — faster than switch, GPU-friendly (no branches)
// Maps ASCII -> 2-bit encoding: A=0, C=1, T=2, G=3, invalid=-1
static constexpr int8_t BASE_LUT[256] = {
    // Generated: only A/a=0, C/c=1, T/t=2, G/g=3 are valid
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 0-15
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 16-31
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 32-47
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 48-63
    -1, 0,-1, 1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1, // 64-79  (A,C,G)
    -1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 80-95  (T)
    -1, 0,-1, 1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1, // 96-111 (a,c,g)
    -1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 112-127 (t)
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 128-143
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 144-159
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 160-175
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 176-191
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 192-207
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 208-223
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 224-239
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 240-255
};

inline int8_t base_to_bits(char base)
{
    return BASE_LUT[static_cast<uint8_t>(base)];
}

template <size_t KMER>
bool kmer_one_hot_raw(
    const char*  sequence,
    size_t       seq_len,
    uint8_t*     out,
    int*         error_flag
)
{
    constexpr size_t N    = kmer_vector_size<KMER>();
    constexpr size_t mask = N - 1;

    std::memset(out, 0, N);

    if (seq_len < KMER)
        return false; // invalid sequence length

    // Seed the first k-mer
    size_t index = 0;
    for (size_t i = 0; i < KMER; ++i)
    {
        int8_t code = base_to_bits(sequence[i]);
        if (code == -1)
        {
            *error_flag = 1;
            return false;
        }
        index = (index << 2) | static_cast<size_t>(code);
    }
    out[index] = 1;

    // Slide the window
    for (size_t i = KMER; i < seq_len; ++i)
    {
        int8_t code = base_to_bits(sequence[i]);
        if (code == -1)
        {
            *error_flag = 1;
            return false;
        }
        index = ((index << 2) & mask) | static_cast<size_t>(code);
        out[index] = 1;
    }

    return true;
}

template <size_t KMER>
void kmer_one_hot_batch(
    const char* const* sequences,  // array of C-strings
    const size_t*      seq_lens,   // length of each sequence
    size_t             num_seqs,
    uint8_t*           out,        // pre-allocated: num_seqs * kmer_vector_size<KMER>()
    int*               error_flags // pre-allocated: num_seqs, zero-initialised by caller
)
{
    constexpr size_t N = kmer_vector_size<KMER>();
    for (size_t s = 0; s < num_seqs; ++s)
    {
        kmer_one_hot_raw<KMER>(
            sequences[s],
            seq_lens[s],
            out + s * N,
            &error_flags[s]
        );
    }
}

template <size_t KMER>
std::vector<uint8_t> kmer_one_hot(const std::string& sequence)
{
    constexpr size_t N = kmer_vector_size<KMER>();
    std::vector<uint8_t> result(N, 0);

    int error_flag = 0;
    kmer_one_hot_raw<KMER>(
        sequence.data(),
        sequence.size(),
        result.data(),
        &error_flag
    );

    if (error_flag)
        throw std::invalid_argument("Invalid base in sequence");

    return result;
}
