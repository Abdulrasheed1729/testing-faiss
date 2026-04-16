#include <string>
#include <vector>

inline int base_to_bits(const char base)
{
    switch (base)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return -1;
    }
}

std::vector<bool> kmer_one_hot(const std::string& kmer, int k = 5)
{
    const int total = 1 << (2 * k);
    std::vector result(total, false);

    if (kmer.size() < k)
    {
        return result;
    }

    const int mask = total - 1;
    int index = 0;
    for (int i = 0; i < k; ++i)
    {
        index = (index << 2) | base_to_bits(kmer[i]);
    }
    result[index] = true;

    for (int i = k; i < kmer.length(); ++i)
    {
        index = ((index << 2) & mask) | base_to_bits(kmer[i]);
        result[index] = true;
    }

    return result;
}
