#include <string>
#include <vector>

inline int base_to_bits(char base) {
  switch (base) {
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

std::vector<bool> kmer_one_hot(const std::string &kmer, int k = 5) {
  int total = 1 << (2 * k);
  std::vector<bool> result(total, false);

  if (kmer.length() != k) {
    return result;
  }

  int mask = total - 1;
  int index = 0;
  for (char base : kmer) {
    index = (index << 2) | base_to_bits(base);
  }
  result[index & mask] = true;

  return result;
}
