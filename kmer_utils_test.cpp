#include "kmer_utils.cpp"
#include <gtest/gtest.h>
#include <algorithm>

// ── helpers ──────────────────────────────────────────────────────────────────

static int true_count(const std::vector<bool>& v)
{
    return std::count(v.begin(), v.end(), true);
}

static int true_index(const std::vector<bool>& v)
{
    for (int i = 0; i < static_cast<int>(v.size()); ++i)
        if (v[i]) return i;
    return -1;
}

// ── base_to_bits ─────────────────────────────────────────────────────────────

TEST(BaseToBitsTest, ValidBases)
{
    EXPECT_EQ(base_to_bits('A'), 0);
    EXPECT_EQ(base_to_bits('C'), 1);
    EXPECT_EQ(base_to_bits('G'), 2);
    EXPECT_EQ(base_to_bits('T'), 3);
}

TEST(BaseToBitsTest, InvalidBasesReturnMinusOne)
{
    EXPECT_EQ(base_to_bits('N'), -1);
    EXPECT_EQ(base_to_bits('X'), -1);
    EXPECT_EQ(base_to_bits('a'), -1); // lowercase not accepted
    EXPECT_EQ(base_to_bits('t'), -1);
    EXPECT_EQ(base_to_bits(' '), -1);
}

// ── kmer_one_hot: vector size ─────────────────────────────────────────────────

TEST(KmerOneHotTest, VectorSizeIsCorrectForDefaultK)
{
    // default k=5 → 4^5 = 1024
    auto result = kmer_one_hot("AAAAA");
    EXPECT_EQ(static_cast<int>(result.size()), 1 << (2 * 5));
}

TEST(KmerOneHotTest, VectorSizeScalesWithK)
{
    for (int k = 1; k <= 6; ++k)
    {
        std::string kmer(k, 'A');
        auto result = kmer_one_hot(kmer, k);
        EXPECT_EQ(static_cast<int>(result.size()), 1 << (2 * k)) << "k=" << k;
    }
}

// ── kmer_one_hot: single bit set ──────────────────────────────────────────────

TEST(KmerOneHotTest, ExactlyOneBitSetForValidKmer)
{
    EXPECT_EQ(true_count(kmer_one_hot("ACGTA")), 1);
    EXPECT_EQ(true_count(kmer_one_hot("TTTTT")), 1);
    EXPECT_EQ(true_count(kmer_one_hot("AAAAA")), 1);
    EXPECT_EQ(true_count(kmer_one_hot("CATGC")), 1);
}

// ── kmer_one_hot: k=1 ────────────────────────────────────────────────────────

TEST(KmerOneHotTest, K1AllBases)
{
    EXPECT_EQ(true_index(kmer_one_hot("A", 1)), 0);
    EXPECT_EQ(true_index(kmer_one_hot("C", 1)), 1);
    EXPECT_EQ(true_index(kmer_one_hot("G", 1)), 2);
    EXPECT_EQ(true_index(kmer_one_hot("T", 1)), 3);
}

// ── kmer_one_hot: k=2 ────────────────────────────────────────────────────────

TEST(KmerOneHotTest, K2KnownIndices)
{
    // AA=0, AC=1, AG=2, AT=3
    // CA=4, CC=5, CG=6, CT=7 ...
    EXPECT_EQ(true_index(kmer_one_hot("AA", 2)), 0);
    EXPECT_EQ(true_index(kmer_one_hot("AC", 2)), 1);
    EXPECT_EQ(true_index(kmer_one_hot("CA", 2)), 4);
    EXPECT_EQ(true_index(kmer_one_hot("TG", 2)), 14);
    EXPECT_EQ(true_index(kmer_one_hot("TT", 2)), 15);
}

// ── kmer_one_hot: k=5 known indices ──────────────────────────────────────────

TEST(KmerOneHotTest, K5AllA_IndexIsZero)
{
    EXPECT_EQ(true_index(kmer_one_hot("AAAAA")), 0);
}

TEST(KmerOneHotTest, K5AllT_IndexIsLast)
{
    // 4^5 - 1 = 1023
    EXPECT_EQ(true_index(kmer_one_hot("TTTTT")), 1023);
}

TEST(KmerOneHotTest, K5MixedAcgta)
{
    // A→0, AC→1, ACG→6, ACGT→27, ACGTA→108
    EXPECT_EQ(true_index(kmer_one_hot("ACGTA")), 108);
}

TEST(KmerOneHotTest, K5MixedCatgc)
{
    // C→1, CA→4, CAT→19, CATG→78, CATGC→313
    EXPECT_EQ(true_index(kmer_one_hot("CATGC")), 313);
}

TEST(KmerOneHotTest, K5DistinctKmersProduceDifferentIndices)
{
    EXPECT_NE(true_index(kmer_one_hot("ACGTA")),
              true_index(kmer_one_hot("ACGTT")));
}

// ── kmer_one_hot: wrong-length input → all false ─────────────────────────────

TEST(KmerOneHotTest, TooShortReturnsAllFalse)
{
    const auto result = kmer_one_hot("ACGT"); // length 4, k=5
    EXPECT_EQ(true_count(result), 0);
}

TEST(KmerOneHotTest, LongerThanK_SetsBitPerKmer)
{
    // k=2, "ACGT" has 3 k-mers: "AC"→1, "CG"→6, "GT"→11
    auto result = kmer_one_hot("ACGT", 2);
    EXPECT_EQ(true_count(result), 3);
    EXPECT_TRUE(result[1]); // "AC"
    EXPECT_TRUE(result[6]); // "CG"
    EXPECT_TRUE(result[11]); // "GT"
}

TEST(KmerOneHotTest, LongerThanK_SlidingWindowIndices)
{
    // k=5, "AAAAAC": "AAAAA"→0, "AAAAC"→1
    auto result = kmer_one_hot("AAAAAC", 5);
    EXPECT_EQ(true_count(result), 2);
    EXPECT_TRUE(result[0]); // "AAAAA"
    EXPECT_TRUE(result[1]); // "AAAAC"
}

TEST(KmerOneHotTest, LongerThanK_DuplicateKmersCounted_Once)
{
    // k=2, "AAAA" → all k-mers are "AA"→0, so only 1 bit set
    auto result = kmer_one_hot("AAAA", 2);
    EXPECT_EQ(true_count(result), 1);
    EXPECT_TRUE(result[0]);
}

TEST(KmerOneHotTest, LongerThanK_BitCountEqualsDistinctKmers)
{
    // k=2, "ACGT" has 3 distinct k-mers → 3 bits set
    EXPECT_EQ(true_count(kmer_one_hot("ACGT", 2)), 3);
    // k=5, "AAAAAC" has 2 distinct k-mers → 2 bits set
    EXPECT_EQ(true_count(kmer_one_hot("AAAAAC", 5)), 2);
}

TEST(KmerOneHotTest, EmptyStringReturnsAllFalse)
{
    auto result = kmer_one_hot("");
    EXPECT_EQ(true_count(result), 0);
}

// ── kmer_one_hot: invalid base ────────────────────────────────────────────────
// base_to_bits returns -1 for non-ACGT characters. The first k-mer's index
// computation has no mask applied, so an invalid base produces an out-of-bounds
// vector access (undefined behaviour). No tests are provided for this case;
// callers are responsible for ensuring inputs contain only A, C, G, T.
