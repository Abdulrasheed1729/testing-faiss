#include "kmer_utils.cpp"
#include <gtest/gtest.h>
#include <algorithm>

// ── helpers ──────────────────────────────────────────────────────────────────

static int true_count(const std::vector<uint8_t>& v)
{
    return static_cast<int>(
        std::count_if(v.begin(), v.end(), [](uint8_t b) { return b != 0; }));
}

static int true_index(const std::vector<uint8_t>& v)
{
    for (int i = 0; i < static_cast<int>(v.size()); ++i)
        if (v[i]) return i;
    return -1;
}

// ── base_to_bits ─────────────────────────────────────────────────────────────

TEST(BaseToBitsTest, ValidBases)
{
    // Encoding: A=0, C=1, T=2, G=3
    EXPECT_EQ(base_to_bits('A'), 0);
    EXPECT_EQ(base_to_bits('C'), 1);
    EXPECT_EQ(base_to_bits('T'), 2);
    EXPECT_EQ(base_to_bits('G'), 3);
}

TEST(BaseToBitsTest, ValidLowercaseBases)
{
    // LUT now accepts lowercase — same mapping as uppercase
    EXPECT_EQ(base_to_bits('a'), 0);
    EXPECT_EQ(base_to_bits('c'), 1);
    EXPECT_EQ(base_to_bits('t'), 2);
    EXPECT_EQ(base_to_bits('g'), 3);
}

TEST(BaseToBitsTest, InvalidBasesReturnMinusOne)
{
    EXPECT_EQ(base_to_bits('N'), -1);
    EXPECT_EQ(base_to_bits('X'), -1);
    EXPECT_EQ(base_to_bits(' '), -1);
    EXPECT_EQ(base_to_bits('\0'), -1);
}

// ── kmer_one_hot: vector size ─────────────────────────────────────────────────

TEST(KmerOneHotTest, VectorSizeIsCorrectForDefaultK)
{
    // k=5 → 4^5 = 1024
    auto result = kmer_one_hot<5>("AAAAA");
    EXPECT_EQ(static_cast<int>(result.size()), 1 << (2 * 5));
}

TEST(KmerOneHotTest, VectorSizeScalesWithK)
{
    EXPECT_EQ(kmer_one_hot<1>("A").size(),      1u << 2);
    EXPECT_EQ(kmer_one_hot<2>("AA").size(),     1u << 4);
    EXPECT_EQ(kmer_one_hot<3>("AAA").size(),    1u << 6);
    EXPECT_EQ(kmer_one_hot<4>("AAAA").size(),   1u << 8);
    EXPECT_EQ(kmer_one_hot<5>("AAAAA").size(),  1u << 10);
    EXPECT_EQ(kmer_one_hot<6>("AAAAAA").size(), 1u << 12);
}

// ── kmer_one_hot: single bit set ──────────────────────────────────────────────

TEST(KmerOneHotTest, ExactlyOneBitSetForValidKmer)
{
    EXPECT_EQ(true_count(kmer_one_hot<5>("ACGTA")), 1);
    EXPECT_EQ(true_count(kmer_one_hot<5>("TTTTT")), 1);
    EXPECT_EQ(true_count(kmer_one_hot<5>("AAAAA")), 1);
    EXPECT_EQ(true_count(kmer_one_hot<5>("CATGC")), 1);
}

// ── kmer_one_hot: k=1 ────────────────────────────────────────────────────────

TEST(KmerOneHotTest, K1AllBases)
{
    // Encoding: A=0, C=1, T=2, G=3
    EXPECT_EQ(true_index(kmer_one_hot<1>("A")), 0);
    EXPECT_EQ(true_index(kmer_one_hot<1>("C")), 1);
    EXPECT_EQ(true_index(kmer_one_hot<1>("T")), 2);
    EXPECT_EQ(true_index(kmer_one_hot<1>("G")), 3);
}

// ── kmer_one_hot: k=2 ────────────────────────────────────────────────────────

TEST(KmerOneHotTest, K2KnownIndices)
{
    // AA=0, AC=1, AT=2, AG=3
    // CA=4, CC=5, CT=6, CG=7 ...
    EXPECT_EQ(true_index(kmer_one_hot<2>("AA")), 0);
    EXPECT_EQ(true_index(kmer_one_hot<2>("AC")), 1);
    EXPECT_EQ(true_index(kmer_one_hot<2>("CA")), 4);
    EXPECT_EQ(true_index(kmer_one_hot<2>("TG")), 11);  // T=2,G=3 → (2<<2)|3
    EXPECT_EQ(true_index(kmer_one_hot<2>("TT")), 10);  // T=2,T=2 → (2<<2)|2
}

// ── kmer_one_hot: k=5 known indices ──────────────────────────────────────────

TEST(KmerOneHotTest, K5AllA_IndexIsZero)
{
    EXPECT_EQ(true_index(kmer_one_hot<5>("AAAAA")), 0);
}

TEST(KmerOneHotTest, K5AllG_IndexIsLast)
{
    // G=3 → GGGGG = 0b1111111111 = 1023
    EXPECT_EQ(true_index(kmer_one_hot<5>("GGGGG")), 1023);
}

TEST(KmerOneHotTest, K5AllT_Index)
{
    // T=2 → TTTTT = (2<<8)|(2<<6)|(2<<4)|(2<<2)|2 = 682
    EXPECT_EQ(true_index(kmer_one_hot<5>("TTTTT")), 682);
}

TEST(KmerOneHotTest, K5MixedAcgta)
{
    // A→0, AC→1, ACG→7, ACGT→30, ACGTA→120
    EXPECT_EQ(true_index(kmer_one_hot<5>("ACGTA")), 120);
}

TEST(KmerOneHotTest, K5MixedCatgc)
{
    // C→1, CA→4, CAT→18, CATG→75, CATGC→301
    EXPECT_EQ(true_index(kmer_one_hot<5>("CATGC")), 301);
}

TEST(KmerOneHotTest, K5DistinctKmersProduceDifferentIndices)
{
    EXPECT_NE(true_index(kmer_one_hot<5>("ACGTA")),
              true_index(kmer_one_hot<5>("ACGTT")));
}

// ── kmer_one_hot: wrong-length input → all false ─────────────────────────────

TEST(KmerOneHotTest, TooShortReturnsAllFalse)
{
    const auto result = kmer_one_hot<5>("ACGT"); // length 4, k=5
    EXPECT_EQ(true_count(result), 0);
}

TEST(KmerOneHotTest, LongerThanK_SetsBitPerKmer)
{
    // k=2, "ACGT" has 3 k-mers: "AC"→1, "CG"→7, "GT"→14
    auto result = kmer_one_hot<2>("ACGT");
    EXPECT_EQ(true_count(result), 3);
    EXPECT_TRUE(result[1]);  // "AC": A=0,C=1 → 1
    EXPECT_TRUE(result[7]);  // "CG": C=1,G=3 → (1<<2)|3 = 7
    EXPECT_TRUE(result[14]); // "GT": G=3,T=2 → (3<<2)|2 = 14
}

TEST(KmerOneHotTest, LongerThanK_SlidingWindowIndices)
{
    // k=5, "AAAAAC": "AAAAA"→0, "AAAAC"→1
    auto result = kmer_one_hot<5>("AAAAAC");
    EXPECT_EQ(true_count(result), 2);
    EXPECT_TRUE(result[0]); // "AAAAA"
    EXPECT_TRUE(result[1]); // "AAAAC"
}

TEST(KmerOneHotTest, LongerThanK_DuplicateKmersCounted_Once)
{
    // k=2, "AAAA" → all k-mers are "AA"→0, so only 1 bit set
    auto result = kmer_one_hot<2>("AAAA");
    EXPECT_EQ(true_count(result), 1);
    EXPECT_TRUE(result[0]);
}

TEST(KmerOneHotTest, LongerThanK_BitCountEqualsDistinctKmers)
{
    // k=2, "ACGT" has 3 distinct k-mers → 3 bits set
    EXPECT_EQ(true_count(kmer_one_hot<2>("ACGT")), 3);
    // k=5, "AAAAAC" has 2 distinct k-mers → 2 bits set
    EXPECT_EQ(true_count(kmer_one_hot<5>("AAAAAC")), 2);
}

TEST(KmerOneHotTest, EmptyStringReturnsAllFalse)
{
    auto result = kmer_one_hot<5>("");
    EXPECT_EQ(true_count(result), 0);
}

// ── kmer_one_hot: invalid base throws ────────────────────────────────────────

TEST(KmerOneHotTest, InvalidBaseThrows)
{
    EXPECT_THROW(kmer_one_hot<5>("ACGNA"), std::invalid_argument);
    EXPECT_THROW(kmer_one_hot<3>("NNN"),   std::invalid_argument);
}

// ── pack_kmer_one_hot ─────────────────────────────────────────────────────────

TEST(PackKmerOneHotTest, OutputSizeIsInputDividedByEight)
{
    EXPECT_EQ(pack_kmer_one_hot(std::vector<uint8_t>(8,   0)).size(), 1u);
    EXPECT_EQ(pack_kmer_one_hot(std::vector<uint8_t>(16,  0)).size(), 2u);
    EXPECT_EQ(pack_kmer_one_hot(std::vector<uint8_t>(1024,0)).size(), 128u);
}

TEST(PackKmerOneHotTest, EmptyInputProducesEmptyOutput)
{
    EXPECT_TRUE(pack_kmer_one_hot({}).empty());
}

TEST(PackKmerOneHotTest, InputSmallerThanEightProducesEmptyOutput)
{
    // 7 bits → 7>>3 = 0 bytes
    EXPECT_TRUE(pack_kmer_one_hot(std::vector<uint8_t>(7, 1)).empty());
}

TEST(PackKmerOneHotTest, TrailingBitsIgnored)
{
    // 9 elements → only first 8 packed, 9th dropped
    std::vector<uint8_t> v(9, 0);
    v[8] = 1; // only the trailing bit is set
    auto packed = pack_kmer_one_hot(v);
    ASSERT_EQ(packed.size(), 1u);
    EXPECT_EQ(packed[0], 0u); // trailing bit was not packed
}

TEST(PackKmerOneHotTest, MSBFirst_FirstBitSet)
{
    // bit 0 → most significant bit of byte 0
    std::vector<uint8_t> v(8, 0);
    v[0] = 1;
    EXPECT_EQ(pack_kmer_one_hot(v)[0], 0b10000000u); // 128
}

TEST(PackKmerOneHotTest, MSBFirst_LastBitSet)
{
    // bit 7 → least significant bit of byte 0
    std::vector<uint8_t> v(8, 0);
    v[7] = 1;
    EXPECT_EQ(pack_kmer_one_hot(v)[0], 0b00000001u); // 1
}

TEST(PackKmerOneHotTest, AllBitsSet)
{
    EXPECT_EQ(pack_kmer_one_hot(std::vector<uint8_t>(8, 1))[0], 0xFFu);
}

TEST(PackKmerOneHotTest, NoBitsSet)
{
    EXPECT_EQ(pack_kmer_one_hot(std::vector<uint8_t>(8, 0))[0], 0u);
}

TEST(PackKmerOneHotTest, AlternatingBits)
{
    // [1,0,1,0,1,0,1,0] → 0b10101010 = 170
    std::vector<uint8_t> v = {1,0,1,0,1,0,1,0};
    EXPECT_EQ(pack_kmer_one_hot(v)[0], 0b10101010u);
}

TEST(PackKmerOneHotTest, NonZeroValuesTreatedAsSet)
{
    // values 5 and 255 are both non-zero → bits 0 and 1 set → 0b11000000 = 192
    std::vector<uint8_t> v = {5, 255, 0, 0, 0, 0, 0, 0};
    EXPECT_EQ(pack_kmer_one_hot(v)[0], 0b11000000u);
}

TEST(PackKmerOneHotTest, MultiByteKnownValues)
{
    // byte 0: [1,0,0,0,0,0,0,0] → 128
    // byte 1: [0,0,0,0,0,0,0,1] → 1
    std::vector<uint8_t> v(16, 0);
    v[0]  = 1;
    v[15] = 1;
    auto packed = pack_kmer_one_hot(v);
    ASSERT_EQ(packed.size(), 2u);
    EXPECT_EQ(packed[0], 128u);
    EXPECT_EQ(packed[1], 1u);
}

TEST(PackKmerOneHotTest, IntegrationWithKmerOneHot_AAAAA)
{
    // kmer_one_hot<5>("AAAAA") sets only index 0
    // → byte 0 has bit 0 (MSB) set = 128, all other bytes = 0
    auto vec    = kmer_one_hot<5>("AAAAA");
    auto packed = pack_kmer_one_hot(vec);
    ASSERT_EQ(packed.size(), 128u);
    EXPECT_EQ(packed[0], 128u);
    for (size_t i = 1; i < packed.size(); ++i)
        EXPECT_EQ(packed[i], 0u) << "byte " << i << " should be 0";
}

TEST(PackKmerOneHotTest, IntegrationWithKmerOneHot_GGGGG)
{
    // kmer_one_hot<5>("GGGGG") sets only index 1023
    // byte_idx = 1023/8 = 127, bit pos = 7 → shift 0 → packed[127] = 1
    auto vec    = kmer_one_hot<5>("GGGGG");
    auto packed = pack_kmer_one_hot(vec);
    ASSERT_EQ(packed.size(), 128u);
    for (size_t i = 0; i < 127; ++i)
        EXPECT_EQ(packed[i], 0u) << "byte " << i << " should be 0";
    EXPECT_EQ(packed[127], 1u);
}