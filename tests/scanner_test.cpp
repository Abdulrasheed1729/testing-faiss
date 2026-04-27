#include <gtest/gtest.h>

#include "../include/scanner.hpp"

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "data"
#endif

static const std::string DATA_DIR = TEST_DATA_DIR;

using namespace scanner;

TEST(FastqScannerTest, ValidFastqFile)
{
    FastqScanner scanner(DATA_DIR + "/test.fq");

    ASSERT_TRUE(scanner.hasNext());

    FastqRecord record = scanner.next();

    ASSERT_EQ(record.header, "simulated.1/1");
    ASSERT_EQ(
      record.sequence,
      "ACGAGTTGCTGATCGAATTCCACACACACCCTCCTCCCATGCCTTCCCAATGGCTGCTGATAGTCT"
      "CTTTGGGGAATCAGCTTGAAATTTGTATTGCTCT");
    ASSERT_EQ(record.quality,
              "IHIIIIHIIIGHHHGHIHHHIEIHIFIIIIHGIIDEIGIIHIIHGIIIFIHF=FI<HHHIGIF?"
              "GAIIIIBIEIIHEIIHIIFIIBBIIII;I3IIHIIF");

    ASSERT_FALSE(scanner.hasNext());
}

TEST(FastqScannerTest, MultipleRecords)
{
    FastqScanner scanner(DATA_DIR + "/test_multi.fq");

    ASSERT_TRUE(scanner.hasNext());
    FastqRecord first = scanner.next();
    ASSERT_EQ(first.header, "simulated.1/1");
    ASSERT_EQ(first.sequence, "ACGAGTTGCTGATCGAATTCCACACACACC");
    ASSERT_EQ(first.quality, "IHIIIIHIIIGHHHGHIHHHIEIHIFIIII");

    ASSERT_TRUE(scanner.hasNext());
    FastqRecord second = scanner.next();
    ASSERT_EQ(second.header, "simulated.2/1");
    ASSERT_EQ(second.sequence, "ATCCAGATGTGCCGGGTTTAAAATCTAGGGC");
    ASSERT_EQ(second.quality, "HHHHIIIHHHIHIHHHHIHHIGIEBGGGEIF");

    ASSERT_FALSE(scanner.hasNext());
}

TEST(FastqScannerTest, EmptyFile)
{
    FastqScanner scanner(DATA_DIR + "/test_empty.fq");
    ASSERT_FALSE(scanner.hasNext());
}

TEST(FastqScannerTest, FileNotFound)
{
    ASSERT_THROW(FastqScanner(DATA_DIR + "/nonexistent.fq"),
                 std::runtime_error);
}

TEST(FastqScannerTest, HasNextIsIdempotent)
{
    FastqScanner scanner(DATA_DIR + "/test.fq");

    // Calling hasNext() repeatedly must not consume any record
    ASSERT_TRUE(scanner.hasNext());
    ASSERT_TRUE(scanner.hasNext());
    ASSERT_TRUE(scanner.hasNext());

    scanner.next();

    ASSERT_FALSE(scanner.hasNext());
    ASSERT_FALSE(scanner.hasNext());
}

TEST(FastqScannerTest, InvalidHeader)
{
    FastqScanner scanner(DATA_DIR + "/test_invalid_header.fq");
    ASSERT_THROW(scanner.next(), std::runtime_error);
}

TEST(FastqScannerTest, InvalidSeparator)
{
    FastqScanner scanner(DATA_DIR + "/test_invalid_separator.fq");
    ASSERT_THROW(scanner.next(), std::runtime_error);
}

TEST(FastqScannerTest, QualityLengthMismatch)
{
    FastqScanner scanner(DATA_DIR + "/test_quality_mismatch.fq");
    ASSERT_THROW(scanner.next(), std::runtime_error);
}

TEST(FastqScannerTest, LargeFile)
{
    FastqScanner scanner(DATA_DIR + "/left.fq");

    std::size_t count = 0;
    while (scanner.hasNext()) {
        FastqRecord record = scanner.next();
        ASSERT_FALSE(record.header.empty());
        ASSERT_EQ(record.sequence.length(), record.quality.length());
        ++count;
    }

    ASSERT_EQ(count, 25000u);
}

// ── FastaScanner tests
// ────────────────────────────────────────────────────────

TEST(FastaScannerTest, SingleRecord)
{
    FastaScanner scanner(DATA_DIR + "/test_single.fa");
    FastaRecord record = scanner.next();

    EXPECT_EQ(record.header, "seq1 simple test sequence");
    EXPECT_EQ(record.sequence, "ACGTACGTTTAAGGCC");
}

TEST(FastaScannerTest, MultiLineSequenceConcatenated)
{
    // Two sequence lines must be joined into one contiguous string.
    FastaScanner scanner(DATA_DIR + "/test_single.fa");
    FastaRecord record = scanner.next();

    // "ACGTACGT" + "TTAAGGCC" with no separator
    EXPECT_EQ(record.sequence.size(), 16u);
    EXPECT_EQ(record.sequence, "ACGTACGTTTAAGGCC");
}

// TEST(FastaScannerTest, FileNotFound)
// {
//     FastaScanner scanner(DATA_DIR + "/nonexistent.fa");
//     EXPECT_THROW((void)scanner.next(), std::runtime_error);
// }

TEST(FastaScannerTest, EmptyFile)
{
    // An empty file should return an empty record without throwing.
    FastaScanner scanner(DATA_DIR + "/test_empty.fa");
    FastaRecord record = scanner.next();

    EXPECT_TRUE(record.header.empty());
    EXPECT_TRUE(record.sequence.empty());
}

TEST(FastaScannerTest, HeaderOnlyNoSequence)
{
    // A record that has a header line but no sequence lines.
    FastaScanner scanner(DATA_DIR + "/test_header_only.fa");
    FastaRecord record = scanner.next();

    EXPECT_EQ(record.header, "seq1 header only no sequence");
    EXPECT_TRUE(record.sequence.empty());
}

// Note:** `MultiRecordUsesLastHeaderAndAllSequences` deliberately pins
// the *current* (arguably quirky) behavior — if the parser is ever
// refactored to return one record per `>` header, that test will
// fail and alert you to update the call sites.
TEST(FastaScannerTest, MultiRecordUsesLastHeaderAndAllSequences)
{
    // Current FastaScanner::next() reads the whole file in one shot:
    // the last '>' line wins as the header, and every non-header line
    // is concatenated into a single sequence string.
    FastaScanner scanner(DATA_DIR + "/test_multi.fa");
    FastaRecord record = scanner.next();

    EXPECT_EQ(record.header, "seq2 second sequence");
    EXPECT_EQ(record.sequence, "AAAACCCCGGGGTTTT");
}

TEST(FastaScannerTest, WhitespaceStrippedFromSequence)
{
    // Sequence must contain no whitespace characters.
    FastaScanner scanner(DATA_DIR + "/test_single.fa");
    FastaRecord record = scanner.next();

    for (char c : record.sequence) {
        EXPECT_FALSE(std::isspace(static_cast<unsigned char>(c)))
          << "Unexpected whitespace character: " << static_cast<int>(c);
    }
}

TEST(FastaScannerTest, LargeFastaFile)
{
    FastaScanner scanner(DATA_DIR + "/dengue_ref_sequences.fasta");
    FastaRecord record = scanner.next();

    EXPECT_EQ(record.header, "NC_001477.1 |Dengue virus 1, complete genome");
    EXPECT_FALSE(record.sequence.empty());

    // Dengue genome is ~10 700 bp; sequence must contain only ACGT.
    EXPECT_GT(record.sequence.size(), 10000u);
    for (char c : record.sequence) {
        EXPECT_NE(std::string("ACGT").find(c), std::string::npos)
          << "Unexpected base character: " << c;
    }
}
