#include "scanner.cpp"
#include <gtest/gtest.h>

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "data"
#endif

static const std::string DATA_DIR = TEST_DATA_DIR;

using namespace scanner;

TEST(FastqParserTest, ValidFastqFile)
{
    FastqScanner parser(DATA_DIR + "/test.fq");

    ASSERT_TRUE(parser.hasNext());

    FastqRecord record = parser.next();

    ASSERT_EQ(record.header, "simulated.1/1");
    ASSERT_EQ(record.sequence,
              "ACGAGTTGCTGATCGAATTCCACACACACCCTCCTCCCATGCCTTCCCAATGGCTGCTGATAGTCTCTTTGGGGAATCAGCTTGAAATTTGTATTGCTCT");
    ASSERT_EQ(record.quality,
              "IHIIIIHIIIGHHHGHIHHHIEIHIFIIIIHGIIDEIGIIHIIHGIIIFIHF=FI<HHHIGIF?GAIIIIBIEIIHEIIHIIFIIBBIIII;I3IIHIIF");

    ASSERT_FALSE(parser.hasNext());
}

TEST(FastqParserTest, MultipleRecords)
{
    FastqScanner parser(DATA_DIR + "/test_multi.fq");

    ASSERT_TRUE(parser.hasNext());
    FastqRecord first = parser.next();
    ASSERT_EQ(first.header, "simulated.1/1");
    ASSERT_EQ(first.sequence, "ACGAGTTGCTGATCGAATTCCACACACACC");
    ASSERT_EQ(first.quality, "IHIIIIHIIIGHHHGHIHHHIEIHIFIIII");

    ASSERT_TRUE(parser.hasNext());
    FastqRecord second = parser.next();
    ASSERT_EQ(second.header, "simulated.2/1");
    ASSERT_EQ(second.sequence, "ATCCAGATGTGCCGGGTTTAAAATCTAGGGC");
    ASSERT_EQ(second.quality, "HHHHIIIHHHIHIHHHHIHHIGIEBGGGEIF");

    ASSERT_FALSE(parser.hasNext());
}

TEST(FastqParserTest, EmptyFile)
{
    FastqScanner parser(DATA_DIR + "/test_empty.fq");
    ASSERT_FALSE(parser.hasNext());
}

TEST(FastqParserTest, FileNotFound)
{
    ASSERT_THROW(FastqScanner(DATA_DIR + "/nonexistent.fq"), std::runtime_error);
}

TEST(FastqParserTest, HasNextIsIdempotent)
{
    FastqScanner parser(DATA_DIR + "/test.fq");

    // Calling hasNext() repeatedly must not consume any record
    ASSERT_TRUE(parser.hasNext());
    ASSERT_TRUE(parser.hasNext());
    ASSERT_TRUE(parser.hasNext());

    parser.next();

    ASSERT_FALSE(parser.hasNext());
    ASSERT_FALSE(parser.hasNext());
}

TEST(FastqParserTest, InvalidHeader)
{
    FastqScanner parser(DATA_DIR + "/test_invalid_header.fq");
    ASSERT_THROW(parser.next(), std::runtime_error);
}

TEST(FastqParserTest, InvalidSeparator)
{
    FastqScanner parser(DATA_DIR + "/test_invalid_separator.fq");
    ASSERT_THROW(parser.next(), std::runtime_error);
}

TEST(FastqParserTest, QualityLengthMismatch)
{
    FastqScanner parser(DATA_DIR + "/test_quality_mismatch.fq");
    ASSERT_THROW(parser.next(), std::runtime_error);
}

TEST(FastqParserTest, LargeFile)
{
    FastqScanner parser(DATA_DIR + "/left.fq");

    std::size_t count = 0;
    while (parser.hasNext())
    {
        FastqRecord record = parser.next();
        ASSERT_FALSE(record.header.empty());
        ASSERT_EQ(record.sequence.length(), record.quality.length());
        ++count;
    }

    ASSERT_EQ(count, 25000u);
}
