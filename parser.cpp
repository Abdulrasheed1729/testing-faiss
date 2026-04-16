#include <fstream>
#include <iostream>
#include <utility>

namespace parser
{
    struct FastaRecord
    {
        std::string header;
        std::string sequence;
    };

    struct FastqRecord
    {
        std::string header;
        std::string sequence;
        std::string quality;
    };

    class FastaParser
    {
    private:
        std::string filename;

    public:
        explicit FastaParser(std::string filename) : filename(std::move(filename))
        {
        }

        [[nodiscard]] FastaRecord next() const
        {
            std::ifstream file(filename);
            FastaRecord record;

            if (!file.is_open())
            {
                throw std::runtime_error("Failed to open file");
            }

            if (file.is_open())
            {
                std::string line;
                while (std::getline(file, line))
                {
                    if (line[0] == '>')
                    {
                        record.header = line.substr(1);
                    }
                    else
                    {
                        record.sequence += line;
                    }
                    if (line[0] == '>')
                    {
                        break;
                    }
                }
            }

            return record;
        }
    };

    class FastqParser
    {
    private:
        std::string filename;
        std::ifstream file;

    public:
        explicit FastqParser(const std::string& filename) : filename(filename), file(filename)
        {
            if (!file.is_open())
            {
                throw std::runtime_error("Failed to open file: " + filename);
            }
        }

        bool hasNext() { return file.peek() != EOF; }

        FastqRecord next()
        {
            FastqRecord record;
            std::string line;

            // Line 1: header (starts with '@')
            if (!std::getline(file, line) || line.empty() || line[0] != '@')
            {
                throw std::runtime_error(
                    "Invalid FASTQ: expected header line starting with '@'");
            }
            record.header = line.substr(1);

            // Line 2: sequence
            if (!std::getline(file, line))
            {
                throw std::runtime_error("Invalid FASTQ: expected sequence line");
            }
            record.sequence = line;

            // Line 3: '+' separator
            if (!std::getline(file, line) || line.empty() || line[0] != '+')
            {
                throw std::runtime_error("Invalid FASTQ: expected '+' separator line");
            }

            // Line 4: quality scores
            if (!std::getline(file, line))
            {
                throw std::runtime_error("Invalid FASTQ: expected quality line");
            }
            record.quality = line;

            if (record.quality.length() != record.sequence.length())
            {
                throw std::runtime_error("Invalid FASTQ: quality length (" +
                    std::to_string(record.quality.length()) +
                    ") does not match sequence length (" +
                    std::to_string(record.sequence.length()) +
                    ") for record: " + record.header);
            }

            return record;
        }
    };
} // Namespace parser
