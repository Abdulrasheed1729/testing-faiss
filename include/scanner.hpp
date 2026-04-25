#include <fstream>
#include <string>

namespace scanner {

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

class FastaScanner
{

  public:
    explicit FastaScanner(const std::string filename)
      : filename(std::move(filename)) {};
    [[nodiscard]] FastaRecord next() const;

  private:
    std::string filename;
};

class FastqScanner
{
  public:
    explicit FastqScanner(const std::string& filename)
      : filename(filename)
      , file(filename) {
		if (!file.is_open()) {
			throw std::runtime_error("Failed to open file");
		}
	};
    ~FastqScanner() { file.close(); };
    bool hasNext();
    FastqRecord next();

  private:
    std::string filename;
    std::ifstream file;
};

} // namespace scanner
