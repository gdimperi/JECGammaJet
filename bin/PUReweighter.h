#include <string>
#include <TH1.h>
#include <map>

enum class PUProfile : uint8_t {
  S6,
  S7
};

class PUReweighter {
  public:
    PUReweighter(const std::string& dataFile, const std::string& mcFile);
    PUReweighter(const std::string& dataFile, PUProfile profile = PUProfile::S7);

    ~PUReweighter() {
      delete puHisto;
    }

    double weight(float interactions) const;

  private:
    void initPUProfiles();

    TH1* puHisto;

    std::map<PUProfile, std::vector<double>> mPUCoefs;
};

