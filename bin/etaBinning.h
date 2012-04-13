#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <utility>

struct EtaBin {
  std::pair<float, float> bin;
  std::string name;
  std::string title;
};

class EtaBinning {
  public:
    EtaBinning() {
      fillEtaBins();
    }

    int getBin(float eta) const {
      eta = fabs(eta);
      std::vector<EtaBin>::const_iterator it = mEtaBins.begin();
      for (; it != mEtaBins.end(); ++it) {
        EtaBin bin = *it;
        if (eta >= bin.bin.first && eta < bin.bin.second) {
          return it - mEtaBins.begin();
        }
      }

      return -1;
    }

    std::string getBinName(int bin) const {
      return mEtaBins[bin].name;
    }

    std::string getBinTitle(int bin) const {
      return mEtaBins[bin].title;
    }

    size_t size() const {
      return mEtaBins.size();
    }

  private:
    std::vector<EtaBin> mEtaBins;

    void fillEtaBins() {
      EtaBin bin = {std::make_pair(0., 1.1), "eta011", "|#eta| < 1.1"};
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.5, 2.4);
      bin.name = "eta1524";
      bin.title = "1.5 < |#eta| < 2.4";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(2.5, 3.);
      bin.name = "eta243";
      bin.title = "2.4 < |#eta| < 3";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(3., 5.);
      bin.name = "eta35";
      bin.title = "3 < |#eta| < 5";
      mEtaBins.push_back(bin);
    }

};
