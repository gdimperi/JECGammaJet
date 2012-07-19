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
      EtaBin bin = {std::make_pair(0., 0.8), "eta008", "|#eta| < 0.8"};
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(0.8, 1.3);
      bin.name = "eta0813";
      bin.title = "0.8 #leq |#eta| < 1.3";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.3, 1.9);
      bin.name = "eta1319";
      bin.title = "1.3 #leq |#eta| < 1.9";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.9, 2.5);
      bin.name = "eta1925";
      bin.title = "1.9 #leq |#eta| < 2.5";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(2.5, 3.0);
      bin.name = "eta2530";
      bin.title = "2.5 #leq |#eta| < 3.0";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(3.0, 3.2);
      bin.name = "eta3032";
      bin.title = "3.0 #leq |#eta| < 3.2";

      mEtaBins.push_back(bin);
      bin.bin = std::make_pair(3.2, 5.2);
      bin.name = "eta3052";
      bin.title = "3.2 #leq |#eta| < 5.2";
      mEtaBins.push_back(bin);
    }
};
