#pragma once

#include <cmath>
#include <vector>
#include <utility>

class PtBinning {
  public:
    PtBinning() {
      fillPtBins();
    }

    int getPtBin(float pt) {
      std::vector<std::pair<float, float> >::const_iterator it = mPtBins.begin();
      for (; it != mPtBins.end(); ++it) {
        std::pair<float, float> bin = *it;
        if (pt >= bin.first && pt < bin.second) {
          return it - mPtBins.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mPtBins.size();
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mPtBins[bin];
    }

    std::vector<std::pair<float, float> > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<std::pair<float, float> >(mPtBins.begin(), mPtBins.begin() + n);
    }

    std::vector<std::pair<float, float> > getBinning(unsigned int from, unsigned int to) const {
      if (to > size()) {
        to = size();
      }

      return std::vector<std::pair<float, float> >(mPtBins.begin() + from, mPtBins.begin() + to);
    }

  private:
    std::vector<std::pair<float, float> > mPtBins;

    void fillPtBins() {
      mPtBins.push_back(std::make_pair(40., 50.));
      mPtBins.push_back(std::make_pair(50., 60.));
      mPtBins.push_back(std::make_pair(60., 75.));
      mPtBins.push_back(std::make_pair(75., 95.));
      mPtBins.push_back(std::make_pair(95., 125.));
      mPtBins.push_back(std::make_pair(125., 155.));
      mPtBins.push_back(std::make_pair(155., 180.));
      mPtBins.push_back(std::make_pair(180., 210.));
      mPtBins.push_back(std::make_pair(210., 250.));
      mPtBins.push_back(std::make_pair(250., 300.));
      mPtBins.push_back(std::make_pair(300., 350.));
      mPtBins.push_back(std::make_pair(350., 400.));
      mPtBins.push_back(std::make_pair(400., 500.));
      mPtBins.push_back(std::make_pair(500., 600.));
      mPtBins.push_back(std::make_pair(600., 800.));
    }

};
