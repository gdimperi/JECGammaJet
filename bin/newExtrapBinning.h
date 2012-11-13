#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <utility>

class NewExtrapBinning {
  public:
    NewExtrapBinning() {}

    void initialize(float alpha) {
      mAlpha = alpha;
      mSize = alpha / 0.05;

      // Construct binning
      for (int i = 0; i < mSize; i++) {
        std::pair<float, float> bin = std::make_pair(i * 0.05, (i + 1) * 0.05);
        mBinning.push_back(bin);
      }
    }

    int getBin(float ptPhoton, float ptSecondJet) const {

      float alpha = ptSecondJet / ptPhoton;

      std::vector<std::pair<float, float>>::const_iterator it = mBinning.begin();
      for (; it != mBinning.end(); ++it) {
        std::pair<float, float> bin = *it;
        if (alpha >= bin.first && alpha < bin.second) {
          return it - mBinning.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mSize;
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mBinning[bin];
    }

    float getBinWidth() const {
      return 0.05;
    }

  private:
    float mAlpha;
    int mSize;

    std::vector<std::pair<float, float>> mBinning;
};
