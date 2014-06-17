#pragma once

#include <cmath>
#include <vector>
#include <utility>

class VertexBinning {
  public:
    VertexBinning() {
      fillVertexBins();
    }

    int getVertexBin(int n) {
      std::vector<std::pair<int, int> >::const_iterator it = mVertexBins.begin();
      for (; it != mVertexBins.end(); ++it) {
        std::pair<int, int> bin = *it;
        if (n >= bin.first && n < bin.second) {
          return it - mVertexBins.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mVertexBins.size();
    }

    std::pair<int, int> getBinValue(int bin) const {
      return mVertexBins[bin];
    }

    std::vector<std::pair<int, int> > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<std::pair<int, int> >(mVertexBins.begin(), mVertexBins.begin() + n);
    }

    std::vector<std::pair<int, int> > getBinning(unsigned int from, unsigned int to) const {
      if (to > size()) {
        to = size();
      }

      return std::vector<std::pair<int, int> >(mVertexBins.begin() + from, mVertexBins.begin() + to);
    }

  private:
    std::vector<std::pair<int, int> > mVertexBins;

    void fillVertexBins() {
      mVertexBins.push_back(std::make_pair(0, 5));
      mVertexBins.push_back(std::make_pair(5, 8));
      mVertexBins.push_back(std::make_pair(8, 11));
      mVertexBins.push_back(std::make_pair(11, 13));
      mVertexBins.push_back(std::make_pair(13, 15));
      mVertexBins.push_back(std::make_pair(15, 18));
      mVertexBins.push_back(std::make_pair(18, 21));
      mVertexBins.push_back(std::make_pair(21, 23));
      mVertexBins.push_back(std::make_pair(23, 35));
    }

};
