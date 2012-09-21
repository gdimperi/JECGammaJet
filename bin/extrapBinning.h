#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <utility>

#include "ptBinning.h"

class ExtrapBinning {
  public:
    ExtrapBinning() {}

    void initialize(PtBinning ptBinning, const std::string& recoType) {
      size_t s = ptBinning.size();
      for (size_t i = 0; i < s; i++) {

        std::pair<float, float> bin = ptBinning.getBinValue(i);
        float minPt = 0.;
        float ptStep = 0.;
        // In percent
        if (recoType == "PFlow" || recoType == "JPT") {
          if (bin.first <= 80.) {
            minPt = 8.;
            ptStep = 3.;
          } else if (bin.first <= 350) {
            minPt = 5.;
            ptStep = 2.5;
          } else {
            minPt = 2.;
            ptStep = 2.;
          }
        } else {
          if (bin.first <= 60.) {
            minPt = 8.;
            ptStep = 1.5;
          } else if (bin.first <= 350.) {
            minPt = 6.;
            ptStep = 1.5;
          } else {
            minPt = 2.;
            ptStep = 2.;
          }
        }
        //      if (bin.first <= 40.)
        //        minPt += 2 * ptStep;
        //      if (bin.first <= 30.)
        //        minPt += ptStep;
        //
        minPt /= 100.;
        ptStep /= 100.;
        float maxPt = minPt + ptStep;

        mMapping.push_back(std::make_pair(minPt, maxPt));
      }

    }

    int getBin(float ptPhoton, float ptSecondJet, int ptBin) const {

      std::pair<float, float> mapping = mMapping[ptBin]; // mPtBinning.getBinValue(ptBin);

      double ptRel = ptSecondJet / ptPhoton;
      size_t extrapBin = (size_t) floor((ptRel - mapping.first) / (mapping.second - mapping.first));

      return (extrapBin >= size()) ? -1 : extrapBin;
    }

    size_t size() const {
      return mSize;
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mMapping[bin];
    }

  private:
    PtBinning mPtBinning;
    static const int mSize = 10;

    std::vector<std::pair<float, float> > mMapping; // first is minPt, second is maxPt
};
