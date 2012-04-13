// derived class from DrawBase
// draws the response and resolution extrapolations
// to pt (2nd Jet) -> 0

#pragma once

#include "drawBase.h"
#include "TFile.h"

#include "ptBinning.h"
#include "etaBinning.h"
#include "extrapBinning.h"

class drawExtrap : public drawBase {

public:

  drawExtrap(const std::string& analysisType, const std::string& recoType = "", const std::string& jetAlgo = "", bool outputGraphs = false, const std::string& flags = "");
  virtual ~drawExtrap() {};

  void set_FITRMS(const std::string& fit_rms) {
    FIT_RMS_ = fit_rms;
  };
  void set_NOQ(bool noq) {
    NOQ_ = noq;
  };
  void set_INTPERC(float intperc) {
    INTPERC_ = intperc;
  };
  void set_FIXM(bool fixm) {
    FIXM_ = fixm;
  };
  void set_EXCLUDEFIRSTPOINT(bool exclfirstpoint) {
    EXCLUDE_FIRST_POINT_ = exclfirstpoint;
  };

  void drawResponseExtrap(const std::string& etaRegion, bool rawJets);


private:

  void drawExtrapSinglePtBin(int iPtBin, Float_t& DATAReso, Float_t& DATAReso_err, Float_t& extrapReso, Float_t& extrapReso_err, Float_t& trueReso, Float_t& trueReso_err, Float_t& intrReso, Float_t& intrReso_err,
                             Float_t& DATAResp, Float_t& DATAResp_err, Float_t& extrapResp, Float_t& extrapResp_err, Float_t& trueResp, Float_t& trueResp_err, Float_t& intrResp, Float_t& intrResp_err,
                             Float_t& imbalanceResp, Float_t& imbalanceResp_err, const std::string& recoGen) const;
  void getXPoints(int ptBin, Float_t* x, Float_t* x_err) const;
  void getYPoints(TFile * file, const char* yHistoName, Int_t nPoints, Float_t* y_resp, Float_t* y_resp_err,  Float_t* y_reso, Float_t* y_reso_err) const;
  void getYPointsVector(TFile *file, const std::string& yHistoName, Int_t nPoints, std::vector<Float_t>& y_resp, std::vector<Float_t>& y_resp_err, std::vector<Float_t>& y_reso, std::vector<Float_t>& y_reso_err) const;

  bool NOQ_;
  std::string FIT_RMS_;
  float INTPERC_;
  bool FIXM_;
  bool EXCLUDE_FIRST_POINT_;

  PtBinning mPtBinning;
  EtaBinning mEtaBinning;
  ExtrapBinning mExtrapBinning;
};
