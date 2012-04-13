#ifndef fitTools_h
#define fitTools_h


#include "TFile.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TMatrixD.h"


class fitTools {

public:

  static double delta_phi(double phi1, double phi2);

  static float delta_phi(float phi1, float phi2);

  static std::vector<float> getPtPhot_binning();

  static void getBins(int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog = true);

  static void getBins_int(int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog = true);

  static void getPtBins(int nBins, Double_t* Lower, bool plotLog = true);

  static int getNbins_stack(const std::string& varName);

  static void getBins_stack(int nBins, Double_t* Lower, const std::string& varName);

  static void drawSingleGraph(TGraph* gr, const std::string& canvasName);

  static void fitProjection(TH1* h1_projection, TF1* gaussian, Float_t nSigma = 1.5, std::string option = "RQ", bool add = false);

  static void fitProjection_sameArea(TH1* h1_projection, TF1* gaussian, TH1** newhisto, Float_t percIntegral = 0.9, const std::string& option = "RQ", bool useMode = false);

  static void getTruncatedMeanAndRMS(TH1* h1_projection, Float_t& mean, Float_t& mean_err, Float_t& rms, Float_t& rms_err, Double_t percentIntegral_MEAN = 0.9, Double_t percentIntegral_RMS = 0.68);
// static TCanvas* getTruncatedMeanAndRMS(TH1D* h1_projection, Float_t& mean, Float_t& mean_err, Float_t& rms, Float_t& rms_err, Double_t percentIntegral_MEAN=0.9, Double_t percentIntegral_RMS=0.68) {

  static void fillProfile(TH1F* h1_response_FIT, TH1F* h1_resolution_FIT, TH1F* h1_response_MEAN, TH1F* h1_resolution_RMS, TH2D* h2, std::string name = "");

  //new fit reponse:
  static void fitDistribution_TGraph(TH2D* h2, TH2D* genMean, const std::string& varY, const std::string& varX, const std::string& etaRegion, const std::string& flag, const std::string& algoType, const std::string& outFileName, const std::string& name = "", Float_t percIntegral = 0.95, bool use_samearea = false);

  static void fillPositionResolution(TH1F* h1_sigmaEta, TH1F* h1_sigmaPhi, TH2D* h2_deltaEta, TH2D* h2_deltaPhi);

  //used by getEfficiencyHisto (later on):
  static int getEfficiencyUncertainties(int n, int k, double p, double& xmin, double& xmax);

  static TGraphAsymmErrors* getEfficiencyGraph(const std::string& name, TH1F* h1_numerator, TH1F* h1_denominator);

  static TF1* fitResolutionGraph(TGraphErrors* graph , std::string funcType, std::string funcName, const std::string& option = "RQ", float rangeMax = 1000., float rangeMin = 10.);

  static TF1* fitResponseGraph(TGraphErrors* graph , std::string funcType, std::string funcName, const std::string& option = "RQ", float rangeMax = 350., float rangeMin = 10.);

  static TH1D* getBand(TF1* f, const std::string& name);

  static TH1D* getBand(TF1* f, TMatrixD const& m, std::string name, bool getRelativeBand = false, int npx = 100);

  static TGraphErrors* get_graphRatio(TGraphErrors* gr_data, TGraphErrors* gr_MC);

  //static TGraphAsymmErrors* getGraphPoissonErrors( TH1D* histo, const std::string xerrType="binWidth", float cl=0.6826 );
  static TGraphAsymmErrors* getGraphPoissonErrors(TH1* histo, const std::string xerrType = "binWidth", float nSigma = 1.);



private:

};


#endif
