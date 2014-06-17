#include "fitTools.h"
#include <cmath>
#include "TMinuit.h"
#include "RooHistError.h"


Double_t rpf(Double_t* x, Double_t* p);
Double_t NSC(Double_t* x, Double_t* p);
Double_t NSCPF(Double_t* x, Double_t* p);
Double_t powerlaw(Double_t* x, Double_t* p);


double fitTools::delta_phi(double phi1, double phi2) {

  double dphi = fabs(phi1 - phi2);
  return (dphi <= TMath::Pi()) ? dphi : TMath::TwoPi() - dphi;
}


float fitTools::delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}



std::vector<float> fitTools::getPtPhot_binning() {

  std::vector<float> returnVector;

////returnVector.push_back(10.);
//  returnVector.push_back(15.);
//  returnVector.push_back(22.);
//  returnVector.push_back(32.);
//  returnVector.push_back(47.);
//  returnVector.push_back(70.);
//  returnVector.push_back(100.);
//  returnVector.push_back(150.);
//  //returnVector.push_back(220.);
//  returnVector.push_back(320.);
//  returnVector.push_back(470.);
//  returnVector.push_back(3500.);

  returnVector.push_back(15.);
  returnVector.push_back(22.);
  returnVector.push_back(32.);
  returnVector.push_back(53.);
  returnVector.push_back(80.);
  returnVector.push_back(100.);
  returnVector.push_back(150.);
  returnVector.push_back(220.);
  returnVector.push_back(320.);
  returnVector.push_back(470.);
  returnVector.push_back(700.);
  returnVector.push_back(3500.);

//returnVector.push_back(15.);
//returnVector.push_back(18.);
//returnVector.push_back(22.);
//returnVector.push_back(27.);
//returnVector.push_back(32.);
//returnVector.push_back(39.);
//returnVector.push_back(47.);
//returnVector.push_back(57.);
//returnVector.push_back(69.);
//returnVector.push_back(84.);
//returnVector.push_back(100.);
//returnVector.push_back(3500.);

  return returnVector;

}

void fitTools::getBins(int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog) {

  int nBins = nBins_total - 1;
  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
  Lower[0] = xmin;
  for (int i = 1; i != nBins; ++i) {

    if (plotLog) {
      Lower[i] = Lower[i - 1] * dx;
    } else {
      Lower[i] = Lower[i - 1] + dx;
    }


  }

  Lower[nBins] = xmax;

}

void fitTools::getBins_int(int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog) {

  Double_t Lower_exact;
  int nBins = nBins_total - 1;
  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
  Lower[0] = xmin;
  Lower_exact = Lower[0];
  for (int i = 1; i != nBins; ++i) {

    if (plotLog) {
      Lower_exact *= dx;
      Lower[i] = TMath::Ceil(Lower_exact);
    } else {
      Lower[i] = TMath::Ceil(Lower[i - 1] + dx);
    }

  }

  Lower[nBins] = xmax;

}


void fitTools::getPtBins(int nBins, Double_t* Lower, bool plotLog) {

  //hardwired for now:
//Int_t ii=0;
//Lower[ii++] = 20.;
////Lower[ii++] = 26.;
//Lower[ii++] = 34.;
//Lower[ii++] = 44.;
//Lower[ii++] = 57.;
//Lower[ii++] = 74.;
//Lower[ii++] = 96.;
//Lower[ii++] = 124.;
//Lower[ii++] = 181.;
//Lower[ii++] = 230.;
//Lower[ii++] = 312.;
//Lower[ii++] = 457.;
//Lower[ii++] = 700.;
//Lower[ii++] = 902.;
//Lower[ii++] = 1176.;
//Lower[ii++] = 1533.;
//Lower[ii++] = 2000.;

  Int_t ii = 0;
  Lower[ii++] = 20.;
  Lower[ii++] = 42.;
  Lower[ii++] = 62.;
  Lower[ii++] = 87.;
  Lower[ii++] = 108.;
  Lower[ii++] = 140.;
  Lower[ii++] = 183.;
  Lower[ii++] = 239.;
  Lower[ii++] = 312.;
  Lower[ii++] = 407.;
  Lower[ii++] = 530.;
  Lower[ii++] = 691.;
  Lower[ii++] = 902.;
  Lower[ii++] = 1176.;
  Lower[ii++] = 1533.;
  Lower[ii++] = 2000.;

}


int fitTools::getNbins_stack(const std::string& varName) {

  int nBins = 0;

  if (varName == "eta") {
    nBins = 51;
  } else if (varName == "pt" || varName == "ptCorr") {
    nBins = 51;
  } else if (varName == "phi") {
    nBins = 51;
  } else {
    std::cout << "Binning not yet implememented for variable '" << varName << "'. Exiting." << std::endl;
    exit(38112);
  }

  return nBins;

}

void fitTools::getBins_stack(int nBins, Double_t* Lower, const std::string& varName) {

  if (varName == "eta") {
    getBins(nBins, Lower, -5., 5., (bool)false);
  } else if (varName == "pt" || varName == "ptCorr") {
    getBins(nBins, Lower, 5., 100., (bool)true);
  } else if (varName == "phi") {
    getBins(nBins, Lower, -3.142, 3.142, (bool)false);
  } else {
    std::cout << "Binning not yet implememented for variable '" << varName << "'. Exiting." << std::endl;
    exit(98112);
  }

}


void fitTools::drawSingleGraph(TGraph* gr, const std::string& canvasName) {

  TH2F* h = new TH2F("h_ausialiario", "", 10, 20., 1000., 10, 0.2, 1.);
  h->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h->SetStats(0);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h->Draw();
  gr->Draw("*same");
  c1->SaveAs(canvasName.c_str());

  delete c1;
  delete h;
}


void fitTools::fitProjection(TH1* h1_projection, TF1* gaussian, Float_t nSigma, std::string option, bool add) {

  Float_t histMean = h1_projection->GetMean();
  Float_t histRMS = h1_projection->GetRMS();

  gaussian->SetParameter(0, h1_projection->GetMaximum());
  gaussian->SetParameter(1, histMean);
  gaussian->SetParameter(2, histRMS);

  if (histRMS == 0.) {
    return;
  }

  gaussian->SetParLimits(1, 0., 2.*histMean);

  Float_t lowerBound = histMean - nSigma * histRMS;
  Float_t upperBound = histMean + nSigma * histRMS;

  gaussian->SetRange(lowerBound, upperBound);

  h1_projection->Fit(gaussian, option.c_str());

  int n_iter = 3;

  for (int i = 0; i < n_iter; ++i) {

    Float_t lowerBound = gaussian->GetParameter(1) - nSigma * gaussian->GetParameter(2);
    Float_t upperBound = gaussian->GetParameter(1) + nSigma * gaussian->GetParameter(2);

    gaussian->SetRange(lowerBound, upperBound);

    if (add && (i == (n_iter - 1))) {
      option = option + "+";
    }

    h1_projection->Fit(gaussian, option.c_str());

  }

}



void fitTools::fitProjection_sameArea(TH1* h1_projection, TF1* gaussian, TH1** newhisto, Float_t percIntegral, const std::string& option, bool useMode) {


  if (percIntegral < 0. || percIntegral > 1.) {
    std::cout << "WARNING! percIntegral is " << percIntegral << "!! Setting it to 90%." << std::endl;
    percIntegral = 0.9;
  }

  Int_t nBins = h1_projection->GetNbinsX();
  Double_t xMin = h1_projection->GetXaxis()->GetXmin();
  Double_t xMax = h1_projection->GetXaxis()->GetXmax();
  Double_t binWidth = (xMax - xMin) / (Double_t)nBins; //WARNING: this works only if bins are of the same size
  Double_t integral = h1_projection->Integral();

  //first: find maximum

  Int_t maxBin;
  if (useMode) {
    maxBin = h1_projection->GetMaximumBin();
  } else {
    TF1* tmp_gaussian = new TF1("tmp_gaussian", "gaus");
    fitProjection(h1_projection, tmp_gaussian, 2.5, "RQN");
    maxBin = (Int_t)ceil((tmp_gaussian->GetParameter(1) - xMin) / binWidth);
    delete tmp_gaussian;
  }

//  std::cout << "maxBin: " << maxBin << "\tbin center: " << h1_projection->GetXaxis()->GetBinCenter(maxBin) << "\t gauss mu: " << gaussian->GetParameter(1) << std::endl;
  TH1D* newHisto_tmp = new TH1D("newHisto_tmp", "", nBins, xMin, xMax);
  newHisto_tmp->SetBinContent(maxBin, h1_projection->GetBinContent(maxBin));
  newHisto_tmp->SetBinError(maxBin, h1_projection->GetBinError(maxBin));
  Int_t iBin = maxBin;
  Int_t delta_iBin = 1;
  Int_t sign  = 1;
  Float_t xMin_fit = newHisto_tmp->GetXaxis()->GetBinLowEdge(maxBin);
  Float_t xMax_fit = newHisto_tmp->GetXaxis()->GetBinUpEdge(maxBin);


  //add bins till percent area is reached:
  while (newHisto_tmp->Integral() < percIntegral * integral) {

    iBin += sign * delta_iBin;

    newHisto_tmp->SetBinContent(iBin, h1_projection->GetBinContent(iBin));
    newHisto_tmp->SetBinError(iBin, h1_projection->GetBinError(iBin));

    if (newHisto_tmp->GetXaxis()->GetBinLowEdge(iBin) < xMin_fit) {
      xMin_fit = newHisto_tmp->GetXaxis()->GetBinLowEdge(iBin);
    }
    if (newHisto_tmp->GetXaxis()->GetBinUpEdge(iBin)  > xMax_fit) {
      xMax_fit = newHisto_tmp->GetXaxis()->GetBinLowEdge(iBin);
    }

    delta_iBin += 1;
    sign *= -1; //makes it jump from left to right about max

  }

//  std::cout << "done with rms." << std::endl;
//    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
//    c1->cd();
//    h1_projection->Draw();
//    newHisto->SetFillColor(kRed);
//    newHisto->DrawClone("HISTO same");


  //initialize parameters to likely values:
  gaussian->SetParameter(0, newHisto_tmp->Integral());
  gaussian->SetParameter(1, newHisto_tmp->GetMean());
  gaussian->SetParameter(2, newHisto_tmp->GetRMS());

  gaussian->SetRange(xMin_fit, xMax_fit);
  newHisto_tmp->Fit(gaussian, option.c_str());

  *newhisto = newHisto_tmp;
}



void fitTools::getTruncatedMeanAndRMS(TH1* h1_projection, Float_t& mean, Float_t& mean_err, Float_t& rms, Float_t& rms_err, Double_t percentIntegral_MEAN, Double_t percentIntegral_RMS) {
//TCanvas* getTruncatedMeanAndRMS(TH1D* h1_projection, Float_t& mean, Float_t& mean_err, Float_t& rms, Float_t& rms_err, Double_t percentIntegral_MEAN=0.9, Double_t percentIntegral_RMS=0.68) {

  bool useMode = false;


  if (percentIntegral_MEAN < 0. || percentIntegral_MEAN > 1.) {
    std::cout << "WARNING! percentIntegral_MEAN is " << percentIntegral_MEAN << "!! Setting it to 90%." << std::endl;
    percentIntegral_MEAN = 0.9;
  }

  if (percentIntegral_RMS < 0. || percentIntegral_RMS > 1.) {
    std::cout << "WARNING! percentIntegral_RMS is " << percentIntegral_RMS << "!! Setting it to 68%." << std::endl;
    percentIntegral_RMS = 0.68;
  }

  Int_t nBins = h1_projection->GetNbinsX();
  Double_t xMin = h1_projection->GetXaxis()->GetXmin();
  Double_t xMax = h1_projection->GetXaxis()->GetXmax();
  Double_t binWidth = (xMax - xMin) / (Double_t)nBins; //WARNING: this works only if bins are of the same size
  Double_t integral = h1_projection->Integral();
//  std::cout << "xmax: " << xMax << "\txMin: " << xMin << std::endl;

  //first: find maximum
//  std::cout << "N: " << gaussian->GetParameter(0) << "\tmu: " << gaussian->GetParameter(1) << "\tsigma: " << gaussian->GetParameter(2) << std::endl;
  Int_t maxBin;
  if (useMode) {
    maxBin = h1_projection->GetMaximumBin();
  } else {
    TF1* gaussian = new TF1("gaussian", "gaus");
    gaussian->SetLineColor(kGreen);
    fitProjection(h1_projection, gaussian, 1.5, "RQN");
    maxBin = (Int_t)ceil((gaussian->GetParameter(1) - xMin) / binWidth);
    delete gaussian;
  }

//  std::cout << "maxBin: " << maxBin << "\tbin center: " << h1_projection->GetXaxis()->GetBinCenter(maxBin) << "\t gauss mu: " << gaussian->GetParameter(1) << std::endl;
  TH1D* newHisto = new TH1D("newHisto", "", nBins, xMin, xMax);
  newHisto->SetBinContent(maxBin, h1_projection->GetBinContent(maxBin));
  newHisto->SetBinError(maxBin, h1_projection->GetBinError(maxBin));
  Int_t iBin = maxBin;
  Int_t delta_iBin = 1;
  Int_t sign  = 1;
//  std::cout << "iBin: " << iBin << "\tint: " << newHisto->Integral()/integral << std::endl;

  while (newHisto->Integral() < percentIntegral_RMS * integral) {

    iBin += sign * delta_iBin;

//  std::cout << "iBin: " << iBin << "\tint: " << newHisto->Integral()/integral << std::endl;
    newHisto->SetBinContent(iBin, h1_projection->GetBinContent(iBin));
    newHisto->SetBinError(iBin, h1_projection->GetBinError(iBin));

    delta_iBin += 1;
    sign *= -1;

  }

//  std::cout << "done with rms." << std::endl;
//    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
//    c1->cd();
//    h1_projection->Draw();
//    newHisto->SetFillColor(kRed);
//    newHisto->DrawClone("HISTO same");

  rms = newHisto->GetRMS();
  rms_err = newHisto->GetRMSError();

//std::cout << "rms: " << rms << std::endl;
  while (newHisto->Integral() < percentIntegral_MEAN * integral) {
//  std::cout << "iBin: " << iBin << "\tint: " << newHisto->Integral()/integral << std::endl;

    iBin += sign * delta_iBin;

    newHisto->SetBinContent(iBin, h1_projection->GetBinContent(iBin));
    newHisto->SetBinError(iBin, h1_projection->GetBinError(iBin));

    delta_iBin += 1;
    sign *= -1;

  }

//    newHisto->SetFillStyle(3004);
//    newHisto->SetFillColor(kBlue);
//    newHisto->DrawClone("HISTO same");

  mean = newHisto->GetMean();
  mean_err = newHisto->GetMeanError();

  delete newHisto;

//    return c1;
}



void fitTools::fillProfile(TH1F* h1_response_FIT, TH1F* h1_resolution_FIT, TH1F* h1_response_MEAN, TH1F* h1_resolution_RMS, TH2D* h2, std::string name) {

  std::string fileName = "Projections/" + name + ".root";
  TFile* projectionFile;
  if (name != "") {
    projectionFile = TFile::Open(fileName.c_str(), "RECREATE");
    projectionFile->cd();
  }

  for (int iBin = 1; iBin < (h2->GetNbinsX() + 1); ++iBin) {

    char histName[50];
    sprintf(histName, "projection_%d", iBin);
    TH1D* h1_projection = h2->ProjectionY(histName, iBin, iBin);


    TF1* gaussian_LL = new TF1("gaussian_LL", "gaus");
    fitProjection(h1_projection, gaussian_LL, 2., "RQLL");

    TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
    fitProjection(h1_projection, gaussian_chi, 2., "RQO+");

    if (name != "") {
      h1_projection->Write();
    }

    Float_t mu = gaussian_LL->GetParameter(1);
    Float_t mu_err = gaussian_chi->GetParError(1);
    h1_response_FIT->SetBinContent(iBin, mu);
    h1_response_FIT->SetBinError(iBin, mu_err);

    Float_t sigma = gaussian_LL->GetParameter(2);
    Float_t resolution = (mu != 0.) ? sigma / mu : -1.;
    h1_resolution_FIT->SetBinContent(iBin, resolution);

    Float_t sigma_err = gaussian_chi->GetParError(2);
    Float_t res_err = (mu != 0.) ? sqrt(sigma_err * sigma_err / (mu * mu) + mu_err * mu_err * sigma * sigma / (mu * mu * mu * mu)) : 0.;
    h1_resolution_FIT->SetBinError(iBin, res_err);


    Float_t n = h1_projection->GetEntries();
    Float_t mean = h1_projection->GetMean();
    Float_t mean_err = (n != 0) ? h1_projection->GetRMS() / sqrt(n) : 0.;
    h1_response_MEAN->SetBinContent(iBin, mean);
    h1_response_MEAN->SetBinError(iBin, mean_err);

    Float_t rms = h1_projection->GetRMS();
    Float_t rms_err = (n != 0) ? h1_projection->GetRMS() / sqrt(n) : 0.;
    resolution = (mean != 0.) ? rms / mean : -1.;
    res_err = (mean != 0.) ? sqrt(rms_err * rms_err / (mean * mean) + mean_err * mean_err * rms * rms / (mean * mean * mean * mean)) : 0.;
    if (resolution != 0.) {
      h1_resolution_RMS->SetBinContent(iBin, resolution);
      h1_resolution_RMS->SetBinError(iBin, res_err);
    }

    h1_projection = 0;

  } //for bins

  if (name != "") {
    projectionFile->Write();
    projectionFile->Close();
    delete projectionFile;
  }
  projectionFile = 0;

} //fill profile



//new fit reponse:
void fitTools::fitDistribution_TGraph(TH2D* h2, TH2D* genMean, const std::string& varY, const std::string& varX, const std::string& etaRegion, const std::string& flag, const std::string& algoType, const std::string& outFileName, const std::string& name, Float_t percIntegral, bool use_samearea) {


  Int_t nBins = h2->GetNbinsX();

  Float_t response_FIT_x[nBins];
  Float_t response_FIT_y[nBins];
  Float_t response_FIT_xerr[nBins];
  Float_t response_FIT_yerr[nBins];

  Float_t response_MEAN_x[nBins];
  Float_t response_MEAN_y[nBins];
  Float_t response_MEAN_xerr[nBins];
  Float_t response_MEAN_yerr[nBins];

  Float_t resolution_FIT_x[nBins];
  Float_t resolution_FIT_y[nBins];
  Float_t resolution_FIT_xerr[nBins];
  Float_t resolution_FIT_yerr[nBins];

  Float_t resolution_RMS_x[nBins];
  Float_t resolution_RMS_y[nBins];
  Float_t resolution_RMS_xerr[nBins];
  Float_t resolution_RMS_yerr[nBins];

  std::string fileName = "Projections/" + name + ".root";
  TFile* projectionFile;
  if (name != "") {
    projectionFile = TFile::Open(fileName.c_str(), "RECREATE");
    projectionFile->cd();
  }



  for (int iBin = 1; iBin < nBins + 1; ++iBin) {

    char projName[100];
    sprintf(projName, "%sGenMean_%dbin", varX.c_str(), iBin);
    TH1* h1_proj = genMean->ProjectionY(projName, iBin, iBin);

    Float_t proj_mean = h1_proj->GetMean();
    Float_t proj_rms = h1_proj->GetRMS();
    Float_t proj_entries = h1_proj->GetEntries();

    response_FIT_x[iBin - 1] =  proj_mean;
    response_MEAN_x[iBin - 1] = proj_mean;
    resolution_FIT_x[iBin - 1] =  proj_mean;
    resolution_RMS_x[iBin - 1] =  proj_mean;

    response_FIT_xerr[iBin - 1] =    proj_rms / sqrt(proj_entries);
    response_MEAN_xerr[iBin - 1] =   proj_rms / sqrt(proj_entries);
    resolution_FIT_xerr[iBin - 1] =  proj_rms / sqrt(proj_entries);
    resolution_RMS_xerr[iBin - 1] =  proj_rms / sqrt(proj_entries);

    //response_FIT_xerr[iBin-1] =    ( proj_entries>1. ) ? proj_rms/sqrt(proj_entries) : proj_rms;
    //response_MEAN_xerr[iBin-1] =   ( proj_entries>1. ) ? proj_rms/sqrt(proj_entries) : proj_rms;
    //resolution_FIT_xerr[iBin-1] =  ( proj_entries>1. ) ? proj_rms/sqrt(proj_entries) : proj_rms;
    //resolution_RMS_xerr[iBin-1] =  ( proj_entries>1. ) ? proj_rms/sqrt(proj_entries) : proj_rms;


    char histName[50];
    sprintf(histName, "projection_%d", iBin);
    TH1* h1_projection = h2->ProjectionY(histName, iBin, iBin);

    TH1* h1_samearea = NULL;
    TF1* gaussian = new TF1("gaussian", "gaus");
    fitProjection_sameArea(h1_projection, gaussian, &h1_samearea, 0.95);
    //TF1* gaussian_LL = new TF1("gaussian_LL", "gaus");
    //fitProjection(h1_projection, gaussian_LL, nSigma, "RQLL");

    //TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
    //fitProjection(h1_projection, gaussian_chi, nSigma, "RQO+");

    if (name != "") {
      h1_proj->Write();
      h1_projection->Write();
    }

    Float_t mu = gaussian->GetParameter(1);
    Float_t mu_err = gaussian->GetParError(1);
    //Float_t mu = gaussian_chi->GetParameter(1);
    //Float_t mu_err = gaussian_chi->GetParError(1);
    response_FIT_y[iBin - 1] = mu;
    response_FIT_yerr[iBin - 1] = mu_err;

    Float_t sigma = gaussian->GetParameter(2);
    //Float_t sigma = gaussian_chi->GetParameter(2);
    Float_t resolution = (mu != 0.) ? sigma / mu : -1.;
    resolution_FIT_y[iBin - 1] = (varY == "response") ? resolution : sigma;

    Float_t sigma_err = gaussian->GetParError(2);
    //Float_t sigma_err = gaussian_chi->GetParError(2);
    Float_t res_err = (mu != 0.) ? sqrt(sigma_err * sigma_err / (mu * mu) + mu_err * mu_err * sigma * sigma / (mu * mu * mu * mu)) : 0.;
    resolution_FIT_yerr[iBin - 1] = (varY == "response") ? res_err : sigma_err;

    Float_t mean, rms, mean_err, rms_err;

    if (use_samearea) {

      mean = h1_projection->GetMean();
      mean_err = h1_projection->GetMeanError();
      rms = h1_projection->GetRMS();
      rms_err = h1_projection->GetRMSError();

    } else { //truncated mean/rms:

      mean = h1_samearea->GetMean();
      mean_err = h1_samearea->GetMeanError();
      rms = h1_samearea->GetRMS();
      rms_err = h1_samearea->GetRMSError();

    }

    //getTruncatedMeanAndRMS(h1_projection, mean, mean_err, rms, rms_err, 1., 0.90);

    response_MEAN_y[iBin - 1] = mean;
    response_MEAN_yerr[iBin - 1] = mean_err;

    resolution = (mean != 0.) ? rms / mean : -1.;
    res_err = (mean != 0.) ? sqrt(rms_err * rms_err / (mean * mean) + mean_err * mean_err * rms * rms / (mean * mean * mean * mean)) : 0.;
    if (rms != 0.) {
      resolution_RMS_y[iBin - 1] = (varY == "response") ? resolution : rms; // resolutions on other variables dont have to be normalized (e.g. deltaPhi)
      resolution_RMS_yerr[iBin - 1] = (varY == "response") ? res_err : rms_err;
    }

    h1_projection = 0;

    delete gaussian;
    //delete gaussian_LL;
    //delete gaussian_chi;

  } //for bins

  TGraphErrors* gr_response_FIT = new TGraphErrors(nBins, response_FIT_x, response_FIT_y, response_FIT_xerr, response_FIT_yerr);
  TGraphErrors* gr_response_MEAN = new TGraphErrors(nBins, response_MEAN_x, response_MEAN_y, response_MEAN_xerr, response_MEAN_yerr);
  TGraphErrors* gr_resolution_FIT = new TGraphErrors(nBins, resolution_FIT_x, resolution_FIT_y, resolution_FIT_xerr, resolution_FIT_yerr);
  TGraphErrors* gr_resolution_RMS = new TGraphErrors(nBins, resolution_RMS_x, resolution_RMS_y, resolution_RMS_xerr, resolution_RMS_yerr);

  std::string gr_name;

  gr_name = "gr_" + varY + "_vs_" + varX + "_FIT_" + etaRegion;
  gr_name = gr_name + "_" + algoType;
  if (flag != "") {
    gr_name = gr_name + "_" + flag;
  }
  gr_response_FIT->SetName(gr_name.c_str());

  gr_name = "gr_" + varY + "_vs_" + varX + "_MEAN_" + etaRegion;
  gr_name = gr_name + "_" + algoType;
  if (flag != "") {
    gr_name = gr_name + "_" + flag;
  }
  gr_response_MEAN->SetName(gr_name.c_str());

  gr_name = "gr_" + varY + "Res_vs_" + varX + "_FIT_" + etaRegion;
  gr_name = gr_name + "_" + algoType;
  if (flag != "") {
    gr_name = gr_name + "_" + flag;
  }
  gr_resolution_FIT->SetName(gr_name.c_str());

  gr_name = "gr_" + varY + "Res_vs_" + varX + "_RMS_" + etaRegion;
  gr_name = gr_name + "_" + algoType;
  if (flag != "") {
    gr_name = gr_name + "_" + flag;
  }
  gr_resolution_RMS->SetName(gr_name.c_str());


  TFile* outFile = TFile::Open(outFileName.c_str(), "update");
  outFile->cd();
  gr_response_FIT->Write();
  gr_response_MEAN->Write();
  gr_resolution_FIT->Write();
  gr_resolution_RMS->Write();

  outFile->Write();
  outFile->Close();

  if (name != "") {
    projectionFile->Write();
    projectionFile->Close();
    delete projectionFile;
  }
  projectionFile = 0;

} //fill profile


//old fit response:
//now deprecated
/*
void fitDistribution_TGraph(TH2D* h2, TProfile* genMean, const std::string& varX, const std::string& etaRegion, const std::string& flag, const std::string& algoType, const std::string& outFileName, const std::string& name="") {


  Int_t nBins = h2->GetNbinsX();

  Float_t response_FIT_x[nBins];
  Float_t response_FIT_y[nBins];
  Float_t response_FIT_xerr[nBins];
  Float_t response_FIT_yerr[nBins];

  Float_t response_MEAN_x[nBins];
  Float_t response_MEAN_y[nBins];
  Float_t response_MEAN_xerr[nBins];
  Float_t response_MEAN_yerr[nBins];

  Float_t resolution_FIT_x[nBins];
  Float_t resolution_FIT_y[nBins];
  Float_t resolution_FIT_xerr[nBins];
  Float_t resolution_FIT_yerr[nBins];

  Float_t resolution_RMS_x[nBins];
  Float_t resolution_RMS_y[nBins];
  Float_t resolution_RMS_xerr[nBins];
  Float_t resolution_RMS_yerr[nBins];

  std::string fileName = "Projections/"+name+".root";
  TFile* projectionFile;
  if( name!= "" ) {
    projectionFile = TFile::Open(fileName.c_str(), "RECREATE");
    projectionFile->cd();
  }



  for(int iBin=1; iBin<nBins+1; ++iBin) {

    response_FIT_x[iBin-1] =  genMean->GetBinContent(iBin);
    response_MEAN_x[iBin-1] = genMean->GetBinContent(iBin);
    resolution_FIT_x[iBin-1] =  genMean->GetBinContent(iBin);
    resolution_RMS_x[iBin-1] =  genMean->GetBinContent(iBin);

    //response_FIT_xerr[iBin-1] =  genMean->GetBinError(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));
    //response_MEAN_xerr[iBin-1] = genMean->GetBinError(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));
    //resolution_FIT_xerr[iBin-1] =  genMean->GetBinError(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));
    //resolution_RMS_xerr[iBin-1] =  genMean->GetBinError(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));

    response_FIT_xerr[iBin-1] =  genMean->GetBinError(iBin);
    response_MEAN_xerr[iBin-1] = genMean->GetBinError(iBin);
    resolution_FIT_xerr[iBin-1] =  genMean->GetBinError(iBin);
    resolution_RMS_xerr[iBin-1] =  genMean->GetBinError(iBin);

    //response_FIT_xerr[iBin-1] =  genMean->GetBinContent(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));
    //response_MEAN_xerr[iBin-1] = genMean->GetBinContent(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));
    //resolution_FIT_xerr[iBin-1] =  genMean->GetBinContent(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));
    //resolution_RMS_xerr[iBin-1] =  genMean->GetBinContent(iBin)/sqrt((Float_t)genMean->GetBinEntries(iBin));


    char histName[50];
    sprintf(histName, "projection_%d",iBin);
    TH1D* h1_projection = h2->ProjectionY(histName, iBin, iBin);

    TF1* gaussian_LL = new TF1("gaussian_LL", "gaus");
    fitProjection(h1_projection, gaussian_LL, 2., "RQLL");

    TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
    fitProjection(h1_projection, gaussian_chi, 2., "RQN");

    if( name!="" ) {
        h1_projection->Write();
    }

    Float_t mu = gaussian_LL->GetParameter(1);
    Float_t mu_err = gaussian_chi->GetParError(1);
    response_FIT_y[iBin-1] = mu;
    response_FIT_yerr[iBin-1] = mu_err;

    Float_t sigma = gaussian_LL->GetParameter(2);
    Float_t resolution = (mu!=0.) ? sigma/mu : -1.;
    resolution_FIT_y[iBin-1] = resolution;

    Float_t sigma_err = gaussian_chi->GetParError(2);
    Float_t res_err = (mu!=0.) ? sqrt( sigma_err*sigma_err/(mu*mu) + mu_err*mu_err*sigma*sigma/(mu*mu*mu*mu) ) : 0.;
    resolution_FIT_yerr[iBin-1] = res_err;


    Float_t n = h1_projection->GetEntries();
    Float_t mean = h1_projection->GetMean();
    Float_t mean_err = (n!=0) ? h1_projection->GetRMS()/sqrt(n) : 0.;
    response_MEAN_y[iBin-1] = mean;
    response_MEAN_yerr[iBin-1] = mean_err;

    Float_t rms = h1_projection->GetRMS();
    Float_t rms_err = (n!=0) ? h1_projection->GetRMS()/sqrt(n) : 0.;
    resolution = (mean!=0.) ? rms/mean : -1.;
    res_err = (mean!=0.) ? sqrt( rms_err*rms_err/(mean*mean) + mean_err*mean_err*rms*rms/(mean*mean*mean*mean) ) : 0.;
    if( resolution != 0. ) {
      resolution_RMS_y[iBin-1] = resolution;
      resolution_RMS_yerr[iBin-1] = res_err;
    }

    h1_projection = 0;

    delete gaussian_LL;
    delete gaussian_chi;

  } //for bins

  TGraphErrors* gr_response_FIT = new TGraphErrors(nBins, response_FIT_x, response_FIT_y, response_FIT_xerr, response_FIT_yerr);
  TGraphErrors* gr_response_MEAN = new TGraphErrors(nBins, response_MEAN_x, response_MEAN_y, response_MEAN_xerr, response_MEAN_yerr);
  TGraphErrors* gr_resolution_FIT = new TGraphErrors(nBins, resolution_FIT_x, resolution_FIT_y, resolution_FIT_xerr, resolution_FIT_yerr);
  TGraphErrors* gr_resolution_RMS = new TGraphErrors(nBins, resolution_RMS_x, resolution_RMS_y, resolution_RMS_xerr, resolution_RMS_yerr);

  std::string gr_name;

  gr_name = "gr_response_vs_" + varX + "_FIT_" + etaRegion;
  gr_name = gr_name + "_" + algoType;
  if(flag!="")
    gr_name = gr_name + "_" + flag;
  gr_response_FIT->SetName(gr_name.c_str());

  gr_name = "gr_response_vs_" + varX + "_MEAN_"+etaRegion;
  gr_name = gr_name + "_" + algoType;
  if(flag!="")
    gr_name = gr_name + "_" + flag;
  gr_response_MEAN->SetName(gr_name.c_str());

  gr_name = "gr_resolution_vs_" + varX + "_FIT_"+etaRegion;
  gr_name = gr_name + "_" + algoType;
  if(flag!="")
    gr_name = gr_name + "_" + flag;
  gr_resolution_FIT->SetName(gr_name.c_str());

  gr_name = "gr_resolution_vs_" + varX + "_RMS_"+etaRegion;
  gr_name = gr_name + "_" + algoType;
  if(flag!="")
    gr_name = gr_name + "_" + flag;
  gr_resolution_RMS->SetName(gr_name.c_str());


  TFile* outFile = TFile::Open(outFileName.c_str(), "update");
  outFile->cd();
  gr_response_FIT->Write();
  gr_response_MEAN->Write();
  gr_resolution_FIT->Write();
  gr_resolution_RMS->Write();

  outFile->Write();
  outFile->Close();

  if(name!="") {
    projectionFile->Write();
    projectionFile->Close();
    delete projectionFile;
  }
  projectionFile = 0;

} //fill profile
*/

void fitTools::fillPositionResolution(TH1F* h1_sigmaEta, TH1F* h1_sigmaPhi, TH2D* h2_deltaEta, TH2D* h2_deltaPhi) {

  for (int iBin = 1; iBin < (h2_deltaEta->GetNbinsX() + 1); ++iBin) {

    TH1D* h1_projection = h2_deltaEta->ProjectionY("projectiony", iBin, iBin);


    TF1* gaussian = new TF1("gaussian", "gaus");
    fitProjection(h1_projection, gaussian, 2., "RQLL");

    h1_sigmaEta->SetBinContent(iBin, gaussian->GetParameter(2));
    h1_sigmaEta->SetBinError(iBin, gaussian->GetParError(2));

    delete gaussian;
    gaussian = 0;

    h1_projection = 0;

  } //for bins eta

  for (int iBin = 1; iBin < (h2_deltaPhi->GetNbinsX() + 1); ++iBin) {

    TH1D* h1_projection = h2_deltaPhi->ProjectionY("projectiony", iBin, iBin);


    TF1* gaussian = new TF1("gaussian", "gaus");
    fitProjection(h1_projection, gaussian, 2., "RQLL");

    h1_sigmaPhi->SetBinContent(iBin, gaussian->GetParameter(2));
    h1_sigmaPhi->SetBinError(iBin, gaussian->GetParError(2));

    delete gaussian;
    gaussian = 0;

    h1_projection = 0;

  } //for bins

} //fill position resolution



//used by getEfficiencyHisto (later on):
int fitTools::getEfficiencyUncertainties(int n, int k, double p, double& xmin, double& xmax) {

  // create a histogram with binomial distribution
  TH1D* h1_binomial = new TH1D("hist_binomial", "", 1000, 0., 1.);

  // loop over bins and fill histogram
  for (int i = 1; i <= 1000; ++i) {
    double x   = h1_binomial -> GetBinCenter(i);
    double val = TMath::Binomial(n, k) * TMath::Power(x, double(k)) * TMath::Power(1 - x, double(n - k));
    h1_binomial -> SetBinContent(i, val);
  }

  // normalize
  h1_binomial -> Scale(1.0 / h1_binomial -> Integral());

  // calculate quantiles
  int nprobSum = 4;
  double q[4];
  double probSum[4];
  probSum[0] = (1. - p) / 2.;
  probSum[1] = 1. - (1. - p) / 2.;
  probSum[2] = 0.05;
  probSum[3] = 0.95;

  h1_binomial -> GetQuantiles(nprobSum, q, probSum);
//  delete h1_binomial;
//  h1_binomial=0;

  double xexp = double(k) / double(n);
  if (xexp > q[1])  {
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();
    h1_binomial->Draw();
    c1->SaveAs("prova.eps");
    exit(11);
  }

  // calculate uncertainties
  if (n == 0) {
    xmin = 0.0;
    xmax = 0.0;
    return -3;
  } else if (xexp < q[0]) {
    xmin = 0;
    xmax = q[3];
    return -2;
  }

  else if (xexp > q[1]) {
    xmin = q[2];
    xmax = 1.0;
    return -1;
  } else {
    xmin = q[0];
    xmax = q[1];
    return 1;
  }

}


TGraphAsymmErrors* fitTools::getEfficiencyGraph(const std::string& name, TH1F* h1_numerator, TH1F* h1_denominator) {

  TGraphAsymmErrors* gr_returnGraph = new TGraphAsymmErrors();
  gr_returnGraph->SetName(name.c_str());

  int npoints = 0;

  // set points
  for (int i = 1; i <= h1_denominator -> GetNbinsX(); ++i) {

    // calculate uncertainties
    double xmin;
    double xmax;
    int flag = fitTools::getEfficiencyUncertainties(
                 int(h1_denominator -> GetBinContent(i)),
                 int(h1_numerator -> GetBinContent(i)),
                 0.68, xmin, xmax);

    if (flag == 1) {
      gr_returnGraph -> SetPoint(
        npoints,
        h1_denominator -> GetBinCenter(i),
        h1_numerator -> GetBinContent(i) / h1_denominator -> GetBinContent(i));
      // set uncertainties
      gr_returnGraph -> SetPointEXhigh(npoints, 0.);
      gr_returnGraph -> SetPointEXlow(npoints, 0.);
      gr_returnGraph -> SetPointEYhigh(npoints, xmax - h1_numerator -> GetBinContent(i) / h1_denominator -> GetBinContent(i));
      gr_returnGraph -> SetPointEYlow(npoints++, h1_numerator -> GetBinContent(i) / h1_denominator -> GetBinContent(i) - xmin);
    } else if (flag == -2) {
      gr_returnGraph -> SetPoint(npoints, h1_denominator -> GetBinCenter(i), 0.);
      // set uncertainties
      gr_returnGraph -> SetPointEXhigh(npoints, 0.);
      gr_returnGraph -> SetPointEXlow(npoints, 0.);
      gr_returnGraph -> SetPointEYhigh(npoints, xmax);
      gr_returnGraph -> SetPointEYlow(npoints++, 0.);
    } else if (flag == -1) {
      gr_returnGraph -> SetPoint(npoints, h1_denominator -> GetBinCenter(i), 1.);
      // set uncertainties
      gr_returnGraph -> SetPointEXhigh(npoints, 0.);
      gr_returnGraph -> SetPointEXlow(npoints, 0.);
      gr_returnGraph -> SetPointEYhigh(npoints, 0.);
      gr_returnGraph -> SetPointEYlow(npoints++, 1. - xmin);
    }
  } //for bins

  return gr_returnGraph;

}




TF1* fitTools::fitResponseGraph(TGraphErrors* graph, std::string funcType, std::string funcName, const std::string& option, float rangeMax, float rangeMin) {


  TF1* fitFunction;

  if (funcType == "rpf") {
    fitFunction = new TF1(funcName.c_str(), rpf, rangeMin, rangeMax, 6); //will have to fix the range issue!
    fitFunction->SetParameters(100, 0.85, 4.2, 80, 250, 1.);
    fitFunction->SetParLimits(1, 0.5, 1.0);
    fitFunction->SetParLimits(2, 1., 10.);
    fitFunction->SetParLimits(4, 100., 500.);
    fitFunction->SetParameter(0, 0.6);
    fitFunction->SetParameter(5, 1.);
    //fitFunction->SetParameter(1, -0.6);
  } else if (funcType == "powerlaw") {
    fitFunction = new TF1(funcName.c_str(), "[0] - [1]/(pow(x, [2]))");
    fitFunction->SetRange(rangeMin, rangeMax);
    fitFunction->SetParameters(1., 1., 0.3);
  } else if (funcType == "powerlawL2L3") {
    fitFunction = new TF1(funcName.c_str(), "1. + [0]/(pow(x, [1]))");
    fitFunction->SetRange(rangeMin, rangeMax);
    fitFunction->SetParameters(0.7, 0.7);
    //fitFunction->SetParLimits(0, 0.5, 2.);
    //fitFunction->SetParLimits(1, 0.5, 2.);
  } else if (funcType == "powerlaw_corr") {
    fitFunction = new TF1(funcName.c_str(), "[0] - [1]/(pow(x, [2])) + [3]/x");
    fitFunction->SetRange(rangeMin, rangeMax);
    fitFunction->SetParameters(1., 1., 0.3, 1.);
  } else {
    std::cout << "Function '" << funcType << "' not implemented yet for fitResponseGraph. Exiting." << std::endl;
    exit(119);
  }

  graph->Fit(fitFunction, option.c_str());

  return fitFunction;

}



TF1* fitTools::fitResolutionGraph(TGraphErrors* graph, std::string funcType, std::string funcName, const std::string& option, float rangeMax, float rangeMin) {


  TF1* fitFunction;

  if (funcType == "NSC") {
    fitFunction = new TF1(funcName.c_str(), NSC, rangeMin, rangeMax, 3);
    fitFunction->SetParameters(0., 1., 0.05);
    fitFunction->SetParLimits(0, -2., 11.);
    fitFunction->SetParLimits(1, -2., 2.);
    fitFunction->SetParLimits(2, -0.3, 0.3);
  } else if (funcType == "NSCPF") {
    fitFunction = new TF1(funcName.c_str(), NSCPF, rangeMin, rangeMax, 4);
    fitFunction->SetParameters(0., 1., 0.05, 0.);
    fitFunction->SetParLimits(0, -2., 5.);
    fitFunction->SetParLimits(1, 0., 1.1);
    fitFunction->SetParLimits(2, 0., 0.12);
    fitFunction->SetParLimits(3, -1., 1.);
  } else {
    std::cout << "Function '" << funcType << "' not implemented yet for fitResolutionGraph. Exiting." << std::endl;
    exit(119);
  }

  graph->Fit(fitFunction, option.c_str());

  return fitFunction;

}


Double_t rpf(Double_t* x, Double_t* p) {

  double pt = x[0];
  double a = p[0];
  double m = p[1];
  double r = p[2];
  double x1 = std::max(std::min(p[3], p[4]), 1.);
  double x2 = std::max(std::max(p[3], p[4]), 1.);
  double r_inf = p[5];

  double f = a * pow(pt, -m);

  // Transition roughly from 80 to 250
  double xmid = 0.5 * (x2 + x1);
  double xwid = 2 * 0.5 * std::max(x2 - x1, 1.);
  double w = 0.5 * (1 + TMath::Erf((pt - xmid) / xwid));

  // Logarithmic transition. Linear seems to work better, though
  //double xmid = 0.5*(log(x2) + log(x1));
  //double xwid = log(x2/sqrt(x1*x2));
  //double w = 0.5 * (1 + TMath::Erf( (log(pt) - xmid) / xwid ));

  //double w = (pt-x1)/max(x2-x1,1.);
  //double w = log(pt/x1)/log(x2/x1);
  if (w < 0.) {
    w = 0.;
  }
  if (w > 1.) {
    w = 1.;
  }

  return ((1. - w) * (r_inf - f) + w * (r_inf - r * f));

}


Double_t NSC(Double_t* x, Double_t* p) {

  double pt = x[0];
  double N = p[0];
  double S = p[1];
  double C = p[2];

  return sqrt(N * N / pt / pt + S * S / pt + C * C);

}

Double_t NSCPF(Double_t* x, Double_t* p) {

  /*
    double pt = x[0];
    double N = 0.;
    double S = p[0];
    double C = p[1];

    return sqrt( N*N/pt/pt + S*S/pt + C*C );
  */

  double pt = x[0];
  double N = p[0];
  double S = p[1];
  double C = p[2];
  double m = p[3];

  return sqrt(N * fabs(N) / pt / pt + S * S * pow(pt, m - 1.) + C * C);

}



Double_t powerlaw(Double_t* x, Double_t* p) {

  Double_t pt = x[0];

  Double_t value = p[0] - p[1] / (pow(pt, p[2])) + p[3] / pt;

  return value;

}




TH1D* fitTools::getBand(TF1* f, const std::string& name) {

  const int ndim_resp_q = f->GetNpar();
  TMatrixD emat_resp_q(ndim_resp_q, ndim_resp_q);
  gMinuit->mnemat(&emat_resp_q[0][0], ndim_resp_q);

  return getBand(f, emat_resp_q, name);

}



// Create uncertainty band (histogram) for a given function and error matrix
// in the range of the function.
TH1D* fitTools::getBand(TF1* f, TMatrixD const& m, std::string name, bool getRelativeBand, int npx) {

  Bool_t islog = true;
  //double xmin = f->GetXmin()*0.9;
  //double xmax = f->GetXmax()*1.1; //fixes problem in drawing with c option
  double xmin = f->GetXmin();
  double xmax = f->GetXmax() * 1.1; //fixes problem in drawing with c option
  int npar = f->GetNpar();
  //TString formula = f->GetExpFormula();

  // Create binning (linear or log)
  Double_t xvec[npx];
  xvec[0] = xmin;
  double dx = (islog ? pow(xmax / xmin, 1. / npx) : (xmax - xmin) / npx);
  for (int i = 0; i != npx; ++i) {
    xvec[i + 1] = (islog ? xvec[i] * dx : xvec[i] + dx);
  }


  //
  // Compute partial derivatives numerically
  // can be used with any fit function
  //
  Double_t sigmaf[npx];
  TH1D* h1_band = new TH1D(name.c_str(), "", npx, xvec);

  for (int ipx = 0; ipx < npx; ++ipx) {

    sigmaf[ipx] = 0.;
    Double_t partDeriv[npar];

    //compute partial derivatives of f wrt its parameters:
    for (int ipar = 0; ipar < npar; ++ipar) {

      Float_t pi = f->GetParameter(ipar);
      Float_t dpi = sqrt(m[ipar][ipar]) * 0.01; //small compared to the par sigma
      f->SetParameter(ipar, pi + dpi);
      Float_t fplus = f->Eval(xvec[ipx]);
      f->SetParameter(ipar, pi - dpi);
      Float_t fminus = f->Eval(xvec[ipx]);
      f->SetParameter(ipar, pi); //put it back as it was

      partDeriv[ipar] = (fplus - fminus) / (2.*dpi);

    } //for params

    //compute sigma(f) at x:
    for (int ipar = 0; ipar < npar; ++ipar) {
      for (int jpar = 0; jpar < npar; ++jpar) {
        sigmaf[ipx] += partDeriv[ipar] * partDeriv[jpar] * m[ipar][jpar];
      }
    }
    sigmaf[ipx] = sqrt(sigmaf[ipx]); //absolute band

    h1_band->SetBinContent(ipx, f->Eval(xvec[ipx]));
    if (getRelativeBand) {
      h1_band->SetBinError(ipx, sigmaf[ipx] / f->Eval(xvec[ipx]));
    } else {
      h1_band->SetBinError(ipx, sigmaf[ipx]);
    }

  } //for points

  h1_band->SetMarkerStyle(20);
  h1_band->SetMarkerSize(0);
  h1_band->SetFillColor(kYellow - 9);


  //TGraph* h1_statError = new TGraph(npx, xvec, sigmaf);
//TH2D* h2_axesStat = new TH2D("axesStat", "", 10, 20., 1400., 10, 0., 10.);
//h2_axesStat->GetXaxis()->SetNoExponent();
//h2_axesStat->GetXaxis()->SetMoreLogLabels();
//TCanvas* cStat = new TCanvas("cStat", "cStat", 600, 600);
//cStat->cd();
//cStat->SetLogx();
//h2_axesStat->Draw();
//h1_band->Draw("psame");
//std::string canvasName = "stat/" + name + ".eps";
//cStat->SaveAs(canvasName.c_str());

//delete h2_axesStat;
//delete cStat;

  return h1_band;

} //getband



TGraphErrors* fitTools::get_graphRatio(TGraphErrors* gr_data, TGraphErrors* gr_MC) {

  TGraphErrors* gr_ratio = new TGraphErrors(0);


  for (int i = 0; i < gr_data->GetN(); ++i) {

    Double_t datax, datay;
    gr_data->GetPoint(i, datax, datay);
    //Double_t dataxerr = gr_data->GetErrorX(i);
    Double_t datayerr = gr_data->GetErrorY(i);

    Double_t mcx, mcy;
    gr_MC->GetPoint(i, mcx, mcy);
    Double_t mcxerr = gr_MC->GetErrorX(i);
    Double_t mcyerr = gr_MC->GetErrorY(i);

    Double_t ratiox = mcx;
    Double_t ratioxerr = mcxerr;

    Double_t ratioy = (mcy > 0.) ? datay / mcy : 0.;
    Double_t ratioyerr = (mcy > 0.) ? sqrt(datayerr * datayerr / (mcy * mcy) + datay * datay * mcyerr * mcyerr / (mcy * mcy * mcy * mcy)) : 0.;


    if (ratioyerr > 0.) {
      gr_ratio->SetPoint(i, ratiox, ratioy);
      gr_ratio->SetPointError(i, ratioxerr, ratioyerr);
    }

  } //for points

  return gr_ratio;

}



TGraphAsymmErrors* fitTools::getGraphPoissonErrors(TH1* histo, const std::string xerrType, float nSigma) {


  TGraphAsymmErrors* graph = new TGraphAsymmErrors(0);

  for (int iBin = 1; iBin < (histo->GetXaxis()->GetNbins() + 1); ++iBin) {

    int y; // these are data histograms, so y has to be integer
    double x, xerr, yerrplus, yerrminus;
    //xerr = 0.; //no xerr for now (maybe binwidth / sqrt(12)?)

    x = histo->GetBinCenter(iBin);
    if (xerrType == "0") {
      xerr = 0.;
    } else if (xerrType == "binWidth") {
      xerr = histo->GetBinWidth(iBin) / 2.;
    } else if (xerrType == "sqrt12") {
      xerr = histo->GetBinWidth(iBin) / sqrt(12.);
    } else {
      std::cout << "Unkown xerrType '" << xerrType << "'. Setting to bin width." << std::endl;
      xerr = histo->GetBinWidth(iBin);
    }
    y = (int)histo->GetBinContent(iBin);

    double ym, yp;
    RooHistError::instance().getPoissonInterval(y, ym, yp, nSigma);

    yerrplus = yp - y;
    yerrminus = y - ym;

    /*
        // and now poissonian errors (bayes flat prior):
        if( y==0 ) {
          yerrminus = 0.;
          yerrplus = 0.5*TMath::ChisquareQuantile(cl, 2.*(y+1.) );
        } else {
          //float lowL = 0.5*TMath::ChisquareQuantile(1.-cl/2., 2.*y);
          //float upL = 0.5*TMath::ChisquareQuantile(cl/2., 2.*(y+1.) );
          float lowL = 0.5*TMath::ChisquareQuantile(1.-cl, 2.*y);
          float upL = 0.5*TMath::ChisquareQuantile(cl, 2.*(y+1.) );
          yerrminus = y - lowL;
          yerrplus = upL - y;
        }
    */
    int thisPoint = graph->GetN();
    graph->SetPoint(thisPoint, x, y);
    graph->SetPointError(thisPoint, xerr, xerr, yerrminus, yerrplus);

  }

  return graph;

}
