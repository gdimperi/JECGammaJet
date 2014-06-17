#include "GaussianProfile.h"

#include <sstream>
#include <TH1D.h>
#include <TF1.h>


void GaussianProfile::createProfiles(TFileDirectory& dir) {
  for (int i = 0; i < m_nXBins; i++) {

    std::stringstream ss;
    ss << m_name << "_" << m_prefix << "_" << getBinLowEdge(i) << "_" << getBinHighEdge(i);

    int nBins = m_nYBins;
    double min = m_YMin, max = m_YMax;

    if (m_autoBinning) {
      nBins = 100;
      min = getBinLowEdge(i) * (1 - m_autoBinningLowPercent);
      max = getBinHighEdge(i) * (1 + m_autoBinningHighPercent);
    }

    TH1* object = dir.make<TH1D>(ss.str().c_str(), ss.str().c_str(), nBins, min, max);
    m_profiles.push_back(object);
  }
}

void GaussianProfile::createGraph() {

  if (m_profiles.size() == 0 || (m_graph.get() && !m_dirty))
    return;

  std::stringstream ss;
  ss << m_name << "_graph";

  m_graph.reset(new TGraphErrors(m_nXBins));
  m_graph->SetName(ss.str().c_str());

  // Create gaussian for fitting
  //TF1* gauss = new TF1("g", "gaus");

  for (int i = 0; i < m_nXBins; i++) {
    TH1* hist = m_profiles[i];
    double x_mean = (getBinLowEdge(i) + getBinHighEdge(i)) / 2.;
    //double min = hist->GetXaxis()->GetBinLowEdge(1);
    //double max = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetLast());

    float mean = 0, mean_error = 0;
    float rms  = 0, rms_error  = 0;
    getTruncatedMeanRMS(hist, mean, mean_error, rms, rms_error);

    /*
    double min = hist->GetXaxis()->GetXmin();
    double max = hist->GetXaxis()->GetXmax();
    if (m_autoBinning) {
      min = m_XBins[i] * 0.90;
      max = m_XBins[i + 1] * 1.10;
    }

    gauss->SetRange(min, max);
    gauss->SetParameter(1, (min + max) / 2.);
    gauss->SetParLimits(1, min, max);
    gauss->SetParameter(2, 20);

    hist->Fit(gauss, "QR");
    */

    m_graph->SetPoint(i, x_mean, mean);
    m_graph->SetPointError(i, 0, mean_error);
  }

  m_dirty = false;
}

void GaussianProfile::getTruncatedMeanRMS(TH1* hist, float& mean, float& mean_error, float& rms, float& rms_error) {
  int nBins = hist->GetNbinsX();
  double xMin = hist->GetXaxis()->GetXmin();
  double xMax = hist->GetXaxis()->GetXmax();
  //double binWidth = (xMax - xMin) / (double) nBins; //WARNING: this works only if bins are of the same size
  double integral = hist->Integral();

  int maxBin = 0;
  TF1* gaussian = new TF1("gaussian", "gaus");
  fitProjection(hist, gaussian, 1.5, "RQN");
  //maxBin = (int) ceil((gaussian->GetParameter(1) - xMin) / binWidth);
  maxBin = hist->FindBin(gaussian->GetParameter(1));
  delete gaussian;

  TH1D* newHisto = new TH1D("newHisto", "", nBins, xMin, xMax);
  newHisto->SetBinContent(maxBin, hist->GetBinContent(maxBin));
  newHisto->SetBinError(maxBin, hist->GetBinError(maxBin));
  int iBin = maxBin;
  int delta_iBin = 1;
  int sign  = 1;

  while (newHisto->Integral() < 0.99 * integral) {
    iBin += sign * delta_iBin;

    newHisto->SetBinContent(iBin, hist->GetBinContent(iBin));
    newHisto->SetBinError(iBin, hist->GetBinError(iBin));

    delta_iBin += 1;
    sign *= -1;
  }

  rms = newHisto->GetRMS();
  rms_error = newHisto->GetRMSError();

  while (newHisto->Integral() < 0.99 * integral) {
    iBin += sign * delta_iBin;

    newHisto->SetBinContent(iBin, hist->GetBinContent(iBin));
    newHisto->SetBinError(iBin, hist->GetBinError(iBin));

    delta_iBin += 1;
    sign *= -1;
  }

  mean = newHisto->GetMean();
  mean_error = newHisto->GetMeanError();

  delete newHisto;
}

void GaussianProfile::fitProjection(TH1* hist, TF1* gaussian, float nSigma, const std::string& option) {
  float histMean = hist->GetMean();
  float histRMS = hist->GetRMS();

  gaussian->SetParameter(0, hist->GetMaximum());
  gaussian->SetParameter(1, histMean);
  gaussian->SetParameter(2, histRMS);

  if (histRMS == 0.) {
    return;
  }

  gaussian->SetParLimits(1, 0., 2.*histMean);

  float lowerBound = histMean - nSigma * histRMS;
  float upperBound = histMean + nSigma * histRMS;

  gaussian->SetRange(lowerBound, upperBound);

  hist->Fit(gaussian, option.c_str());

  int n_iter = 3;

  for (int i = 0; i < n_iter; ++i) {
    float lowerBound = gaussian->GetParameter(1) - nSigma * gaussian->GetParameter(2);
    float upperBound = gaussian->GetParameter(1) + nSigma * gaussian->GetParameter(2);

    gaussian->SetRange(lowerBound, upperBound);

    hist->Fit(gaussian, option.c_str());
  }
}
