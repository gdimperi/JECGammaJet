#include "PUReweighter.h"

#include <TFile.h>
#include <iostream>

PUReweighter::PUReweighter(const std::string& dataFilePath, const std::string& mcFilePath):
  puHisto(NULL) {

    TFile* dataFile = TFile::Open(dataFilePath.c_str());
    TFile* mcFile = TFile::Open(mcFilePath.c_str());

    if (! dataFile) {
      std::cerr << "Error: can't open " << dataFilePath << ". No PU reweighting." << std::endl;
      return;
    }

    if (! mcFile) {
      std::cerr << "Error: can't open " << mcFilePath << ". No PU reweighting." << std::endl;
      return;
    }

    TH1* dataHisto = static_cast<TH1*>(dataFile->Get("pileup"));
    TH1* mcHisto = static_cast<TH1*>(mcFile->Get("pileup"));

    //TODO: Check for NULL ptr

    // Normalize
    dataHisto->Scale(1.0 / dataHisto->Integral());
    mcHisto->Scale(1.0 / mcHisto->Integral());

    // MC * data / MC = data, so the weights are data/MC:
    puHisto = static_cast<TH1*>(dataHisto->Clone());
    puHisto->Divide(mcHisto);
    puHisto->SetDirectory(0); // "detach" the histo from the file

    /*
    std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;

    int NBins = puHisto->GetNbinsX();

    for (int ibin = 1; ibin < NBins + 1; ++ibin) {
      std::cout << "   " << ibin - 1 << " " << puHisto->GetBinContent(ibin) << std::endl;
    }
    */

    dataFile->Close();
    mcFile->Close();

    delete dataFile;
    delete mcFile;
  }

PUReweighter::PUReweighter(const std::string& dataFilePath, PUProfile profile/* = PUProfile::S10*/):
  puHisto(NULL) {

    TFile* dataFile = TFile::Open(dataFilePath.c_str());

    if (! dataFile) {
      std::cerr << "Error: can't open " << dataFilePath << ". No PU reweighting." << std::endl;
      return;
    }

    initPUProfiles();
    std::vector<double>& profile_coefs = mPUCoefs[profile];

    TH1* dataHisto = static_cast<TH1*>(dataFile->Get("pileup"));

    // Create MC PU histogram
    TH1* mcHisto = static_cast<TH1*>(dataHisto->Clone("mc_pileup"));
    mcHisto->Reset();
    mcHisto->SetDirectory(NULL);

    for (int i = 1; i <= dataHisto->GetNbinsX(); i++) {
      int index = static_cast<int>(dataHisto->GetBinLowEdge(i));
      double coef = (index - 1) < (int)profile_coefs.size() ? profile_coefs[index - 1] : 0.;
      if (profile == PUProfile::S7 && index <= 4)
        coef = 0; // For low PU runs

      mcHisto->SetBinContent(i, coef);
    }

    //TODO: Check for NULL ptr

    // Normalize
    dataHisto->Scale(1.0 / dataHisto->Integral());
    mcHisto->Scale(1.0 / mcHisto->Integral());

    // MC * data / MC = data, so the weights are data/MC:
    puHisto = static_cast<TH1*>(dataHisto->Clone());
    puHisto->Divide(mcHisto);
    puHisto->SetDirectory(NULL); // "detach" the histo from the file

    /*
    std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;

    int NBins = puHisto->GetNbinsX();

    for (int ibin = 1; ibin < NBins + 1; ++ibin) {
      std::cout << "   " << ibin - 1 << " " << puHisto->GetBinContent(ibin) << std::endl;
    }
    */

    dataFile->Close();

    /*
    static int i = 1;
    TString tmp = TString::Format("mc_pileup_%d.root", i);
    TFile* f = TFile::Open(tmp, "recreate");
    mcHisto->Write();
    puHisto->Write();
    f->Close();
    delete f;
    i++;
    */

    delete dataFile;
    delete mcHisto;
  }

double PUReweighter::weight(float interactions) const {
  if (!puHisto) {
    return 1.;
  }

  int bin = puHisto->GetXaxis()->FindBin(interactions);
  return puHisto->GetBinContent(bin);
} 

void PUReweighter::initPUProfiles() {

  mPUCoefs[PUProfile::S6] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };

  mPUCoefs[PUProfile::S7] = {
    2.344E-05,
    2.344E-05,
    2.344E-05,
    2.344E-05,
    4.687E-04,
    4.687E-04,
    7.032E-04,
    9.414E-04,
    1.234E-03,
    1.603E-03,
    2.464E-03,
    3.250E-03,
    5.021E-03,
    6.644E-03,
    8.502E-03,
    1.121E-02,
    1.518E-02,
    2.033E-02,
    2.608E-02,
    3.171E-02,
    3.667E-02,
    4.060E-02,
    4.338E-02,
    4.520E-02,
    4.641E-02,
    4.735E-02,
    4.816E-02,
    4.881E-02,
    4.917E-02,
    4.909E-02,
    4.842E-02,
    4.707E-02,
    4.501E-02,
    4.228E-02,
    3.896E-02,
    3.521E-02,
    3.118E-02,
    2.702E-02,
    2.287E-02,
    1.885E-02,
    1.508E-02,
    1.166E-02,
    8.673E-03,
    6.190E-03,
    4.222E-03,
    2.746E-03,
    1.698E-03,
    9.971E-04,
    5.549E-04,
    2.924E-04,
    1.457E-04,
    6.864E-05,
    3.054E-05,
    1.282E-05,
    5.081E-06,
    1.898E-06,
    6.688E-07,
    2.221E-07,
    6.947E-08,
    2.047E-08
  };

  mPUCoefs[PUProfile::S10] = {
    2.560E-06,
    5.239E-06,
    1.420E-05,
    5.005E-05,
    1.001E-04,
    2.705E-04,
    1.999E-03,
    6.097E-03,
    1.046E-02,
    1.383E-02,
    1.685E-02,
    2.055E-02,
    2.572E-02,
    3.262E-02,
    4.121E-02,
    4.977E-02,
    5.539E-02,
    5.725E-02,
    5.607E-02,
    5.312E-02,
    5.008E-02,
    4.763E-02,
    4.558E-02,
    4.363E-02,
    4.159E-02,
    3.933E-02,
    3.681E-02,
    3.406E-02,
    3.116E-02,
    2.818E-02,
    2.519E-02,
    2.226E-02,
    1.946E-02,
    1.682E-02,
    1.437E-02,
    1.215E-02,
    1.016E-02,
    8.400E-03,
    6.873E-03,
    5.564E-03,
    4.457E-03,
    3.533E-03,
    2.772E-03,
    2.154E-03,
    1.656E-03,
    1.261E-03,
    9.513E-04,
    7.107E-04,
    5.259E-04,
    3.856E-04,
    2.801E-04,
    2.017E-04,
    1.439E-04,
    1.017E-04,
    7.126E-05,
    4.948E-05,
    3.405E-05,
    2.322E-05,
    1.570E-05,
    5.005E-06
  };
}
