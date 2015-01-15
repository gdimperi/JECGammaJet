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

    
    std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;

    int NBins = puHisto->GetNbinsX();

    for (int ibin = 1; ibin < NBins + 1; ++ibin) {
      std::cout << "   " << ibin - 1 << " " << puHisto->GetBinContent(ibin) << std::endl;
    }
    

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
  mPUCoefs[PUProfile::RDAB] = {
                  5.99688E-12,
                  3.9153E-10,
                  3.71925E-07,
                  3.03413E-05,
                  7.34144E-05,
                  0.000347794,
                  0.001975866,
                  0.00517831,
                  0.010558132,
                  0.018312126,
                  0.028852472,
                  0.040969532,
                  0.050957346,
                  0.058116316,
                  0.063832662,
                  0.067924013,
                  0.068664032,
                  0.06634696,
                  0.06258012,
                  0.058673398,
                  0.055389375,
                  0.05287113,
                  0.050854502,
                  0.048745749,
                  0.045577918,
                  0.040642954,
                  0.033995191,
                  0.026305836,
                  0.018543156,
                  0.011719526,
                  0.006558513,
                  0.003227318,
                  0.00139409,
                  0.000529933,
                  0.000178348,
                  5.36943E-05,
                  1.46867E-05,
                  3.72253E-06,
                  8.91697E-07,
                  2.04492E-07,
                  4.49872E-08,
                  9.4322E-09,
                  1.86873E-09,
                  3.48142E-10,
                  6.10474E-11,
                  1.01401E-11,
                  1.61054E-12,
                  2.4678E-13,
                  3.66567E-14,
                  5.27305E-15,
                  7.30257E-16,
                  9.65776E-17,
                  1.20959E-17,
                  1.41538E-18,
                  1.59156E-19,
                  1.00245E-20,
                  0,
                  0,
                  0,
                  0
  };
  mPUCoefs[PUProfile::RDC] = {
                      2.30103E-06,
                      6.10952E-06,
                      1.20637E-05,
                      2.55143E-05,
                      3.22897E-05,
                      0.000120076,
                      0.000800727,
                      0.003176761,
                      0.008762044,
                      0.018318462,
                      0.032467397,
                      0.050617606,
                      0.064967885,
                      0.072194727,
                      0.076272353,
                      0.079908524,
                      0.08275316,
                      0.083946479,
                      0.083972048,
                      0.082763689,
                      0.080056866,
                      0.076035269,
                      0.071231913,
                      0.065785365,
                      0.059393488,
                      0.052004578,
                      0.044060277,
                      0.036161288,
                      0.028712778,
                      0.021943183,
                      0.016047885,
                      0.011184911,
                      0.007401996,
                      0.004628635,
                      0.002717483,
                      0.001488303,
                      0.000757141,
                      0.000357766,
                      0.000157852,
                      6.5729E-05,
                      2.62131E-05,
                      1.01753E-05,
                      3.89905E-06,
                      1.48845E-06,
                      5.67959E-07,
                      2.16348E-07,
                      8.19791E-08,
                      3.07892E-08,
                      1.1436E-08,
                      4.19872E-09,
                      1.52513E-09,
                      5.48804E-10,
                      1.95797E-10,
                      6.92424E-11,
                      2.42383E-11,
                      8.37978E-12,
                      2.85421E-12,
                      9.55583E-13,
                      3.13918E-13,
                      1.01066E-13
  };
  mPUCoefs[PUProfile::RDD] = {
                2.54656e-11,
                6.70051e-11,
                2.74201e-06,
                6.9111e-06,
                5.00919e-06,
                6.24538e-05,
                0.000338679,
                0.000892795,
                0.00237358,
                0.00686023,
                0.0144954,
                0.026012,
                0.0360377,
                0.0420151,
                0.0457901,
                0.0482319,
                0.0503176,
                0.052569,
                0.0546253,
                0.0561205,
                0.0568903,
                0.0570889,
                0.0566598,
                0.0553747,
                0.0531916,
                0.0501454,
                0.0463101,
                0.0417466,
                0.0364842,
                0.0306443,
                0.0245417,
                0.0186276,
                0.0133446,
                0.00900314,
                0.00571947,
                0.00342706,
                0.00194292,
                0.00104671,
                0.000538823,
                0.000266973,
                0.000128572,
                6.09778e-05,
                2.89549e-05,
                1.40233e-05,
                7.04619e-06,
                3.71289e-06,
                2.055e-06,
                1.18713e-06,
                7.08603e-07,
                4.32721e-07,
                2.6817e-07,
                1.67619e-07,
                1.05157e-07,
                6.59446e-08,
                4.11915e-08,
                2.55494e-08
  };

}
