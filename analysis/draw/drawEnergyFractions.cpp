#include <stdlib.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"

#include "etaBinning.h"
#include "ptBinning.h"

#include <TColor.h>
#include <TStyle.h>

int main(int argc, char* argv[]) {

//TStyle  *st_style = new TStyle("trystandardstyle", "");
//st_style->cd();
gROOT->SetStyle("Plain");
gROOT->ForceStyle();

  if (argc != 7 && argc != 8) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')] [flags=\"\"]" << std::endl;
    exit(23);
  }
  std::string data_dataset(argv[1]);
  std::string mc_photonjet(argv[2]);
  std::string mc_QCD(argv[3]);
  std::string recoType(argv[4]);
  std::string jetAlgo(argv[5]);
  std::string norm(argv[6]);
  if (norm != "LUMI" && norm != "SHAPE") {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }
  std::string flags = "";
  if (argc == 8) {
    std::string flags_str(argv[7]);
    flags = flags_str;
  }


  std::string algoType;
  if (recoType == "calo") {
    algoType = jetAlgo;
  } else {
    algoType = recoType + jetAlgo;
  }
  if (recoType == "jpt" && jetAlgo == "akt5") {
    algoType = "jptak5";
  }

  jetAlgo = (jetAlgo == "ak5") ? "AK5" : "AK7";
  recoType = (recoType == "pf") ? "PFlow" : "Calo";
  std::string postFix = recoType + jetAlgo;

  postFix += "chs";

  TString dataFileName;
  if (flags.length() > 0) {
    dataFileName = TString::Format("PhotonJet_%s_%s_%s.root", data_dataset.c_str(), postFix.c_str(), flags.c_str());
  } else {
    dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
  }

  TFile* dataFile = TFile::Open(dataFileName);

  if (dataFile) {
    std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
  }

  TString mc1FileName;
  if (flags.length() > 0) {
    mc1FileName = TString::Format("PhotonJet_%s_%s_%s.root", mc_photonjet.c_str(), postFix.c_str(), flags.c_str());
  } else {
    mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str());
  }
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;

//ok that's not VERY nice, but i'm tired and want these ploooots!
 TFile * inputfile = mcPhotonJetFile; 
 TString type("mc");
// TFile * inputfile = dataFile;
// TString type("data");

  PtBinning ptBinning;
  std::vector<std::pair<float, float> > ptBins = ptBinning.getBinning();
  size_t ptBinningSize = ptBinning.size();
  std::pair<float, float> currentBin;

  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();
  TString histoName;
  TString histoName_tmp;
  TH1F *h_tmp;
  TString stackName;
//
  TCanvas *c = new TCanvas("c","c",600,600);

  Float_t meanTruncFraction = 0.99;
  Float_t rmsTruncFraction = 0.99;  
  Float_t dataResponse = 0.;
  Float_t dataResponseErr = 0.;
  Float_t dataRMS = 0.;
  Float_t dataRMSErr = 0.;

  Float_t ptphot_bins[ptBinningSize+1];

  TH1F *hCHFrac[etaBinningSize];
  TH1F *hNHFrac[etaBinningSize];
  TH1F *hElFrac[etaBinningSize];
  TH1F *hPhFrac[etaBinningSize];
  TH1F *hMuFrac[etaBinningSize];
  THStack *hsum[etaBinningSize];

   for(size_t j = 0; j<ptBinningSize; j++) {
     currentBin = ptBinning.getBinValue(j);
     ptphot_bins[j] = currentBin.first;
     if(j == ptBinningSize-1) ptphot_bins[j+1] = currentBin.second;
   }

  for (size_t i = 0; i < etaBinningSize; i++) {
   for(int frac=0; frac<5; frac++) {
if(frac==0){
//    histoName = TString::Format("ChHadronFraction_%s", etaBinning.getBinName(i).c_str());
    histoName = TString::Format("ChHadron_realFractionRaw_%s", etaBinning.getBinName(i).c_str());
    hCHFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
} else if(frac==1){
//    histoName = TString::Format("NHadronFraction_%s", etaBinning.getBinName(i).c_str());
    histoName = TString::Format("NHadron_realFractionRaw_%s", etaBinning.getBinName(i).c_str());
    hNHFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
} else if(frac==2) {
//    histoName = TString::Format("ElFraction_%s", etaBinning.getBinName(i).c_str());
    histoName = TString::Format("El_realFractionRaw_%s", etaBinning.getBinName(i).c_str());
    hElFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
} else if(frac==3) {
//    histoName = TString::Format("PhFraction_%s", etaBinning.getBinName(i).c_str());
    histoName = TString::Format("Ph_realFractionRaw_%s", etaBinning.getBinName(i).c_str());
    hPhFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
} else {
//    histoName = TString::Format("MuFraction_%s", etaBinning.getBinName(i).c_str());
    histoName = TString::Format("Mu_realFractionRaw_%s", etaBinning.getBinName(i).c_str());
    hMuFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
}

   for(size_t j = 0; j<ptBinningSize; j++) {
    histoName_tmp=histoName;
//     histoName.Append(TString::Format("_ptPhot_%s", ptBinning.getBinName(i).c_str()));
//     histoName.Append(TString::Format("_%s", ptBinning.getBinName(i+1).c_str()));
     currentBin = ptBinning.getBinValue(j);
     histoName_tmp.Append(TString::Format("_ptPhot_%i",int(currentBin.first)));
     histoName_tmp.Append(TString::Format("_%i", int(currentBin.second)));
     h_tmp=(TH1F*)inputfile->Get("analysis/ecomposition/"+histoName_tmp);
//     fitTools::getTruncatedMeanAndRMS(h_tmp, dataResponse, dataResponseErr, dataRMS, dataRMSErr, meanTruncFraction, rmsTruncFraction);
dataResponse = h_tmp->GetMean();
dataResponseErr = h_tmp->GetMeanError();
dataRMS = h_tmp->GetRMS();
dataRMSErr = h_tmp->GetRMSError();

if(frac==0) {
     hCHFrac[i]->SetBinContent(j+1,dataResponse);
     hCHFrac[i]->SetBinError(j+1,dataResponseErr);
} else if(frac==1) {
     hNHFrac[i]->SetBinContent(j+1,dataResponse);
     hNHFrac[i]->SetBinError(j+1,dataResponseErr);
} else if(frac==2) {
     hElFrac[i]->SetBinContent(j+1,dataResponse);
     hElFrac[i]->SetBinError(j+1,dataResponseErr);
} else if(frac==3){
     hPhFrac[i]->SetBinContent(j+1,dataResponse);
     hPhFrac[i]->SetBinError(j+1,dataResponseErr);
} else {
     hMuFrac[i]->SetBinContent(j+1,dataResponse);
     hMuFrac[i]->SetBinError(j+1,dataResponseErr);
}

     }
   }


gROOT->SetStyle("Plain");
gROOT->ForceStyle();

//stack them
  stackName = TString::Format("JetEnergyComposition_%s", etaBinning.getBinName(i).c_str());

  hsum[i] = new THStack(stackName,stackName);

  hsum[i]->UseCurrentStyle();
  hCHFrac[i]->UseCurrentStyle();
  hNHFrac[i]->UseCurrentStyle();
  hElFrac[i]->UseCurrentStyle();
  hPhFrac[i]->UseCurrentStyle();
  hMuFrac[i]->UseCurrentStyle();

  hCHFrac[i]->SetFillColor(2);
  hCHFrac[i]->SetFillStyle(1001);
  hCHFrac[i]->SetMarkerStyle(20);
  hNHFrac[i]->SetFillColor(8);
  hNHFrac[i]->SetFillStyle(1001);
  hNHFrac[i]->SetMarkerStyle(21);
  hElFrac[i]->SetFillColor(4);
  hElFrac[i]->SetFillStyle(1001);
  hElFrac[i]->SetMarkerStyle(22);
  hPhFrac[i]->SetFillColor(6);
  hPhFrac[i]->SetFillStyle(1001);
  hPhFrac[i]->SetMarkerStyle(23);
  hMuFrac[i]->SetFillColor(7);
  hMuFrac[i]->SetFillStyle(1001);
  hMuFrac[i]->SetMarkerStyle(29);

  hsum[i]->Add(hMuFrac[i]);
  hsum[i]->Add(hElFrac[i]);
  hsum[i]->Add(hNHFrac[i]);
  hsum[i]->Add(hPhFrac[i]);
  hsum[i]->Add(hCHFrac[i]);




/*  c->cd();
  hsum[i]->Draw();
  c->SaveAs(stackName+".eps");
*/
}



TFile out(type+"testEComposition.root","recreate");
out.cd();
  for (size_t i = 0; i < etaBinningSize; i++) {
  hCHFrac[i]->Write();
  hNHFrac[i]->Write();
  hElFrac[i]->Write();
  hPhFrac[i]->Write();
  hMuFrac[i]->Write();
  hsum[i]->Write();
  }
out.Close();

 }


