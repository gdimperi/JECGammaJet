#include <stdlib.h>
#include "drawExtrap.h"
#include "fitTools.h"


bool useMCassoc_ = false;
bool NOQ = false;
bool FIXM = false;
bool EXCLUDE_FIRST_POINT = false;
bool OUTPUT_GRAPHS = true;

int main(int argc, char* argv[]) {

  if (argc != 8 && argc != 9 && argc != 10) {
    std::cout << "USAGE: ./drawPhotonJetExtrap [data_dataset] [mc_dataset] [mc2_dataset] [recoType] [jetAlgo] [FIT_RMS] [flags=\"\"] [GENJETS=\"\"]" << std::endl;
    exit(23);
  }

  std::string data_dataset(argv[1]);
  std::string mc_dataset(argv[2]);
  std::string mc2_dataset(argv[3]);
  std::string recoType(argv[4]);
  std::string jetAlgo(argv[5]);
  std::string FIT_RMS(argv[6]);
  std::string flags = "";
  if (argc == 8) {
    std::string flags_str(argv[7]);
    flags = flags_str;
    std::cout << "flags set." << std::endl;
  }
  bool useGenJets = false;
  if (argc == 9) {
    std::string flags_str(argv[8]);
    if (flags_str == "GENJETS") useGenJets = true;
  }


  std::string algoType;
  if (recoType == "calo")
    algoType = jetAlgo;
  else
    algoType = recoType + jetAlgo;
  if (recoType == "jpt" && jetAlgo == "akt5") algoType = "jptak5";

  jetAlgo = (jetAlgo == "ak5") ? "AK5" : "AK7";
  recoType = (recoType == "pf") ? "PFlow" : "Calo";
  std::string postFix = jetAlgo + recoType;

  std::string mcFlags = (useGenJets) ? "GENJETS_" + flags : flags;

  drawExtrap* db = new drawExtrap("PhotonJet", recoType, jetAlgo, true, mcFlags);
  db->set_pdf_aussi((bool)false);
  db->set_isCMSArticle((bool)false);

  db->set_FITRMS(FIT_RMS);

  std::string NOQtext = (NOQ) ? "_NOQ" : "";
  std::string FIXMtext = (FIXM) ? "_FIXM" : "";
  std::string NOFIRSTPtext = (EXCLUDE_FIRST_POINT) ? "_NOFIRSTP" : "";

  char outputdir_char[200];
  if (mcFlags != "") {
    sprintf(outputdir_char, "PhotonJetExtrapPlots_%s_vs_%s_%s_%s_%s%s%s%s", data_dataset.c_str(), mc_dataset.c_str(), algoType.c_str(), mcFlags.c_str(), FIT_RMS.c_str(), FIXMtext.c_str(), NOFIRSTPtext.c_str(), NOQtext.c_str());
  } else {
    sprintf(outputdir_char, "PhotonJetExtrapPlots_%s_vs_%s_%s_%s%s%s%s", data_dataset.c_str(), mc_dataset.c_str(), algoType.c_str(), FIT_RMS.c_str(), FIXMtext.c_str(), NOFIRSTPtext.c_str(), NOQtext.c_str());
  }
  std::string outputdir_str(outputdir_char);

  //std::vector< float > ptPhot_binning = fitTools::getPtPhot_binning();

  db->set_outputdir(outputdir_str);

  TString dataFileName;
  if (flags.length() > 0) {
    dataFileName = TString::Format("PhotonJet_%s_%s_%s.root", data_dataset.c_str(), postFix.c_str(), flags.c_str());
  } else {
    dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
  }

  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  db->add_dataFile(dataFile, data_dataset);

  TString mc1FileName;
  if (flags.length() > 0) {
    mc1FileName = TString::Format("PhotonJet_%s_%s_%s.root", mc_dataset.c_str(), postFix.c_str(), flags.c_str());
  } else {
    mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_dataset.c_str(), postFix.c_str());
  }
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;

  db->add_mcFile(mcPhotonJetFile, mc_dataset, "#gamma+jet MC", 46);

  TString mc2FileName;
  if (flags.length() > 0) {
    mc2FileName = TString::Format("PhotonJet_%s_%s_%s.root", mc2_dataset.c_str(), postFix.c_str(), flags.c_str());
  } else {
    mc2FileName = TString::Format("PhotonJet_%s_%s.root", mc2_dataset.c_str(), postFix.c_str());
  }
  TFile* mcQCDFile = TFile::Open(mc2FileName);
  std::cout << "Opened mc file '" << mc2FileName << "'." << std::endl;

  if (mc_dataset != mc2_dataset) {
    db->add_mcFile(mcQCDFile, mc2_dataset, "QCD MC", 38);
  } 

  db->set_lumiNormalization(5.09868 * 1e3);

  db->set_NOQ(NOQ);
  db->set_FIXM(FIXM);
  db->set_EXCLUDEFIRSTPOINT(EXCLUDE_FIRST_POINT);

  //db->drawResponseExtrap("eta011", false, "RecoRel", true);
  db->drawResponseExtrap("eta011", false);
  db->drawResponseExtrap("eta011", true);
  //db->drawResponseExtrap("eta011", (bool)true, "RecoRelRaw");

  //db->drawResponseExtrap("eta013", false, "RecoRel", true);
  db->drawResponseExtrap("eta013", false);
  db->drawResponseExtrap("eta013", true);
  //db->drawResponseExtrap("eta013", (bool)true, "RecoRelRaw");

  //db->drawResponseExtrap("eta1524", false, "RecoRel", true);
  db->drawResponseExtrap("eta1524", false);
  db->drawResponseExtrap("eta1524", true);

  //db->drawResponseExtrap("eta243", false, "RecoRel", true);
  db->drawResponseExtrap("eta243", false);
  db->drawResponseExtrap("eta243", true);

  //  db->drawResponseExtrap("eta35");

  delete db;

  return 0;
}


