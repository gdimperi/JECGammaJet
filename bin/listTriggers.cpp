#include <TFile.h>

#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTree.h>
#include <TParameter.h>
#include <TH2D.h>

#include <fstream>
#include <sstream>

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <FWCore/FWLite/interface/AutoLibraryLoader.h>
#include <FWCore/Framework/interface/Event.h>

#include <DataFormats/Common/interface/Handle.h>
#include <DataFormats/FWLite/interface/Event.h>
#include <DataFormats/FWLite/interface/ChainEvent.h>
#include <DataFormats/Math/interface/deltaPhi.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Photon.h>

#include <PhysicsTools/FWLite/interface/TFileService.h>

#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include "tclap/CmdLine.h"

#include "gammaJetFinalizer.h"
#include "PUReweighter.h"
#include "JECReader.h"

#include <boost/regex.hpp>

#define RESET_COLOR "\033[m"
#define MAKE_RED "\033[31m"
#define MAKE_BLUE "\033[34m"

#define ADD_TREES false

#define DELTAPHI_CUT (M_PI - 1)

#define TRIGGER_OK                    0
#define TRIGGER_NOT_FOUND            -1
#define TRIGGER_FOUND_BUT_PT_OUT     -2

bool EXIT = false;

AnalysisTree analysis;
std::vector<std::string> mInputFiles;

std::map<std::string, int64_t> mTriggers;

void loadFiles(TChain& chain) {
  for (std::vector<std::string>::const_iterator it = mInputFiles.begin(); it != mInputFiles.end(); ++it) {
    chain.Add(it->c_str());
  }
}

std::string removeTriggerVersion(const std::string& trigger) {
  static boost::regex r(R"(_v[0-9]*)");
  return boost::regex_replace(trigger, r, "");
}

void listTriggers() {

  std::cout << "Opening files ..." << std::endl;

  TChain analysisChain("gammaJet/analysis");
/*  TChain photonChain("gammaJet/photon");
  TChain genPhotonChain("gammaJet/photon_gen");
  TChain muonsChain("gammaJet/muons");
  TChain electronsChain("gammaJet/electrons");

  TString treeName = TString::Format("gammaJet/%s/first_jet", postFix.c_str());
  TChain firstJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/first_jet_raw", postFix.c_str());
  TChain firstRawJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/first_jet_gen", postFix.c_str());
  TChain firstGenJetChain(treeName);

  treeName = TString::Format("gammaJet/%s/second_jet", postFix.c_str());
  TChain secondJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/second_jet_raw", postFix.c_str());
  TChain secondRawJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/second_jet_gen", postFix.c_str());
  TChain secondGenJetChain(treeName);

  treeName = TString::Format("gammaJet/%s/met", postFix.c_str());
  TChain metChain(treeName);
  treeName = TString::Format("gammaJet/%s/met_raw", postFix.c_str());
  TChain rawMetChain(treeName);
  treeName = TString::Format("gammaJet/%s/met_gen", postFix.c_str());
  TChain genMetChain(treeName);

  treeName = TString::Format("gammaJet/%s/misc", postFix.c_str());
  TChain miscChain(treeName);
*/

  loadFiles(analysisChain);
/*  loadFiles(photonChain);
  loadFiles(genPhotonChain);
  loadFiles(muonsChain);
  loadFiles(electronsChain);

  loadFiles(firstJetChain);
  loadFiles(firstRawJetChain);
  loadFiles(firstGenJetChain);

  loadFiles(secondJetChain);
  loadFiles(secondGenJetChain);
  loadFiles(secondRawJetChain);

  loadFiles(metChain);
  loadFiles(genMetChain);
  loadFiles(rawMetChain);

  loadFiles(miscChain);
*/

  analysis.Init(&analysisChain);
/*  photon.Init(&photonChain);
  genPhoton.Init(&genPhotonChain);
  muons.Init(&muonsChain);
  electrons.Init(&electronsChain);

  firstJet.Init(&firstJetChain);
  firstRawJet.Init(&firstRawJetChain);
  firstGenJet.Init(&firstGenJetChain);

  secondJet.Init(&secondJetChain);
  secondRawJet.Init(&secondRawJetChain);
  secondGenJet.Init(&secondGenJetChain);

  MET.Init(&metChain);
  genMET.Init(&genMetChain);
  rawMET.Init(&rawMetChain);

  misc.Init(&miscChain);
*/

  std::cout << "done." << std::endl;

  std::cout << "Processing..." << std::endl;

  uint64_t totalEvents = analysisChain.GetEntries();
  uint64_t from = 0;
  uint64_t to = totalEvents;

  for (uint64_t i = from; i < to; i++) {

    if ((i - from) % 50000 == 0) {
      std::cout << "Processing event #" << (i - from + 1) << " of " << (to - from) << " (" << (float) (i - from) / (to - from) * 100 << "%)" << std::endl;
    }

    if (EXIT) {
      break;
    }

    analysis.GetEntry(i);
/*    photon.GetEntry(i);
    genPhoton.GetEntry(i);
    muons.GetEntry(i);
    electrons.GetEntry(i);

    firstJet.GetEntry(i);
    firstRawJet.GetEntry(i);
    firstGenJet.GetEntry(i);

    secondJet.GetEntry(i);
    secondRawJet.GetEntry(i);
    secondGenJet.GetEntry(i);

    MET.GetEntry(i);
    genMET.GetEntry(i);
    rawMET.GetEntry(i);

    misc.GetEntry(i);
*/
    size_t size = analysis.trigger_names->size();
    for (int i = size - 1; i >= 0; i--) {
      bool passed = analysis.trigger_results->at(i);
      if (! passed)
        continue;

      /*TString foo = analysis.trigger_names->at(i);
      if (! foo.Contains("Photon"))
        continue;*/

      std::string t = removeTriggerVersion(analysis.trigger_names->at(i));

      mTriggers[t]++;
    }
  }

  std::cout << std::endl << std::endl;

  std::cout << "Trigger summary" << std::endl;
  std::cout << "===============" << std::endl << std::endl;

  for (auto i: mTriggers) {
    std::cout << std::left << std::setw (85) << i.first << i.second << std::endl;
  }
}


void GammaJetFinalizer::checkInputFiles() {
  for (std::vector<std::string>::iterator it = mInputFiles.begin(); it != mInputFiles.end();) {
    TFile* f = TFile::Open(it->c_str());
    if (! f) {
      std::cerr << "Error: can't open '" << it->c_str() << "'. Removed from input files." << std::endl;
      it = mInputFiles.erase(it);
      continue;
    }

    TTree* analysis = static_cast<TTree*>(f->Get("gammaJet/analysis"));
    if (! analysis || analysis->GetEntry(0) == 0) {
      std::cerr << "Error: Trees inside '" << it->c_str() << "' were empty. Removed from input files." << std::endl;
      it = mInputFiles.erase(it);

      f->Close();
      delete f;

      continue;
    }

    f->Close();
    delete f;

    ++it;
  }
}

std::vector<std::string> readInputFiles(const std::string& list) {
  std::ifstream f(list.c_str());
  std::string line;
  std::vector<std::string> files;
  while (std::getline(f, line)) {
    boost::algorithm::trim(line);
    if (line.length() == 0 || line[0] == '#')
      continue;

    files.push_back(line);
  }

  if (files.size() == 0) {
    throw new TCLAP::ArgException("No input files found in " + list);
  }

  return files;
}

void handleCtrlC(int s){
  EXIT = true;
}

int main(int argc, char** argv) {
  struct sigaction sigIntHandler;

  sigIntHandler.sa_handler = handleCtrlC;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;

  sigaction(SIGINT, &sigIntHandler, NULL);

  try {
    TCLAP::CmdLine cmd("List all triggers passed in analysis trees", ' ', "0.1");

    TCLAP::MultiArg<std::string> inputArg("i", "in", "Input file", true, "string");
    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "Text file containing input files", true, "input.list", "string");
    cmd.xorAdd(inputArg, inputListArg);

    cmd.parse(argc, argv);

    std::vector<std::string> files;
    if (inputArg.isSet()) {
      files = inputArg.getValue();
    } else {
      files = readInputFiles(inputListArg.getValue());
    }

    mInputFiles = files;
    listTriggers();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }
}
