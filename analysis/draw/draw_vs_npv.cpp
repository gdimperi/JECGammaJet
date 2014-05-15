#include <stdlib.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"

#include "etaBinning.h"
#include "ptBinning.h"
#include "vertexBinning.h"


bool useMCassoc_ = false;
bool ONEVTX = false;
bool OUTPUT_GRAPHS = true;

int main(int argc, char* argv[]) {

  if (argc != 7 && argc != 8) {
    std::cout << "USAGE: ./draw_vs_npv [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')] [flags=\"\"]" << std::endl;
    exit(23);
  }

  gROOT->SetBatch();

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

  //jetAlgo = (jetAlgo == "ak5") ? "AK5" : "AK7";
  if (jetAlgo == "ak5") jetAlgo = "AK5" ;
  if (jetAlgo == "ak7") jetAlgo = "AK7" ;
  if (jetAlgo == "ca8") jetAlgo = "CA8" ;
  recoType = (recoType == "pf") ? "PFlow" : "Calo";
  std::string postFix = recoType + jetAlgo;

  postFix += "chs";

  //giulia debug
  //postFix = "PFlowCA8";

  //Float_t etamax = 3.;
  //bool sameEvents = false; //until njets histos have no overflows... or maybe use GetEntries instead of integral?

  drawBase* db = new drawBase("PhotonJet", recoType, jetAlgo, OUTPUT_GRAPHS);
  db->set_pdf_aussi((bool)false);
  db->set_flags(flags);
  db->set_isCMSArticle((bool)true);

  std::cout << "flags set." << std::endl;

  TString dataFileName;
  TFile* dataFile;

  if(data_dataset!=""){
    if (flags.length() > 0) {
      
      dataFileName = TString::Format("PhotonJet_%s_%s_%s.root", data_dataset.c_str(), postFix.c_str(), flags.c_str());
      //dataFileName = TString::Format("PhotonJet_%s_%s_%s.root", data_dataset.c_str(), "PFlowCA8", flags.c_str());
    } else {
      dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
      //dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), "PFlowCA8");
    }
    dataFile = TFile::Open(dataFileName);
  
    if (dataFile) {
      std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
      db->add_dataFile(dataFile, data_dataset);
    }
  }

  TString mc1FileName;
  
  if (flags.length() > 0) {
    mc1FileName = TString::Format("PhotonJet_%s_%s_%s.root", mc_photonjet.c_str(), postFix.c_str(), flags.c_str());
  } else {
    mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str());
  }

  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  //giulia debug  
  //TFile* mcPhotonJetFile = TFile::Open("../../bin/PhotonJet_Gj-Pt120to170_PFlowAK5chs.root");


  std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;

  if (mcPhotonJetFile) {
    db->add_mcFile(mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);
  }

  if (mc_QCD != "") {
    TString mc2FileName;
    if (flags.length() > 0) {
      mc2FileName = TString::Format("PhotonJet_%s_%s_%s.root", mc_QCD.c_str(), postFix.c_str(), flags.c_str());
    } else {
      mc2FileName = TString::Format("PhotonJet_%s_%s.root", mc_QCD.c_str(), postFix.c_str());
    }
    TFile* mcQCDFile = TFile::Open(mc2FileName);
    std::cout << "Opened mc file '" << mc2FileName << "'." << std::endl;

    if (mcQCDFile && mc_QCD != mc_photonjet) {
      db->add_mcFile(mcQCDFile, mc_QCD, "QCD MC", 38);
    }
  }

  // MC should already be normalized to a lumi of 1 pb-1
  // Read luminosity
  double dLumi = 1e6;
  if (dataFile) {
    TParameter<double>* lumi = static_cast<TParameter<double>*>(dataFile->Get("analysis/luminosity"));
    dLumi = lumi->GetVal();
  }

  db->set_lumi(dLumi * 1e-6);
  //db->set_lumi(3000);
  if (norm == "LUMI") {
    db->set_lumiNormalization();
  } else {
    db->set_shapeNormalization();
  }

  db->setFolder("analysis");
  std::string outputDir = "PhotonJetPlots_" + db->get_fullSuffix() + "/vs_npv";
  db->set_outputdir(outputDir);

  bool log = true;
  gErrorIgnoreLevel = kWarning;

  db->setOutputGraphs(OUTPUT_GRAPHS);

  VertexBinning vertexBinning;
  std::vector<std::pair<int, int> > vertexBins = vertexBinning.getBinning();

  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();

  db->set_rebin(2);

  // Balancing
  db->setFolder("analysis/vertex");
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("resp_balancing_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_vertex(vertexBins, responseName.Data(), "Balancing Response", "", "Events", log);

    // Raw jets
    //responseName = TString::Format("resp_balancing_raw_%s", etaBinning.getBinName(i).c_str());
    //db->drawHisto_vs_vertex(ptBins, responseName.Data(), "Balancing Response (raw jets)", "", "Events", log);

    //giulia true response
    responseName = TString::Format("true_resp_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_vertex(vertexBins, responseName.Data(), "True Response", "", "Events", log);
    responseName = TString::Format("true_resp_raw_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_vertex(vertexBins, responseName.Data(), "True Response (raw jets)", "", "Events", log);

  }
  // Special case eta < 1.3

  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_vertex(vertexBins, "resp_balancing_eta013", "Balancing Response", "", "Events", log);
  //db->drawHisto_vs_pt(ptBins, "resp_balancing_raw_eta013", "Balancing Response (raw jets)", "", "Events", log);

    //giulia true response
  db->drawHisto_vs_vertex(vertexBins, "true_resp_eta013", "True Response", "", "Events", log);
  db->drawHisto_vs_vertex(vertexBins, "true_resp_raw_eta013", "True Response (raw jets)", "", "Events", log);

  // MPF
  //db->setFolder("analysis/vertex");
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("resp_mpf_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_vertex(vertexBins, responseName.Data(), "MPF Response", "", "Events", log);

    // Raw jets
    //responseName = TString::Format("resp_mpf_raw_%s", etaBinning.getBinName(i).c_str());
    //db->drawHisto_vs_pt(ptBins, responseName.Data(), "MPF Response (raw ME_{T})", "", "Events", log);

  }
  // Special case eta < 1.3

  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_vertex(vertexBins, "resp_mpf_eta013", "MPF Response", "", "Events", log);
  //db->drawHisto_vs_pt(ptBins, "resp_mpf_raw_eta013", "MPF Response (raw ME_{T})", "", "Events", log);

  delete db;
  db = NULL;

  return 0;

}


