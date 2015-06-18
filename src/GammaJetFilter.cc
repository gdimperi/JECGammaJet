//
// Package:    GammaJetFilter
// Class:      GammaJetFilter
// 
/**\class GammaJetFilter GammaJetFilter.cc JetMETCorrections/GammaJetFilter/src/GammaJetFilter.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Sébastien Brochet
//         Created:  Thu Mar 15 11:27:48 CET 2012
// $Id$
//
//


// system include files
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <unordered_map>
#include <memory>
#include <string>
#include "Math/GenVector/LorentzVector.h"

// Boost
#include "boost/shared_ptr.hpp"

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/GammaJetFilter/interface/json/json.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "JetMETCorrections/GammaJetFilter/interface/EnergyScaleCorrection_class.h"

#include <TParameter.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <boost/regex.hpp>

//
// class declaration
//

enum JetAlgorithm {
  AK4,
  AK8
};

struct JetInfos {
  JetAlgorithm algo;
  edm::InputTag inputTag;
};

#define FOREACH(x) for (std::vector<std::string>::const_iterator it = x.begin(); it != x.end(); ++it)

class GammaJetFilter : public edm::EDFilter {
  public:
    explicit GammaJetFilter(const edm::ParameterSet&);
    ~GammaJetFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    virtual bool beginRun(edm::Run&, edm::EventSetup const&);
    virtual bool endRun(edm::Run&, edm::EventSetup const&);
    virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    //giulia--- correct photon not necessary until we have data @13 TeV 
    void correctPhoton(pat::Photon& photon, edm::Event& iEvent, int isData, int nPV);
    void correctJets(pat::JetCollection& jets, edm::Event& iEvent, const edm::EventSetup& iSetup);
    void extractRawJets(pat::JetCollection& jets);
    //giulia --- comment qg tagging stuff
    void processJets(pat::Photon* photon, pat::JetCollection& jets, const JetAlgorithm algo, /* edm::Handle<edm::ValueMap<float>>& qgTagMLP, edm::Handle<edm::ValueMap<float>>& qgTagLikelihood,*/ const edm::Handle<pat::JetCollection>& handleForRef, std::vector<TTree*>& trees);

    void correctMETWithTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, edm::Event& event);
   void correctMETWithRegressionAndTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets,  edm::Event& event, pat::Photon& photon, const pat::PhotonRef& photonRef);
   //giulia --- comment footprint stuff
   //void correctMETWithFootprintAndTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets,  edm::Event& event, pat::Photon& photon, const pat::PhotonRef& photonRef);
//(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, edm::Event& event,const pat::PhotonRef& photonRef, float regressionCorr);

    bool isValidPhotonEB(const pat::Photon& photon, const double rho, const EcalRecHitCollection* recHits, const CaloTopology& topology);
    bool isValidPhotonEB2012(const pat::PhotonRef& photonRef, edm::Event& event);
    bool isValidPhotonEB_PHYS14(const pat::PhotonRef& photonRef, edm::Event& event);
    bool isValidJet(const pat::Jet& jet);

    void readJSONFile();
    void readCSVFile();
    void updateLuminosity(const edm::LuminosityBlock& lumiBlock);

    // ----------member data ---------------------------
    bool mIsMC;
    bool mFilterData;
    std::string mJSONFile;
    std::string mCSVFile;
    boost::shared_ptr<Json::Value> mValidRuns;
    boost::shared_ptr<Json::Value> mCurrentRunValidLumis;
    std::map<std::pair<unsigned int, unsigned int>, double> mLumiByLS;
    std::map<std::pair<unsigned int, unsigned int>, double> mTruePUByLS;
    bool mIsValidLumiBlock;
    double mCurrentTruePU;

    // Photon ID
    PFIsolationEstimator mPFIsolator;
    bool mCorrPhotonWRegression;

    bool mRedoTypeI;
    bool mDoFootprint;
    bool mDoJEC;
    bool mJECFromRaw;
    std::string mCorrectorLabel;
    GreaterByPt<pat::Jet> mSorter;

    bool mFirstJetPtCut;
    double mFirstJetThreshold;

    std::vector<std::string> mJetCollections;
    std::map<std::string, JetInfos> mJetCollectionsData;

    // Input Tags
    edm::InputTag mPhotonsIT;
    edm::InputTag mJetsAK4PFlowIT;
    edm::InputTag mJetsAK8PFlowIT;
    edm::InputTag mJetsAK4CaloIT;
    edm::InputTag mJetsAK8CaloIT;

    boost::shared_ptr<JetIDSelectionFunctor> mCaloJetID;
    pat::strbitset mCaloJetIDRet;

    double mPtHatMin;
    double mPtHatMax;
    //double sum_GenWeights;

    // Trees
    void createTrees(const std::string& rootName, TFileService& fs);
    TTree* mGenParticlesTree;
    TTree* mPhotonTree;
    TTree* mPhotonGenTree;
    TTree* mAnalysisTree;
    TTree* mElectronsTree;
    TTree* mMuonsTree;
    TParameter<double>*    mTotalLuminosity;
    float                  mEventsWeight;
    TParameter<long long>* mProcessedEvents;
    TParameter<long long>* mSelectedEvents;

    std::map<std::string, std::vector<TTree*> > mJetTrees;
    std::map<std::string, std::vector<TTree*> > mMETTrees;
    std::map<std::string, TTree*>               mMiscTrees;
    std::map<std::string, TTree*> mMETNFTrees;

    // TParameters for storing current config (JEC, correctorLabel, Treshold, etc...
    TParameter<bool>*             mJECRedone;
    TParameter<bool>*             mJECFromRawParameter;
    TNamed*                       mJECCorrectorLabel;
    TParameter<bool>*             mFirstJetPtCutParameter;
    TParameter<double>*           mFirstJetThresholdParameter;

    // to store sum of event weights
    TH1F* h_sumW;
    
    // DEBUG
    TH1F* mFirstJetPhotonDeltaPhi;
    TH1F* mFirstJetPhotonDeltaR;
    TH1F* mFirstJetPhotonDeltaPt;
    TH2F* mFirstJetPhotonDeltaPhiDeltaR;

    TH1F* mSelectedFirstJetIndex;
    TH1F* mSelectedSecondJetIndex;

    TH1F* mSecondJetPhotonDeltaPhi;
    TH1F* mSecondJetPhotonDeltaR;
    TH1F* mSecondJetPhotonDeltaPt;

    TH1F* mSelectedFirstJetPhotonDeltaPhi;
    TH1F* mSelectedFirstJetPhotonDeltaR;

    TH1F* mSelectedSecondJetPhotonDeltaPhi;
    TH1F* mSelectedSecondJetPhotonDeltaR;

    // For B / C jets neutrinos
    TClonesArray* mNeutrinos;
    TClonesArray* mNeutrinosPDG;

    // Cache for MC particles
    bool mDumpAllMCParticles;
    std::unordered_map<const reco::Candidate*, int> mParticlesIndexes;

    void particleToTree(const reco::Candidate* particle, TTree* t, std::vector<boost::shared_ptr<void> >& addresses);
    
    void updateBranch(TTree* tree, void* address, const std::string& name, const std::string& type = "F");
    template<typename U>
      void updateBranch(TTree* tree, std::vector<U>*& address, const std::string& name);

    void updateBranchArray(TTree* tree, void* address, const std::string& name, const std::string& size, const std::string& type = "F");

    void photonToTree(const pat::PhotonRef& photonRef, pat::Photon& photon, const edm::Event& event);
    void metsToTree(const pat::MET& met, const pat::MET& rawMet, const std::vector<TTree*>& trees);
    void metToTree(const pat::MET* met, TTree* tree, TTree* genTree);
    void jetsToTree(const pat::Jet* firstJet, const pat::Jet* secondJet, const std::vector<TTree*>& trees);
    void jetToTree(const pat::Jet* jet, bool findNeutrinos, TTree* tree, TTree* genTree);
    void electronsToTree(const edm::Handle<pat::ElectronCollection>& electrons, const reco::Vertex& pv);
    void muonsToTree(const edm::Handle<pat::MuonCollection>& muons, const reco::Vertex& pv);

    int getMotherIndex(const edm::Handle<reco::GenParticleCollection>& genParticles, const reco::Candidate* mother);
    void genParticlesToTree(const edm::Handle<reco::GenParticleCollection>& genParticles);
  //FactorizedJetCorrector
  FactorizedJetCorrector *jetCorrector;
  FactorizedJetCorrector *jetCorrectorForTypeI;
  FactorizedJetCorrector *jetCorrectorForTypeIL1;
  std::vector<JetCorrectorParameters> vPar;
  std::vector<JetCorrectorParameters> vParTypeI;
  std::vector<JetCorrectorParameters> vParTypeIL1;
//giulia ---- regression no more necessary in 73X
  //define (once for all) corrector for regression
  //EnergyScaleCorrection_class *RegressionCorrector;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GammaJetFilter::GammaJetFilter(const edm::ParameterSet& iConfig):
  mIsMC(false), mIsValidLumiBlock(false)
{

  mIsMC = iConfig.getUntrackedParameter<bool>("isMC", "false");

  //sum_GenWeights = 0;

  //giulia --- no more necessary, regression integrated in release 73X
  //Photon energy regression corrector (need to define it once for data once for mc)
//if (mIsMC) {
//  RegressionCorrector = new EnergyScaleCorrection_class("",edm::FileInPath("JetMETCorrections/GammaJetFilter/data/step8-stochasticSmearing-invMass_SC_regrCorrSemiParV5_pho-loose-Et_20-trigger-noPF-HggRunEtaR9Et.dat").fullPath());
//  } else {
//   RegressionCorrector = new EnergyScaleCorrection_class(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/22Jan2012-runDepMCAll_v3-noR9shift-step8-invMass_SC_regrCorrSemiParV5_pho-loose-Et_20-trigger-noPF-HggRunEtaR9Et.dat").fullPath(),"");
//  }
//

  if (! mIsMC) {
    mJSONFile = iConfig.getParameter<std::string>("json");
    mCSVFile = iConfig.getParameter<std::string>("csv");
    mFilterData = iConfig.getUntrackedParameter<bool>("filterData", true);
    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    //CHS
   
    JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt").fullPath());    
    JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L2Relative_AK4PFchs.txt").fullPath());
    JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt").fullPath());
    //giulia--- corrector for type I is the same now 
    JetCorrectorParameters *L1JetParForTypeI = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt").fullPath());


    //JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V6_DATA_L2L3Residual_AK5PFchs.txt").fullPath());
    //JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V5_MC_L3Absolute_AK5PFchs.txt").fullPath());
    //JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V5_MC_L2Relative_AK5PFchs.txt").fullPath());
    //  //  JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V5_DATA_L1FastJet_AK5PFchs.txt").fullPath());
    //JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V6_DATA_L1FastJet_AK5PFchs.txt").fullPath());
    //Wi//nter14_V1_DATA_L1FastJet_AK5PFchs.txt
    ////txt file to use for L1 only for typeI
    //JetCorrectorParameters *L1JetParForTypeI = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V0_DATA_L1FastJetPU_AK5PFchs_pt.txt").fullPath()); 
    /*
    //NO CHS
    JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V5_DATA_L2L3Residual_AK5PF.txt").fullPath());
    JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V5_MC_L3Absolute_AK5PF.txt").fullPath());
    JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V5_MC_L2Relative_AK5PF.txt").fullPath());
    JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V5_DATA_L1FastJet_AK5PF.txt").fullPath());
    //txt file to use for L1 only for typeI
    JetCorrectorParameters *L1JetParForTypeI = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/Winter14_V0_DATA_L1FastJetPU_AK5PF_pt.txt").fullPath());
    */
    //
    //
    // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);
    //vPar.push_back(*ResJetPar); //comment if you dont want residuals
    jetCorrector = new FactorizedJetCorrector(vPar);
    //FAKE vPar for typeI fix
    vParTypeI.push_back(*L1JetPar);
    vParTypeI.push_back(*L2JetPar);
    vParTypeI.push_back(*L3JetPar);
    //vParTypeI.push_back(*ResJetPar); //comment if you dont want residuals
    jetCorrectorForTypeI = new FactorizedJetCorrector(vParTypeI);
    //FAKE vPar for typeI fix only L1
    vParTypeIL1.push_back(*L1JetParForTypeI);
    jetCorrectorForTypeIL1 = new FactorizedJetCorrector(vParTypeIL1);
    //
    //delete ResJetPar;
    delete L3JetPar;
    delete L2JetPar;
    delete L1JetPar;
    delete L1JetParForTypeI;
  } else {
    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved

    //CHS
    JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt").fullPath());    
    JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L2Relative_AK4PFchs.txt").fullPath());
    JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt").fullPath());
    //giulia--- corrector for type I is the same now 
    JetCorrectorParameters *L1JetParForTypeI = new JetCorrectorParameters(edm::FileInPath("JetMETCorrections/GammaJetFilter/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt").fullPath());


    //
    // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
    //When i thought it was L1(V0)L2L3 - L1(V0)
    //    vParTypeI.push_back(*L1JetParForTypeI);
    vParTypeI.push_back(*L1JetPar);
    vParTypeI.push_back(*L2JetPar);
    vParTypeI.push_back(*L3JetPar);
    jetCorrectorForTypeI = new FactorizedJetCorrector(vParTypeI);
    //FAKE vPar for typeI fix only L1
    vParTypeIL1.push_back(*L1JetParForTypeI);
    jetCorrectorForTypeIL1 = new FactorizedJetCorrector(vParTypeIL1);
    //
    delete L3JetPar;
    delete L2JetPar;
    delete L1JetParForTypeI;
  }

  mPhotonsIT = iConfig.getUntrackedParameter<edm::InputTag>("photons", edm::InputTag("slimmedPhotons"));
  mCorrPhotonWRegression = iConfig.getUntrackedParameter<bool>("doPhotonRegression", false);
  mJetsAK4PFlowIT = iConfig.getUntrackedParameter<edm::InputTag>("jetsAK4PFlow", edm::InputTag("slimmedJets"));
  mJetsAK8PFlowIT = iConfig.getUntrackedParameter<edm::InputTag>("jetsAK8PFlow", edm::InputTag("slimmedJetsAK8"));
  //mJetsAK4CaloIT = iConfig.getUntrackedParameter<edm::InputTag>("jetsAK4Calo", edm::InputTag("selectedPatJets"));
  //mJetsAK8CaloIT = iConfig.getUntrackedParameter<edm::InputTag>("jetsAK8Calo", edm::InputTag("selectedPatJetsCaloAK8"));
  mDoJEC         = iConfig.getUntrackedParameter<bool>("doJetCorrection", false);
  mRedoTypeI     = iConfig.getUntrackedParameter<bool>("redoTypeIMETCorrection", false);
  mDoFootprint     = iConfig.getUntrackedParameter<bool>("doFootprintMETCorrection", false);

  if (mDoJEC) {
    mJECFromRaw = iConfig.getUntrackedParameter<bool>("correctJecFromRaw", false);
    mCorrectorLabel = iConfig.getUntrackedParameter<std::string>("correctorLabel", "ak5PFResidual");
  }

  mFirstJetPtCut = iConfig.getUntrackedParameter<bool>("firstJetPtCut", true);
  mFirstJetThreshold = iConfig.getUntrackedParameter<double>("firstJetThreshold", 0.3);

  bool runOnCHS    = iConfig.getUntrackedParameter<bool>("runOnCHS", true);
  bool runOnNonCHS = iConfig.getUntrackedParameter<bool>("runOnNonCHS", true);

  bool runOnPFAK4    = iConfig.getUntrackedParameter<bool>("runOnPFAK4", true);
  bool runOnPFAK8    = iConfig.getUntrackedParameter<bool>("runOnPFAK8", false);
  bool runOnCaloAK4  = iConfig.getUntrackedParameter<bool>("runOnCaloAK4", false);
  bool runOnCaloAK8  = iConfig.getUntrackedParameter<bool>("runOnCaloAK8", false);
  //giulia --- default is chs, so the inputtag is the same
  //run module to do chs ab initio!
  edm::InputTag jetsAK4PFlowITchs = iConfig.getUntrackedParameter<edm::InputTag>("jetsAK4PFlowchs", edm::InputTag("slimmedJets"));
  edm::InputTag jetsAK8PFlowITchs = iConfig.getUntrackedParameter<edm::InputTag>("jetsAK8PFlowchs", edm::InputTag("slimmedJetsAK8"));

  if (! mIsMC && mFilterData) {
    // Load JSON file of good runs
    readJSONFile();
    readCSVFile();
  }

  edm::Service<TFileService> fs;
  mPhotonTree = fs->make<TTree>("photon", "photon tree");

  if (mIsMC){
    mPhotonGenTree = fs->make<TTree>("photon_gen", "photon gen tree");
  }  else{
    mPhotonGenTree = nullptr;
  }
  mAnalysisTree = fs->make<TTree>("analysis", "analysis tree");
  mMuonsTree = fs->make<TTree>("muons", "muons tree");
  mElectronsTree = fs->make<TTree>("electrons", "electrons tree");

  mTotalLuminosity = fs->make<TParameter<double> >("total_luminosity", 0.);

  mEventsWeight = 1.;
  mPtHatMin     = -1.;
  mPtHatMax     = -1.;

  if (mIsMC) {
    // Read cross section and number of generated events
    double crossSection = iConfig.getParameter<double>("crossSection");
    unsigned long long generatedEvents = iConfig.getParameter<unsigned long long>("generatedEvents");
    mEventsWeight = crossSection / (float) generatedEvents;

    mPtHatMin = iConfig.getUntrackedParameter<double>("ptHatMin", -1.);
    mPtHatMax = iConfig.getUntrackedParameter<double>("ptHatMax", -1.);
  }

  if (runOnNonCHS) {
    if (runOnPFAK4) {
      mJetCollections.push_back("PFlowAK4");
      mJetCollectionsData["PFlowAK4"] = {AK4, mJetsAK4PFlowIT};
    }
    if (runOnPFAK8) {
      mJetCollections.push_back("PFlowAK8");
      mJetCollectionsData["PFlowAK8"] = {AK8, mJetsAK8PFlowIT};
    }
  }

  if (runOnCHS) {
    if (runOnPFAK4) {
      mJetCollections.push_back("PFlowAK4chs");
      mJetCollectionsData["PFlowAK4chs"] = {AK4, jetsAK4PFlowITchs};
    }
    if (runOnPFAK8) {
      mJetCollections.push_back("PFlowAK8chs");
      mJetCollectionsData["PFlowAK8chs"] = {AK8, jetsAK8PFlowITchs};
    }
  }

  if (runOnCaloAK4) {
    mJetCollections.push_back("CaloAK4");
    mJetCollectionsData["CaloAK4"]  = {AK4, mJetsAK4CaloIT};
  }

  if (runOnCaloAK8) {
    mJetCollections.push_back("CaloAK8");
    mJetCollectionsData["CaloAK8"]  = {AK8, mJetsAK8CaloIT};
  }

  FOREACH(mJetCollections) {
    createTrees(*it, *fs);
  }

  mProcessedEvents = fs->make<TParameter<long long> >("total_events", 0);
  mSelectedEvents = fs->make<TParameter<long long> >("passed_events", 0);

  mJECRedone = fs->make<TParameter<bool> >("jec_redone", mDoJEC, '*');
  mFirstJetPtCutParameter = fs->make<TParameter<bool> >("cut_on_first_jet_pt", mFirstJetPtCut, '*');
  if (mDoJEC) {
    mJECFromRawParameter = fs->make<TParameter<bool> >("jec_from_raw_jet", mJECFromRaw, '*');
    mJECCorrectorLabel = fs->make<TNamed>("jec_corrector_label", mCorrectorLabel);
  }

  if (mFirstJetPtCut) {
    mFirstJetThresholdParameter = fs->make<TParameter<double> >("cut_on_first_jet_treshold", mFirstJetThreshold);
  }

  if (mIsMC){
    // to store the sum of weights
    h_sumW = fs->make<TH1F>("h_sumW", "h_sumW", 1, -0.5, 5.5);
    h_sumW->Sumw2();
  }  else{
    h_sumW = nullptr;
  }
  mFirstJetPhotonDeltaPhi = fs->make<TH1F>("firstJetPhotonDeltaPhi", "firstJetPhotonDeltaPhi", 50, 0., M_PI);
  mFirstJetPhotonDeltaR = fs->make<TH1F>("firstJetPhotonDeltaR", "firstJetPhotonDeltaR", 80, 0, 10);
  mFirstJetPhotonDeltaPt = fs->make<TH1F>("firstJetPhotonDeltaPt", "firstJetPhotonDeltaPt", 100, 0, 50);
  mFirstJetPhotonDeltaPhiDeltaR = fs->make<TH2F>("firstJetPhotonDeltaPhiDeltaR", "firstJetPhotonDeltaPhiDeltaR", 50, 0, M_PI, 80, 0, 10);

  mSelectedFirstJetIndex = fs->make<TH1F>("selectedFirstJetIndex", "selectedFirstJetIndex", 20, 0, 20);
  mSelectedSecondJetIndex = fs->make<TH1F>("selectedSecondJetIndex", "selectedSecondJetIndex", 20, 0, 20);

  mSecondJetPhotonDeltaPhi = fs->make<TH1F>("secondJetPhotonDeltaPhi", "secondJetPhotonDeltaPhi", 50, 0., M_PI);
  mSecondJetPhotonDeltaR = fs->make<TH1F>("secondJetPhotonDeltaR", "secondJetPhotonDeltaR", 80, 0, 10);
  mSecondJetPhotonDeltaPt = fs->make<TH1F>("secondJetPhotonDeltaPt", "secondJetPhotonDeltaPt", 100, 0, 50);

  mSelectedFirstJetPhotonDeltaPhi = fs->make<TH1F>("selectedFirstJetPhotonDeltaPhi", "selectedFirstJetPhotonDeltaPhi", 50, 0., M_PI);
  mSelectedFirstJetPhotonDeltaR = fs->make<TH1F>("selectedFirstJetPhotonDeltaR", "selectedFirstJetPhotonDeltaR", 80, 0, 10);

  mSelectedSecondJetPhotonDeltaPhi = fs->make<TH1F>("selectedSecondJetPhotonDeltaPhi", "selectedSecondJetPhotonDeltaPhi", 50, 0., M_PI);
  mSelectedSecondJetPhotonDeltaR = fs->make<TH1F>("selectedSecondJetPhotonDeltaR", "selectedSecondJetPhotonDeltaR", 80, 0, 10);

  mPFIsolator.initializePhotonIsolation(true);
  mPFIsolator.setConeSize(0.3);

  mNeutrinos = NULL;
  mNeutrinosPDG = NULL;
  if (mIsMC) {
    mNeutrinos = new TClonesArray("TLorentzVector", 3);
    mNeutrinosPDG = new TClonesArray("TParameter<int>", 3);
  }
}


GammaJetFilter::~GammaJetFilter()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  delete mNeutrinos;
  delete mNeutrinosPDG;
}

void GammaJetFilter::createTrees(const std::string& rootName, TFileService& fs) {

  TFileDirectory dir = fs.mkdir(rootName);
  std::vector<TTree*>& trees = mJetTrees[rootName];

  trees.push_back(dir.make<TTree>("first_jet", "first jet tree"));
  trees.push_back(dir.make<TTree>("second_jet", "second jet tree"));

  trees.push_back(dir.make<TTree>("first_jet_raw", "first raw jet tree"));
  trees.push_back(dir.make<TTree>("second_jet_raw", "second raw jet tree"));

  if (mIsMC) {
    trees.push_back(dir.make<TTree>("first_jet_gen", "first gen jet tree"));
    trees.push_back(dir.make<TTree>("second_jet_gen", "second gen jet tree"));
  } else {
    trees.push_back(nullptr);
    trees.push_back(nullptr);
  }

  // MET
  std::vector<TTree*>& met = mMETTrees[rootName];
  met.push_back(dir.make<TTree>("met", "met tree"));
  met.push_back(dir.make<TTree>("met_raw", "met raw tree"));

  if (mIsMC)
    met.push_back(dir.make<TTree>("met_gen", "met gen tree"));
  else
    met.push_back(nullptr);

  // Misc
  mMiscTrees[rootName] = dir.make<TTree>("misc", "misc tree");
}

void GammaJetFilter::updateBranch(TTree* tree, void* address, const std::string& name, const std::string& type/* = "F"*/) {
  TBranch* branch = tree->GetBranch(name.c_str());
  if (branch == NULL) {
    branch = tree->Branch(name.c_str(), address, std::string(name + "/" + type).c_str()); 
  } else {
    branch->SetAddress(address);
  }
}

template<typename U> void GammaJetFilter::updateBranch(TTree* tree, std::vector<U>*& address, const std::string& name) {
  TBranch* branch = tree->GetBranch(name.c_str());
  if (branch == NULL) {
    branch = tree->Branch(name.c_str(), &address); 
  } else {
    branch->SetAddress(&address);
  }
}

void GammaJetFilter::updateBranchArray(TTree* tree, void* address, const std::string& name, const std::string& size, const std::string& type/* = "F"*/) {
  TBranch* branch = tree->GetBranch(name.c_str());
  if (branch == NULL) {
    branch = tree->Branch(name.c_str(), address, std::string(name + "[" + size + "]/" +type).c_str()); 
  } else {
    branch->SetAddress(address);
  }
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool GammaJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  mProcessedEvents->SetVal(mProcessedEvents->GetVal() + 1);


  if (! mIsMC && mFilterData && ! mIsValidLumiBlock) {
    return false;
  }

  // Vertex
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices", vertices);


  // Keep events with at least one vertex
  if (!vertices.isValid() || vertices->size() == 0 || vertices->front().isFake())
    return false;

  const reco::Vertex& primaryVertex = vertices->at(0);

  double generatorWeight = 1.;

  if (mIsMC) {
    edm::Handle<GenEventInfoProduct> eventInfos;
    iEvent.getByLabel("generator", eventInfos);
    if (eventInfos.isValid() && eventInfos->hasBinningValues()) {
      double genPt = eventInfos->binningValues()[0];

      if (mPtHatMin >= 0. && genPt < mPtHatMin)
      	return false;

      if (mPtHatMax >= 0. && genPt > mPtHatMax)
	return false;
    }

    generatorWeight = eventInfos->weight();
    h_sumW->Fill(0.,generatorWeight);
    if (generatorWeight == 0.) {
      generatorWeight = 1.;
    }
    //sum_GenWeights += generatorWeight;
  }

  edm::Handle<double> pFlowRho;
  iEvent.getByLabel(edm::InputTag("offlineSlimmedPrimaryVertices"), pFlowRho); // For photon ID

  // Necesseray collection for calculate sigmaIPhiIPhi
  // 2011 Photon ID
  edm::ESHandle<CaloTopology> topology;
  iSetup.get<CaloTopologyRecord>().get(topology);

  edm::Handle<pat::PhotonCollection> photons; //handle for pat::photonref
  iEvent.getByLabel(mPhotonsIT, photons);
  pat::PhotonCollection photons_nonconst = *photons; //handle for pat::photon
  pat::PhotonCollection::iterator it = photons_nonconst.begin();
  std::vector<pat::Photon> photonsVec;
  pat::Photon pho_tmp;

  uint32_t index = 0;
  uint32_t goodPhoIndex = -1;
  for (; it != photons_nonconst.end(); ++it, index++) { 
    pho_tmp=*it;
    if (fabs(it->eta()) <= 1.3) {
      pat::PhotonRef PhotonReftmp(photons, index);
      if (isValidPhotonEB_PHYS14(PhotonReftmp, iEvent)) {
      	photonsVec.push_back(*it);
	goodPhoIndex=index;
      }
    }
  }

  // Only one good photon per event
  if (photonsVec.size() != 1)    return false;
  pat::Photon photon = photonsVec[0];
  pat::PhotonRef GoodphotonRef(photons, goodPhoIndex);

  //for technical reasons i need a photonref and a photon. 
  //Since there is only one photon in these events, we are sur that the 
  //goodPhoIndex is referring to the same photon as photonsVec[0];
  //
  //float regressionCorr=1.;


  //giulia --- regression integrated in 73X , putting FALSE by hand!
  mCorrPhotonWRegression = false;
  if (mCorrPhotonWRegression) {
    
    //calculate the regression energy using photonRef and getting the reco object
    edm::Handle<edm::ValueMap<float>> regressionEnergyHandle;
    iEvent.getByLabel(edm::InputTag("eleNewEnergiesProducer", "energySCEleJoshPhoSemiParamV5ecorr", "PAT"),regressionEnergyHandle);
    edm::Ptr<reco::Candidate> GoodrecoObject = GoodphotonRef->originalObjectRef();
    //float GoodregressionEnergy = (*regressionEnergyHandle)[GoodrecoObject] ;
    //correct this january 15
    //  regressionCorr = GoodregressionEnergy/(GoodphotonRef->energy());
    //rescale the photon to the regression energy (rescale the whole p4 by the ratio of regression energy over uncorrected energy)
    //correct this january 15
    //  photon.setP4(photon.p4()*regressionCorr);
    //now apply the additional correction (data) or smearing (mc)
    //
    int processingdata=1;
    if(mIsMC) processingdata=0;
      correctPhoton(photon,iEvent, processingdata, int(vertices->size()));
  }



    // Process jets
    edm::Handle<pat::JetCollection> jetsHandle;

    FOREACH(mJetCollections) {

      JetInfos infos = mJetCollectionsData[*it];

      iEvent.getByLabel(infos.inputTag, jetsHandle);
      pat::JetCollection jets = *jetsHandle;
      if (mDoJEC) {
	correctJets(jets, iEvent, iSetup);
      } else {
	extractRawJets(jets);
      }
      //giulia --- comment QG taggiung stuff
      //edm::Handle<edm::ValueMap<float>>  qgTagHandleMLP;
      //edm::Handle<edm::ValueMap<float>>  qgTagHandleLikelihood;
      //iEvent.getByLabel("QGTagger" + *it,"qgMLP", qgTagHandleMLP);
      //iEvent.getByLabel("QGTagger" + *it,"qgLikelihood", qgTagHandleLikelihood);
      processJets(&photon, jets, infos.algo, /*qgTagHandleMLP, qgTagHandleLikelihood, */jetsHandle, mJetTrees[*it]);

      // MET
      // giulia --- taking MT from miniAOD --- assuming is it corrected -- ??
      // which algo is used to compute MET? -- AK4?
      edm::Handle<pat::METCollection> metsHandle;
      //iEvent.getByLabel(std::string("patMETs" + ((*it == "AK4Calo") ? "" : *it)), metsHandle);
      iEvent.getByLabel(edm::InputTag("slimmedMETs"), metsHandle);

      //giulia - raw MET TO BE CHANGED		
      edm::Handle<pat::METCollection> rawMets;
      //iEvent.getByLabel(std::string("patPFMet" + ((*it == "AK4Calo") ? "" : *it)), rawMets);
      iEvent.getByLabel(edm::InputTag("slimmedMETs"),  rawMets);

      pat::METCollection mets = *metsHandle;
      pat::MET& met = mets[0];
      const pat::MET& rawMet = rawMets->at(0);

      if (mDoJEC || mRedoTypeI) {
      	//cout << "WARNING: keeping MET from miniAOD... not clear wether/how it is corrected..." << endl;

        if (mDoFootprint) {
	  //giulia --- comment footprint stuff
	  //correctMETWithFootprintAndTypeI(rawMet, met, jets, iEvent, photon, GoodphotonRef);
	  cout << "WARNING: Removed footprint stuff for debugging..." << endl;
	} else {
	  if (mCorrPhotonWRegression) {
	    correctMETWithRegressionAndTypeI(rawMet, met, jets, iEvent, photon, GoodphotonRef);
	  } else {
	    correctMETWithTypeI(rawMet, met, jets, iEvent);
	  }
	}
      }

      if (rawMets.isValid())
	metsToTree(met, rawMet, mMETTrees[*it]);
      else {
	pat::MET emptyRawMet = pat::MET();
	metsToTree(met, emptyRawMet, mMETTrees[*it]);
      }
      //

      // Rho
      edm::Handle<double> rhos;
      if (it->find("Calo") != std::string::npos)
	//iEvent.getByLabel(edm::InputTag("kt6CaloJets", "rho"), rhos);
	iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetAllCalo"), rhos);

      else
	//iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho"), rhos);
	iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetAll"), rhos);

      double rho = *rhos;
      updateBranch(mMiscTrees[*it], &rho, "rho", "D");

      mMiscTrees[*it]->Fill();
    }

    // Number of vertices for pu reweighting
    edm::Handle<std::vector<PileupSummaryInfo> > puInfos;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), puInfos);

    float nTrueInteractions = -1;
    int nPUVertex = -1;
    unsigned int nVertex = vertices->size();

    edm::EventID eventId = iEvent.id();
    EventNumber_t event = eventId.event();
    RunNumber_t run = eventId.run();
    LuminosityBlockNumber_t lumiBlock = eventId.luminosityBlock();

    if (mIsMC) {
      for (std::vector<PileupSummaryInfo>::const_iterator it = puInfos->begin(); it != puInfos->end();
	   ++it) {

  	int BX = it->getBunchCrossing();
	if (BX == 0) {
	  nPUVertex = it->getPU_NumInteractions();
	  nTrueInteractions = it->getTrueNumInteractions();
	  break;
	}
      }

      if (nPUVertex < 0) {
	throw cms::Exception("PUReweighting") << "No in-time beam crossing found!" << std::endl;
      }
    } else {
      nTrueInteractions = mCurrentTruePU;
    }

    updateBranch(mAnalysisTree, &run, "run", "i");
    updateBranch(mAnalysisTree, &lumiBlock, "lumi_block", "i");
    updateBranch(mAnalysisTree, &event, "event", "i");
    updateBranch(mAnalysisTree, &nVertex, "nvertex", "i");
    updateBranch(mAnalysisTree, &nTrueInteractions, "ntrue_interactions");
    updateBranch(mAnalysisTree, &nPUVertex, "pu_nvertex", "I");
    updateBranch(mAnalysisTree, &mEventsWeight, "event_weight"); // Only valid for binned samples
    updateBranch(mAnalysisTree, &generatorWeight, "generator_weight", "D"); // Only valid for flat samples
    //updateBranch(mAnalysisTree, &sum_GenWeights, "sum_GenWeights", "D"); // Only valid for flat samples

    // Triggers
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), triggerResults);

    std::vector<std::string>* trigNames = new std::vector<std::string>();
    std::vector<bool>* trigResults = new std::vector<bool>();

    if (triggerResults.isValid()) {
      static std::vector<boost::regex> validTriggers = { boost::regex("HLT_.*Photon.*", boost::regex_constants::icase) };
      const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);

      size_t size = triggerResults->size();

      for (size_t i = 0; i < size; i++) {
	std::string triggerName = triggerNames.triggerName(i);
	bool isValid = false;
	for (boost::regex& validTrigger: validTriggers) {
	  if (boost::regex_match(triggerName, validTrigger)) {
	    isValid = true;
	    break;
	  }
	}

  	if (!isValid)
	  continue;

  	unsigned int index = triggerNames.triggerIndex(triggerName);
	bool passed = triggerResults->accept(index);

  	trigResults->push_back(passed);
	trigNames->push_back(triggerName);
      }
    }

    // Create branches, even if they're empty
    updateBranch(mAnalysisTree, trigNames, "trigger_names");
    updateBranch(mAnalysisTree, trigResults, "trigger_results");

    mAnalysisTree->Fill();

    delete trigNames;
    delete trigResults;

    photonToTree(GoodphotonRef, photon, iEvent);

    // Electrons
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByLabel("slimmedElectrons", electrons);
    electronsToTree(electrons, primaryVertex);

    // Muons
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByLabel("slimmedMuons", muons);
    muonsToTree(muons, primaryVertex);

    mSelectedEvents->SetVal(mSelectedEvents->GetVal() + 1);
    return true;
}



// giulia -- correct energy of the photon already implemented in 73X
// smearing - scale factors not necessary until we have new data @13 TeV
//DO NOTHING!
void GammaJetFilter::correctPhoton(pat::Photon& photon, edm::Event& iEvent, int isData, int nPV) {
//  edm::EventID eventId = iEvent.id();

//  if(isData==1) {
//    
//    float scalecorr = RegressionCorrector->ScaleCorrection(eventId.run(),true,photon.r9(),photon.eta(),photon.pt(),nPV,15.); //nPVmean is not actually used in the function, so it is set to a dummy value
//    
//    
//    //ScaleCorrection(int runNumber, bool isEBEle, double R9Ele, double etaSCEle, double EtEle, int nPV, float nPVmean=0);
//    photon.setP4(photon.p4()*scalecorr);
//  } else {
//    reco::SuperClusterRef superCluster = photon.superCluster();
//    float smearcorr = RegressionCorrector->getSmearing(eventId.run(),photon.energy(),true,photon.r9(),superCluster->eta());
//    //getSmearing(int runNumber, float energy, bool isEBEle, float R9Ele, float etaSCEle);
//    photon.setP4(photon.p4()*smearcorr);
//  }
}


void GammaJetFilter::correctJets(pat::JetCollection& jets, edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get Jet corrector
  const JetCorrector* corrector = JetCorrector::getJetCorrector(mCorrectorLabel, iSetup);

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    // Store raw jet, it's not possible to get it after corrections
    pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true); // Store raw jet inside our jet. This allow us to correctly sort the resulting collection
    pat::Jet L1Jet  = jet.correctedJet("L1FastJet");
    jet.addUserData("L1Jet", L1Jet, true); // Embed L1 corrected jet for TypeI correction

    if (mJECFromRaw) {
      double toRaw = jet.jecFactor("Uncorrected");
      jet.setP4(jet.p4() * toRaw); // It's now a raw jet
    }

    double corrections =1.;
    if(mIsMC)    {
      corrections = corrector->correction(jet, iEvent, iSetup);
    } else {
      edm::Handle<double> rho_;
      iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho"), rho_);
      jetCorrector->setJetEta(jet.eta());
      jetCorrector->setJetPt(jet.pt());
      jetCorrector->setJetA(jet.jetArea());
      jetCorrector->setRho(*rho_);
      corrections = jetCorrector->getCorrection();
    }
    jet.scaleEnergy(corrections);
  }

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}


void GammaJetFilter::correctMETWithTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, edm::Event& event) {
  double deltaPx = 0., deltaPy = 0.;
  // See https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=174324 slide 4
  // and http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h?revision=1.8&view=markup
  for (pat::JetCollection::const_iterator it = jets.begin(); it != jets.end(); ++it) {
    const pat::Jet& jet = *it;

    if (jet.pt() > 10) {

      const pat::Jet* rawJet = jet.userData<pat::Jet>("rawJet");
/*//without typei fix
      const pat::Jet* L1Jet  = jet.userData<pat::Jet>("L1Jet");
      reco::Candidate::LorentzVector L1JetP4  = L1Jet->p4();
*/

//with typeI fix
    double corrsForTypeI =1.;
    double corrsForTypeIL1=1.;
    edm::Handle<double> rho_;
    event.getByLabel(edm::InputTag("fixedGridRhoFastjetAllCalo"), rho_);

    jetCorrectorForTypeIL1->setJetEta(rawJet->eta());
    jetCorrectorForTypeIL1->setJetPt(rawJet->pt());
    jetCorrectorForTypeIL1->setJetA(rawJet->jetArea());
    jetCorrectorForTypeIL1->setRho(*rho_);
    corrsForTypeIL1 = jetCorrectorForTypeIL1->getCorrection();

    pat::Jet jetL1 = *rawJet;
    jetL1.scaleEnergy(corrsForTypeIL1);

    jetCorrectorForTypeI->setJetEta(rawJet->eta());
    jetCorrectorForTypeI->setJetPt(rawJet->pt());
    jetCorrectorForTypeI->setJetA(rawJet->jetArea());
    jetCorrectorForTypeI->setRho(*rho_);
    corrsForTypeI = jetCorrectorForTypeI->getCorrection();

    pat::Jet jet = *rawJet;
    jet.scaleEnergy(corrsForTypeI);

//typeI fix
      reco::Candidate::LorentzVector L1JetP4  = jetL1.p4();

      double emEnergyFraction = rawJet->chargedEmEnergyFraction() + rawJet->neutralEmEnergyFraction();
      if (emEnergyFraction > 0.90)
        continue;

////without typei fix
//     //reco::Candidate::LorentzVector rawJetP4 = rawJet->p4();
//      reco::Candidate::LorentzVector L1JetP4  = L1Jet->p4();

      deltaPx += (jet.px() - L1JetP4.px());
      deltaPy += (jet.py() - L1JetP4.py());
    }
  }

  double correctedMetPx = rawMet.px() - deltaPx;
  double correctedMetPy = rawMet.py() - deltaPy;
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);

  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));
}

//giulia --- comment footprint stuff

//void GammaJetFilter::correctMETWithFootprintAndTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets,  edm::Event& event, pat::Photon& photon, const pat::PhotonRef& photonRef) {
////retrieve the footprint corrections to MET vector
//edm::Handle<edm::ValueMap<double>> footpxHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintPx", "PAT"), footpxHandle);
//edm::Handle<edm::ValueMap<double>> footpyHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintPy", "PAT"), footpyHandle);
//
//edm::Handle<double> footprintMExrawHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintMExraw", "PAT"), footprintMExrawHandle);
//edm::Handle<double> footprintMEyrawHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintMEyraw", "PAT"), footprintMEyrawHandle);
//
//edm::Handle<edm::ValueMap<double>> footprintMExCorrHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintMExCorr", "PAT"), footprintMExCorrHandle);
//edm::Handle<edm::ValueMap<double>> footprintMEyCorrHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintMEyCorr", "PAT"), footprintMEyCorrHandle);
//
//double footprintMExCorr = (*footprintMExCorrHandle)[photonRef];
//double footprintMEyCorr = (*footprintMEyCorrHandle)[photonRef];
//
//  double deltaPx = 0., deltaPy = 0.;
//
//  // See https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=174324 slide 4
//  // and http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h?revision=1.8&view=markup
// //TypeI fix : use different L1corrections - only for typeI calculation! - 
//
//  for (pat::JetCollection::const_iterator it = jets.begin(); it != jets.end(); ++it) { 
////embed raw and l1 jet
//     const pat::Jet* rawJet = it->userData<pat::Jet>("rawJet");
////apply the ad hoc corrections
////calculate the corrections
//    double corrsForTypeI =1.;
//    double corrsForTypeIL1=1.;
//    edm::Handle<double> rho_;
//    event.getByLabel(edm::InputTag("kt6PFJets", "rho"), rho_);
////
//    jetCorrectorForTypeIL1->setJetEta(rawJet->eta());
//    jetCorrectorForTypeIL1->setJetPt(rawJet->pt());
//    jetCorrectorForTypeIL1->setJetA(rawJet->jetArea());
//    jetCorrectorForTypeIL1->setRho(*rho_);
//    corrsForTypeIL1 = jetCorrectorForTypeIL1->getCorrection();
//
//    pat::Jet jetL1 = *rawJet;
//    jetL1.scaleEnergy(corrsForTypeIL1);
////
//    jetCorrectorForTypeI->setJetEta(rawJet->eta()); 
//    jetCorrectorForTypeI->setJetPt(rawJet->pt());
//    jetCorrectorForTypeI->setJetA(rawJet->jetArea());
//    jetCorrectorForTypeI->setRho(*rho_);
//    corrsForTypeI = jetCorrectorForTypeI->getCorrection();
//
//    pat::Jet jet = *rawJet;
//    jet.scaleEnergy(corrsForTypeI);
////go ahead with typeI
//    if (jet.pt() > 10) {
//
//      double emEnergyFraction = rawJet->chargedEmEnergyFraction() + rawJet->neutralEmEnergyFraction();
//      if (emEnergyFraction > 0.90)
//        continue;
//
//      reco::Candidate::LorentzVector L1JetP4  = jetL1.p4();
//
//      deltaPx += (jet.px() - L1JetP4.px());
//      deltaPy += (jet.py() - L1JetP4.py());
//    }
//  }
////used for footprint correction
//  double correctedMetPx = footprintMExCorr + photonRef->px() - photon.px() - deltaPx;
//  double correctedMetPy = footprintMEyCorr  + photonRef->py() - photon.py() - deltaPy;
//  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);
///*
//cout<< "old MET = "<<rawMet.pt()  << endl;
//cout <<"recalculated MET "<< sqrt(pow(footprintMExraw,2)+pow(footprintMEyraw,2))<<endl;
//cout <<"corrected MET "<< correctedMetPt << endl;
//cout <<"corrected MET without photon "<< sqrt(pow(rawMet.px() + footprintpx,2)+pow(rawMet.py() + footprintpy,2))<<endl;
//cout <<"corrected MET without photon and with typeI "<< sqrt(pow(rawMet.px() + footprintpx - deltaPx,2)+pow(rawMet.py() + footprintpy - deltaPy,2))<<endl;
//cout <<"corrected MET adding back the photon "<< sqrt(pow(rawMet.px() + footprintpx - photonRef->px();,2)+pow(rawMet.py() + footprintpy - photonRef->py();,2))<<endl;
//cout <<"footprint MET "<< sqrt(pow(footprintMExCorr - deltaPx,2)+pow(footprintMEyCorr - deltaPy,2))<<endl;
//cout<< "" << endl;
//*/
//  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));
//}
//


void GammaJetFilter::correctMETWithRegressionAndTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets,  edm::Event& event, pat::Photon& photon, const pat::PhotonRef& photonRef) {
  //photonRef is the one before regression
  //photon is the one after


 double deltaPx = 0., deltaPy = 0.;
  for (pat::JetCollection::const_iterator it = jets.begin(); it != jets.end(); ++it) {
    const pat::Jet& jet = *it;

    if (jet.pt() > 10) {

      const pat::Jet* rawJet = jet.userData<pat::Jet>("rawJet");
    double corrsForTypeI =1.;
    double corrsForTypeIL1=1.;
    edm::Handle<double> rho_;
    event.getByLabel(edm::InputTag("fixedGridRhoFastjetAllCalo"), rho_);

    jetCorrectorForTypeIL1->setJetEta(rawJet->eta());
    jetCorrectorForTypeIL1->setJetPt(rawJet->pt());
    jetCorrectorForTypeIL1->setJetA(rawJet->jetArea());
    jetCorrectorForTypeIL1->setRho(*rho_);
    corrsForTypeIL1 = jetCorrectorForTypeIL1->getCorrection();

    pat::Jet jetL1 = *rawJet;
    jetL1.scaleEnergy(corrsForTypeIL1);

    jetCorrectorForTypeI->setJetEta(rawJet->eta());
    jetCorrectorForTypeI->setJetPt(rawJet->pt());
    jetCorrectorForTypeI->setJetA(rawJet->jetArea());
    jetCorrectorForTypeI->setRho(*rho_);
    corrsForTypeI = jetCorrectorForTypeI->getCorrection();

    pat::Jet jet = *rawJet;
    jet.scaleEnergy(corrsForTypeI);

    reco::Candidate::LorentzVector L1JetP4  = jetL1.p4();

    double emEnergyFraction = rawJet->chargedEmEnergyFraction() + rawJet->neutralEmEnergyFraction();
    if (emEnergyFraction > 0.90)

      deltaPx += (jet.px() - L1JetP4.px());
    deltaPy += (jet.py() - L1JetP4.py());
    }
  }

  double correctedMetPx = rawMet.px() + photonRef->px() - photon.px() - deltaPx;
  double correctedMetPy = rawMet.py() + photonRef->py() - photon.py() - deltaPy;
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);

  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));

}


void GammaJetFilter::extractRawJets(pat::JetCollection& jets) {

  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it) {
    pat::Jet& jet = *it;

    const pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true);
    const pat::Jet L1Jet  = jet.correctedJet("L1FastJet");
    jet.addUserData("L1Jet", L1Jet, true); // Embed L1 corrected jet for TypeI correction
  }

}

void GammaJetFilter::processJets(pat::Photon* photon, pat::JetCollection& jets, const JetAlgorithm algo,/* edm::Handle<edm::ValueMap<float>>& qgTagMLP, edm::Handle<edm::ValueMap<float>>& qgTagLikelihood,*/ const edm::Handle<pat::JetCollection>& handleForRef, std::vector<TTree*>& trees) {

  pat::JetCollection selectedJets;

  pat::JetCollection::iterator it = jets.begin();
  uint32_t index = 0;
  uint32_t goodJetIndex = -1;
  for (; it != jets.end(); ++it, index++) {

    if (! isValidJet(*it))
      continue;

    goodJetIndex++;

    if (goodJetIndex == 0) {
      mFirstJetPhotonDeltaPhi->Fill(fabs(reco::deltaPhi(*photon, *it)));
      mFirstJetPhotonDeltaR->Fill(reco::deltaR(*photon, *it));
      mFirstJetPhotonDeltaPt->Fill(fabs(photon->pt() - it->pt()));

      mFirstJetPhotonDeltaPhiDeltaR->Fill(fabs(reco::deltaPhi(*photon, *it)), reco::deltaR(*photon, *it));
    } else if (goodJetIndex == 1) {
      mSecondJetPhotonDeltaPhi->Fill(fabs(reco::deltaPhi(*photon, *it)));
      mSecondJetPhotonDeltaR->Fill(reco::deltaR(*photon, *it));
      mSecondJetPhotonDeltaPt->Fill(fabs(photon->pt() - it->pt()));
    }

    //giulia --- comment QG tagging stuff
    // Extract Quark Gluon tagger value
    pat::JetRef jetRef(handleForRef, index);
    //it->addUserFloat("qgTagMLP", (*qgTagMLP)[jetRef]);
    //it->addUserFloat("qgTagLikelihood", (*qgTagLikelihood)[jetRef]);

    const double deltaR_threshold = (algo == AK4) ? 0.4 : 0.8;

    if (selectedJets.size() == 0) {
      // First jet selection

      if (index > 1) {
      	// It's the third jet of the event. We only want to consider the first two jets for our leadind jet,
	// so, throw this event
	break;
      }

      const double deltaPhi = reco::deltaPhi(*photon, *it);
      if (fabs(deltaPhi) < M_PI / 2.)
	continue; // Only back 2 back event are interesting

      const double deltaR = reco::deltaR(*photon, *it);
      if (deltaR < deltaR_threshold) // This jet is inside the photon. This is probably the photon mis-reconstructed as a jet
      	continue;

      // Jet are ordered by pt value.
      // Events are supposed to be balanced between Jet and Gamma
      // If the leading jet has less than 30% of the Photon pt,
      // dump the event as it's not interesting
      if (mFirstJetPtCut && (it->pt() < photon->pt() * mFirstJetThreshold))
      	break;

      mSelectedFirstJetIndex->Fill(goodJetIndex);
      selectedJets.push_back(*it);

    } else {

      // Second jet selection
      const double deltaR = reco::deltaR(*photon, *it);

      if (deltaR > deltaR_threshold) {
	mSelectedSecondJetIndex->Fill(goodJetIndex);
	selectedJets.push_back(*it);
      } else {
      	continue;
      }

      break;
    }

  }

  const pat::Jet* firstJet = NULL;
  const pat::Jet* secondJet = NULL;

  if (selectedJets.size() > 0) {

    firstJet = &selectedJets[0];
    mSelectedFirstJetPhotonDeltaPhi->Fill(fabs(reco::deltaPhi(*photon, *firstJet)));
    mSelectedFirstJetPhotonDeltaR->Fill(reco::deltaR(*photon, *firstJet));

    if (selectedJets.size() > 1) {
      secondJet = &selectedJets[1];

      mSelectedSecondJetPhotonDeltaPhi->Fill(fabs(reco::deltaPhi(*photon, *secondJet)));
      mSelectedSecondJetPhotonDeltaR->Fill(reco::deltaR(*photon, *secondJet));
    }
  }

  jetsToTree(firstJet, secondJet, trees);

  return;
}

// ------------ method called once each job just before starting event loop  ------------
void GammaJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GammaJetFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool GammaJetFilter::beginRun(edm::Run& run, edm::EventSetup const&)
{
  if (! mIsMC && mFilterData) {
    // Check if this run is valid
    std::stringstream stream;
    stream << run.run();

    if (! mValidRuns->isMember(stream.str()))
      return false; // Drop run

    mCurrentRunValidLumis.reset(new Json::Value((*mValidRuns)[stream.str()]));
  }

  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool GammaJetFilter::endRun(edm::Run& run, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool GammaJetFilter::beginLuminosityBlock(edm::LuminosityBlock& lumiBlock, edm::EventSetup const&)
{
  if (! mIsMC && mFilterData) {

    mIsValidLumiBlock = false;

    if (! mCurrentRunValidLumis.get())
      return false;

    // Check if this lumi block is valid
    assert(mCurrentRunValidLumis->isArray());
    for (Json::ArrayIndex i = 0; i < mCurrentRunValidLumis->size(); i++) {
      Json::Value lumiRange = (*mCurrentRunValidLumis)[i];

      assert(lumiRange.size() == 2);
      edm::LuminosityBlockNumber_t lumi = lumiBlock.luminosityBlock();
      if (lumi >= lumiRange[0].asUInt64() && lumi <= lumiRange[1].asUInt64()) {
      	mIsValidLumiBlock = true;

	mCurrentTruePU = mTruePUByLS[std::make_pair(lumiBlock.id().run(), lumiBlock.id().luminosityBlock())];
	return true;
      }
    }


    return false;
  }

  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool GammaJetFilter::endLuminosityBlock(edm::LuminosityBlock& lumiBlock, edm::EventSetup const&)
{
  if (! mIsMC && mFilterData && mIsValidLumiBlock) {
    updateLuminosity(lumiBlock);
  }
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GammaJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool GammaJetFilter::isValidJet(const pat::Jet& jet) {
  // First, check if this pat::Jet has a gen jet
  if (mIsMC && !jet.genJet()) {
    return false;
  }

  if (jet.isPFJet()) {

    // Jet ID
    // From https://twiki.cern.ch/twiki/bin/view/CMS/JetID

    // Jet ID works on uncorrected jets. *EnergyFraction take that into account when calculating the fraction,
    // so there's *NO* need to use an uncorrected jet
    bool isValid = jet.neutralHadronEnergyFraction() < 0.99;
    isValid &= jet.neutralEmEnergyFraction() < 0.99;
    //isValid &= jet.getPFConstituents().size() > 1;
    if (fabs(jet.eta()) < 2.4) {
      isValid &= jet.chargedHadronEnergyFraction() > 0.;
      isValid &= jet.chargedMultiplicity() > 0;
      isValid &= jet.chargedEmEnergyFraction() < 0.99;
    }

    return isValid;

  } else if (jet.isCaloJet() || jet.isJPTJet()) {

    if (! mCaloJetID.get()) {
      mCaloJetID.reset(new JetIDSelectionFunctor(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE));
      mCaloJetIDRet = mCaloJetID->getBitTemplate();
    }

    mCaloJetIDRet.set(false);
    return (*mCaloJetID)(jet, mCaloJetIDRet);

  } else {
    throw cms::Exception("UnsupportedJetType")
      << "Only PF and Calo jets are supported at this time" << std::endl;
  }

  return false;
}

enum class IsolationType {
  CHARGED_HADRONS,
  NEUTRAL_HADRONS,
  PHOTONS
};


//updated effective areas for PHYS14 -> https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2
float getEffectiveArea(float eta, IsolationType type) {
  eta = fabs(eta);
  switch (type) {
    case IsolationType::CHARGED_HADRONS:
      if (eta < 1.0)
        return 0.0130;
      else if (eta < 1.479)
      	return 0.0096;
      else if (eta < 2.0)
        return 0.0107;
      else if (eta < 2.2)
      	return 0.0077;
      else if (eta < 2.3)
      	return 0.0088;
      else if (eta < 2.4)
        return 0.0065;
      else
        return 0.0030;
      break;

    case IsolationType::NEUTRAL_HADRONS:
      if (eta < 1.0)
        return 0.0056;
      else if (eta < 1.479)
        return 0.0107;
      else if (eta < 2.0)
        return 0.0019;
      else if (eta < 2.2)
      	return 0.0011;
      else if (eta < 2.3)
      	return 0.0077;
      else if (eta < 2.4)
	return 0.0178;
      else
      	return 0.1675;
      break;

    case IsolationType::PHOTONS:
      if (eta < 1.0)
        return 0.0896;
      else if (eta < 1.479)
      	return 0.0762;
      else if (eta < 2.0)
      	return 0.0383;
      else if (eta < 2.2)
	return 0.0534;
      else if (eta < 2.3)
      	return 0.0846;
      else if (eta < 2.4)
        return 0.1032;
      else
        return 0.1598;
      break;
  }

  return -1;
}

double getCorrectedPFIsolation(double isolation, double rho, float eta, IsolationType type) {
  float effectiveArea = getEffectiveArea(eta, type);

  return std::max(isolation - rho * effectiveArea, 0.);
}

//-----------------
// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2 -- tight WP 76% 
bool GammaJetFilter::isValidPhotonEB_PHYS14(const pat::PhotonRef& photonRef, edm::Event& event) {
  if (mIsMC && !photonRef->genPhoton())
    return false;

  bool isValid = true;

  isValid &= photonRef->hadTowOverEm() < 0.011;
  isValid &= photonRef->sigmaIetaIeta() < 0.0099;

  if (! isValid)
    return false;

  edm::Handle<double> rhos;
  event.getByLabel(edm::InputTag("fixedGridRhoFastjetAll"), rhos);
  double rho = *rhos;


  //giulia --  use the method " passElectronVeto" for the PAT
  // Isolations are produced at PAT level by the PḧotonPFIsolation producer
  //edm::Handle<edm::ValueMap<bool>> hasMatchedPromptElectronHandle;
  //event.getByLabel(edm::InputTag("photonPFIsolation", "hasMatchedPromptElectron", "PAT"), hasMatchedPromptElectronHandle);

  //isValid &= ! (*hasMatchedPromptElectronHandle)[photonRef];
  isValid &= photonRef->passElectronVeto();


  if (! isValid)
    return false;
  //giulia -- use isolation stored in the pat::photon


 // // Now, isolations
 // edm::Handle<edm::ValueMap<double>> chargedHadronsIsolationHandle;
 // event.getByLabel(edm::InputTag("photonPFIsolation", "chargedHadronsIsolation", "PAT"), chargedHadronsIsolationHandle);

 // edm::Handle<edm::ValueMap<double>> neutralHadronsIsolationHandle;
 // event.getByLabel(edm::InputTag("photonPFIsolation", "neutralHadronsIsolation", "PAT"), neutralHadronsIsolationHandle);

 // edm::Handle<edm::ValueMap<double>> photonIsolationHandle;
 // event.getByLabel(edm::InputTag("photonPFIsolation", "photonIsolation", "PAT"), photonIsolationHandle);
  /*
     edm::Handle<edm::ValueMap<double>> chargedHadronsIsolationHandle;
     event.getByLabel(edm::InputTag("photonPFIsolation", "footchargediso", "PAT"), chargedHadronsIsolationHandle);

     edm::Handle<edm::ValueMap<double>> neutralHadronsIsolationHandle;
     event.getByLabel(edm::InputTag("photonPFIsolation", "footneutraliso", "PAT"), neutralHadronsIsolationHandle);

     edm::Handle<edm::ValueMap<double>> photonIsolationHandle;
     event.getByLabel(edm::InputTag("photonPFIsolation", "footphotoniso", "PAT"), photonIsolationHandle);
     */

//  isValid &= getCorrectedPFIsolation((*chargedHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 1.86;
//  isValid &= getCorrectedPFIsolation((*neutralHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (2.64 + 0.0025 * photonRef->pt());
//  isValid &= getCorrectedPFIsolation((*photonIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) < (1.2 + 0.001 * photonRef->pt());
//
  //giulia -- use isolation stored in the pat::photon

  isValid &= getCorrectedPFIsolation(photonRef->chargedHadronIso() ,rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 1.86;
  isValid &= getCorrectedPFIsolation(photonRef->neutralHadronIso(), rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (2.64 + 0.0025 * photonRef->pt());
  isValid &= getCorrectedPFIsolation(photonRef->photonIso(), rho, photonRef->eta(), IsolationType::PHOTONS) < (1.2 + 0.001 * photonRef->pt());


  return isValid;
}
//---------------------


// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012
bool GammaJetFilter::isValidPhotonEB2012(const pat::PhotonRef& photonRef, edm::Event& event) {
  if (mIsMC && !photonRef->genPhoton())
    return false;

  bool isValid = true;

  isValid &= photonRef->hadTowOverEm() < 0.05;
  isValid &= photonRef->sigmaIetaIeta() < 0.011;

  if (! isValid)
    return false;

  edm::Handle<double> rhos;
  event.getByLabel(edm::InputTag("fixedGridRhoFastjetAll"), rhos);
  double rho = *rhos;


  //giulia --  use the method " passElectronVeto" for the PAT
  // Isolations are produced at PAT level by the PḧotonPFIsolation producer
  //edm::Handle<edm::ValueMap<bool>> hasMatchedPromptElectronHandle;
  //event.getByLabel(edm::InputTag("photonPFIsolation", "hasMatchedPromptElectron", "PAT"), hasMatchedPromptElectronHandle);

  //isValid &= ! (*hasMatchedPromptElectronHandle)[photonRef];
  isValid &= photonRef->passElectronVeto();


  if (! isValid)
    return false;
  //giulia -- use isolation stored in the pat::photon


 // // Now, isolations
 // edm::Handle<edm::ValueMap<double>> chargedHadronsIsolationHandle;
 // event.getByLabel(edm::InputTag("photonPFIsolation", "chargedHadronsIsolation", "PAT"), chargedHadronsIsolationHandle);

 // edm::Handle<edm::ValueMap<double>> neutralHadronsIsolationHandle;
 // event.getByLabel(edm::InputTag("photonPFIsolation", "neutralHadronsIsolation", "PAT"), neutralHadronsIsolationHandle);

 // edm::Handle<edm::ValueMap<double>> photonIsolationHandle;
 // event.getByLabel(edm::InputTag("photonPFIsolation", "photonIsolation", "PAT"), photonIsolationHandle);
  /*
     edm::Handle<edm::ValueMap<double>> chargedHadronsIsolationHandle;
     event.getByLabel(edm::InputTag("photonPFIsolation", "footchargediso", "PAT"), chargedHadronsIsolationHandle);

     edm::Handle<edm::ValueMap<double>> neutralHadronsIsolationHandle;
     event.getByLabel(edm::InputTag("photonPFIsolation", "footneutraliso", "PAT"), neutralHadronsIsolationHandle);

     edm::Handle<edm::ValueMap<double>> photonIsolationHandle;
     event.getByLabel(edm::InputTag("photonPFIsolation", "footphotoniso", "PAT"), photonIsolationHandle);
     */

//  isValid &= getCorrectedPFIsolation((*chargedHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 0.7;
//  isValid &= getCorrectedPFIsolation((*neutralHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (0.4 + 0.04 * photonRef->pt());
//  isValid &= getCorrectedPFIsolation((*photonIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) < (0.5 + 0.005 * photonRef->pt());
//
  //giulia -- use isolation stored in the pat::photon

  isValid &= getCorrectedPFIsolation(photonRef->chargedHadronIso() ,rho, photonRef->eta(), IsolationType::CHARGED_HADRONS) < 0.7;
  isValid &= getCorrectedPFIsolation(photonRef->neutralHadronIso(), rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (0.4 + 0.04 * photonRef->pt());
  isValid &= getCorrectedPFIsolation(photonRef->photonIso(), rho, photonRef->eta(), IsolationType::PHOTONS) < (0.5 + 0.005 * photonRef->pt());


  return isValid;
}
bool GammaJetFilter::isValidPhotonEB(const pat::Photon& photon, const double rho, const EcalRecHitCollection* recHits, const CaloTopology& topology) {
  if (mIsMC && !photon.genPhoton())
    return false;

  bool isValid = ! photon.hasPixelSeed();
  isValid &= photon.hadronicOverEm() < 0.05;
  isValid &= photon.sigmaIetaIeta() < 0.011;
  isValid &= photon.trkSumPtHollowConeDR04() < (2.0 + 0.001 * photon.et() + 0.0167 * rho);
  isValid &= photon.ecalRecHitSumEtConeDR04() < (4.2 + 0.006 * photon.et() + 0.183 * rho);
  isValid &= photon.hcalTowerSumEtConeDR04() < (2.2 + 0.0025 * photon.et() + 0.062 * rho);

  if (isValid) {
    // Spike cleaning

    isValid &= photon.sigmaIetaIeta() > 0.001;

    do {
      // Get sigmaIPhiIPhi for photon
      const float w0 = 4.7; // See http://cmslxr.fnal.gov/lxr/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h#087
      if (! recHits)
      	break;

      reco::SuperClusterRef superCluster = photon.superCluster();
      if (superCluster.isNull())
      	break;

      reco::CaloClusterPtr caloCluster = superCluster->seed();
      if (caloCluster.isNull())
      	break;

      std::vector<float> covariances = EcalClusterTools::localCovariances(*caloCluster, recHits, &topology, w0);
      float sigmaIPhiIPhi = sqrt(covariances[2]);

      isValid &= sigmaIPhiIPhi > 0.001;
    } while (false);
  }

  return isValid;
}

void GammaJetFilter::readJSONFile() {
  Json::Value root;
  Json::Reader reader;
  std::ifstream file(mJSONFile.c_str());
  if (! reader.parse(file, root)) {
    throw cms::Exception("ReadError")
      << "Failed to parse luminosity JSON file '" << mJSONFile << "'" << std::endl;
  }

  mValidRuns.reset(new Json::Value(root));
}

void GammaJetFilter::readCSVFile() {
  FILE* iff = fopen(mCSVFile.c_str(), "r");

  if(iff == 0) {
    throw cms::Exception("ReadError")
      << "Failed to parse luminosity CSV file '" << mCSVFile << "'" << std::endl;
  }

  int run = 0, fill = 0;
  int lumiSection_left = 0, lumiSection_right = 0;
  double lumiRecorded = 0.;
  double true_pu = 0;


  /* lumiCalc2 format :
   * Run:Fill,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub),avgPU
   * use 'lumiCalc2.py -i lumiSummary.json -o output.csv -b stable lumibyls' to generate file
   */

  // Skip header line
  char buffer[1024];
  fgets(buffer, 1024, iff);

  while (fscanf(iff, "%d:%d,%d:%d,%*[^,],%*[^,],%*f,%*f,%lf,%lf", &run, &fill, &lumiSection_left, &lumiSection_right, &lumiRecorded, &true_pu) > 0 ) {

    if (lumiSection_right == 0)
      continue;

    mLumiByLS[std::pair<unsigned int, unsigned int>(run, lumiSection_right)] = lumiRecorded; //in mb^(-1)
    mTruePUByLS[std::pair<unsigned int, unsigned int>(run, lumiSection_right)] = true_pu; //in mb^(-1)
  }

  fclose(iff);

  assert(mLumiByLS.size() > 0);
}

void GammaJetFilter::updateLuminosity(const edm::LuminosityBlock& lumiBlock) {
  double eventLumi = mLumiByLS[std::pair<unsigned int, unsigned int>(lumiBlock.id().run(), lumiBlock.id().luminosityBlock())];
  double newLumi = mTotalLuminosity->GetVal() + eventLumi;
  mTotalLuminosity->SetVal(newLumi);
}

void GammaJetFilter::particleToTree(const reco::Candidate* particle, TTree* t, std::vector<boost::shared_ptr<void> >& addresses) {
  addresses.clear();

  addresses.push_back(boost::shared_ptr<void>(new int((particle) ? 1 : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->et() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->pt() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->eta() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->phi() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->px() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->py() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->pz() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->energy() : 0)));

  updateBranch(t, addresses[0].get(), "is_present", "I");
  updateBranch(t, addresses[1].get(), "et");
  updateBranch(t, addresses[2].get(), "pt");
  updateBranch(t, addresses[3].get(), "eta");
  updateBranch(t, addresses[4].get(), "phi");
  updateBranch(t, addresses[5].get(), "px");
  updateBranch(t, addresses[6].get(), "py");
  updateBranch(t, addresses[7].get(), "pz");
  updateBranch(t, addresses[8].get(), "e");
}


void GammaJetFilter::photonToTree(const pat::PhotonRef& photonRef, pat::Photon& photon, const edm::Event& event) {
  std::vector<boost::shared_ptr<void> > addresses;
//redefine the common variables instead of using particleTotree, because the photonhas been 
//corrected for regression etc.
 int pho_is_present =1;
 float pho_et = photon.et();
 float pho_pt = photon.pt();
 float pho_eta = photon.eta();
 float pho_phi = photon.phi();
 float pho_px = photon.px();
 float pho_py = photon.py();
 float pho_pz = photon.pz();
 float pho_e = photon.energy();
 updateBranch(mPhotonTree,&pho_is_present,"is_present","I");
 updateBranch(mPhotonTree,&pho_et,"et");
 updateBranch(mPhotonTree,&pho_pt,"pt");
 updateBranch(mPhotonTree,&pho_eta,"eta");
 updateBranch(mPhotonTree,&pho_phi,"phi"); 
 updateBranch(mPhotonTree,&pho_px,"px");
 updateBranch(mPhotonTree,&pho_py,"py");
 updateBranch(mPhotonTree,&pho_pz,"pz");
 updateBranch(mPhotonTree,&pho_e,"e");

  bool hasPixelSeed = photonRef->hasPixelSeed();
  updateBranch(mPhotonTree, &hasPixelSeed, "has_pixel_seed", "O");

  // Photon ID related
  float hadTowOverEm = photonRef->hadTowOverEm();
  updateBranch(mPhotonTree, &hadTowOverEm, "hadTowOverEm");

  float sigmaIetaIeta = photonRef->sigmaIetaIeta();
  updateBranch(mPhotonTree, &sigmaIetaIeta, "sigmaIetaIeta");

  edm::Handle<double> rhos;
  event.getByLabel(edm::InputTag("fixedGridRhoFastjetAll"), rhos);
  float rho = *rhos;
  updateBranch(mPhotonTree, &rho, "rho");

  //giulia -- taking iso from pat:Photon
  // Isolations are produced at PAT level by the PḧotonPFIsolation producer
  //edm::Handle<edm::ValueMap<bool>> hasMatchedPromptElectronHandle;
  //event.getByLabel(edm::InputTag("photonPFIsolation", "hasMatchedPromptElectron", "PAT"), hasMatchedPromptElectronHandle);

  //bool hasMatchedPromptElectron = (*hasMatchedPromptElectronHandle)[photonRef];
  
  bool hasMatchedPromptElectron ; 
  if (photonRef->passElectronVeto())
    hasMatchedPromptElectron = false;
  else 
    hasMatchedPromptElectron = true;

  updateBranch(mPhotonTree, &hasMatchedPromptElectron, "hasMatchedPromptElectron", "O");

  // Now, isolations
//  edm::Handle<edm::ValueMap<double>> chargedHadronsIsolationHandle;
//  event.getByLabel(edm::InputTag("photonPFIsolation", "chargedHadronsIsolation", "PAT"), chargedHadronsIsolationHandle);
//
//  edm::Handle<edm::ValueMap<double>> neutralHadronsIsolationHandle;
//  event.getByLabel(edm::InputTag("photonPFIsolation", "neutralHadronsIsolation", "PAT"), neutralHadronsIsolationHandle);
//
//  edm::Handle<edm::ValueMap<double>> photonIsolationHandle;
//  event.getByLabel(edm::InputTag("photonPFIsolation", "photonIsolation", "PAT"), photonIsolationHandle);
//
//giulia ---- regression integrated in 73X  
////regression energy
//  edm::Handle<edm::ValueMap<float>> regressionEnergyHandle;
//  event.getByLabel(edm::InputTag("eleNewEnergiesProducer", "energySCEleJoshPhoSemiParamV5ecorr", "PAT"),regressionEnergyHandle);
//  edm::Ptr<reco::Candidate> recoObject = photonRef->originalObjectRef();
// float regressionEnergy = (*regressionEnergyHandle)[recoObject] ;
// updateBranch(mPhotonTree, &regressionEnergy, "regressionEnergy");
// float originalEnergy = photonRef->energy();
// updateBranch(mPhotonTree, &originalEnergy, "originalEnergy");

//giulia --- comment footprint stuff

////retrieve px and py of pfcandidates to exclude from met calculation
//edm::Handle<edm::ValueMap<double>> footprintMExCorrHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintMExCorr", "PAT"), footprintMExCorrHandle);
//edm::Handle<edm::ValueMap<double>> footprintMEyCorrHandle;
//event.getByLabel(edm::InputTag("photonPFIsolation", "footprintMEyCorr", "PAT"), footprintMEyCorrHandle);
////MEx/yCorr are the quantities to compare to rawMex/y
//float footprintMExCorr = (*footprintMExCorrHandle)[photonRef]- photonRef->px();
//float footprintMEyCorr = (*footprintMEyCorrHandle)[photonRef]- photonRef->py();
//
//  updateBranch(mPhotonTree, &footprintMExCorr, "footprintMExCorr");
//  updateBranch(mPhotonTree, &footprintMEyCorr, "footprintMEyCorr");
//
//  float chargedHadronsIsolation = getCorrectedPFIsolation((*chargedHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS);
//  float neutralHadronsIsolation = getCorrectedPFIsolation((*neutralHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS);
//  float photonIsolation = getCorrectedPFIsolation((*photonIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS); 
//
  
   float chargedHadronsIsolation = getCorrectedPFIsolation(photonRef->chargedHadronIso(), rho, photonRef->eta(), IsolationType::CHARGED_HADRONS); 
   float neutralHadronsIsolation = getCorrectedPFIsolation(photonRef->neutralHadronIso(), rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS); 
   float photonIsolation = getCorrectedPFIsolation(photonRef->photonIso(), rho, photonRef->eta(), IsolationType::PHOTONS);                         
  
  
  
  
  
  updateBranch(mPhotonTree, &chargedHadronsIsolation, "chargedHadronsIsolation");
  updateBranch(mPhotonTree, &neutralHadronsIsolation, "neutralHadronsIsolation");
  updateBranch(mPhotonTree, &photonIsolation, "photonIsolation");


  mPhotonTree->Fill();

  if (mIsMC) {
    particleToTree(photonRef->genPhoton(), mPhotonGenTree, addresses);
    mPhotonGenTree->Fill();
  }
}

void GammaJetFilter::jetsToTree(const pat::Jet* firstJet, const pat::Jet* secondJet, const std::vector<TTree*>& trees) {
  jetToTree(firstJet, mIsMC, trees[0], trees[4]);
  jetToTree(secondJet, false, trees[1], trees[5]);

  // Raw jets
  const pat::Jet* rawJet = (firstJet) ? firstJet->userData<pat::Jet>("rawJet") : NULL;
  jetToTree(rawJet, false, trees[2], NULL);

  rawJet = (secondJet) ? secondJet->userData<pat::Jet>("rawJet") : NULL;
  jetToTree(rawJet, false, trees[3], NULL);
}

void findNeutrinos(const reco::Candidate* parent, std::vector<const reco::Candidate*>& neutrinos) {

  int pdg_id = abs(parent->pdgId());
  if (pdg_id == 12 || pdg_id == 14 || pdg_id == 16) {
    if (std::find_if(neutrinos.begin(), neutrinos.end(), [parent] (const reco::Candidate* candidate) -> bool {
          static double EPSILON = 0.0001;
          bool same = (candidate->pdgId() == parent->pdgId());
          same &= (fabs(candidate->px() - candidate->px()) < EPSILON);
          same &= (fabs(candidate->py() - candidate->py()) < EPSILON);
          same &= (fabs(candidate->pz() - candidate->pz()) < EPSILON);

	  return same;

	  }) == neutrinos.end()) {
      neutrinos.push_back(parent);
    }
    return;
  }

  for (unsigned int i = 0; i < parent->numberOfDaughters(); i++) {
    const reco::Candidate* d = parent->daughter(i);
    findNeutrinos(d, neutrinos);
  }
}

void GammaJetFilter::jetToTree(const pat::Jet* jet, bool _findNeutrinos, TTree* tree, TTree* genTree) {
  std::vector<boost::shared_ptr<void> > addresses;
  particleToTree(jet, tree, addresses);

  if (mIsMC) {
    mNeutrinos->Clear("C");
    mNeutrinosPDG->Clear("C");
  }

  if (jet) {
    float area = jet->jetArea();
    updateBranch(tree, &area, "jet_area");


    //giulia --- comment b-tagging  and qg tagging stuff
    // // B-Tagging
    // float tcHighEfficiency = jet->bDiscriminator("trackCountingHighEffBJetTags");
    // float tcHighPurity = jet->bDiscriminator("trackCountingHighPurBJetTags");

    // float ssvHighEfficiency = jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    // float ssvHighPurity = jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");

    // float jetProbability = jet->bDiscriminator("jetProbabilityBJetTags");
    // float jetBProbability = jet->bDiscriminator("jetBProbabilityBJetTags");

    // // New 2012
    // float csv = jet->bDiscriminator("combinedSecondaryVertexBJetTags");

    // // Quark Gluon tagging
    // float qgTagMLP = jet->userFloat("qgTagMLP");
    // float qgTagLikelihood = jet->userFloat("qgTagLikelihood");

    // updateBranch(tree, &tcHighEfficiency, "btag_tc_high_eff");
    // updateBranch(tree, &tcHighPurity, "btag_tc_high_pur");
    // updateBranch(tree, &ssvHighEfficiency, "btag_ssv_high_eff");
    // updateBranch(tree, &ssvHighPurity, "btag_ssv_high_pur");
    // updateBranch(tree, &jetProbability, "btag_jet_probability");
    // updateBranch(tree, &jetBProbability, "btag_jet_b_probability");
    // updateBranch(tree, &csv, "btag_csv");
    // updateBranch(tree, &qgTagMLP, "qg_tag_mlp");
    // updateBranch(tree, &qgTagLikelihood, "qg_tag_likelihood");

    //jet energy composition
    float jetCHEn = jet->chargedHadronEnergy();
    float jetNHEn = jet->neutralHadronEnergy();
    float jetCEEn = jet->chargedEmEnergy();
    float jetNEEn = jet->neutralEmEnergy();
    float jetPhEn = jet->photonEnergy();
    float jetElEn = jet->electronEnergy();
    float jetMuEn = jet->chargedMuEnergy();
    //jet constituents multiplicities
    int jetPhMult = jet->photonMultiplicity();
    int jetNHMult = jet->neutralHadronMultiplicity();
    int jetElMult = jet->electronMultiplicity();
    int jetCHMult = jet->chargedHadronMultiplicity();

    updateBranch(tree, &jetCHEn, "jet_CHEn");
    updateBranch(tree, &jetNHEn, "jet_NHEn");
    updateBranch(tree, &jetPhEn, "jet_PhEn");
    updateBranch(tree, &jetElEn, "jet_ElEn");
    updateBranch(tree, &jetMuEn, "jet_MuEn");
    updateBranch(tree, &jetCEEn, "jet_CEEn");
    updateBranch(tree, &jetNEEn, "jet_NEEn");
    updateBranch(tree, &jetPhMult, "jet_PhMult", "I");
    updateBranch(tree, &jetNHMult, "jet_NHMult", "I");
    updateBranch(tree, &jetElMult, "jet_ElMult", "I");
    updateBranch(tree, &jetCHMult, "jet_CHMult", "I");

    tree->Fill(); // This Fill() must be called inside the {} block, otherwise it'll crash. Don't move it!
  } else {
    tree->Fill();
  }

  if (genTree) {
    particleToTree((jet) ? jet->genJet() : NULL, genTree, addresses);

    if (_findNeutrinos) {
      static bool init = false;
      if (! init) {
        genTree->Branch("neutrinos", &mNeutrinos, 32000, 0);
        genTree->Branch("neutrinos_pdg_id", &mNeutrinosPDG, 32000, 0);
        init = true;
      }
    }

    if (jet && _findNeutrinos) {
      const reco::Candidate* parton = (jet) ? jet->genParton() : NULL;

      if (parton) {
        if (abs(parton->pdgId()) == 5 || abs(parton->pdgId()) == 4) {

          std::vector<const reco::Candidate*> neutrinos;
          findNeutrinos(parton, neutrinos);

          if (neutrinos.size() > 0) {
            // Build TCloneArray of TLorentzVector
            unsigned int index = 0;
            for (const reco::Candidate* neutrino: neutrinos) {
              TLorentzVector* p4 = (TLorentzVector*) mNeutrinos->ConstructedAt(index);
              p4->SetPxPyPzE(neutrino->px(), neutrino->py(), neutrino->pz(), neutrino->energy());

              TParameter<int>* pdg_id = (TParameter<int>*) mNeutrinosPDG->ConstructedAt(index++);
              pdg_id->SetVal(neutrino->pdgId());
            }
          }
        }

      }
    }

    // Add parton id and pt
    const reco::Candidate* parton = (jet) ? jet->genParton() : NULL;
    int pdgId = (parton) ? parton->pdgId() : 0;
    updateBranch(genTree, &pdgId, "parton_pdg_id", "I");

    TLorentzVector parton_p4;
    if (parton) {
      parton_p4.SetPxPyPzE(parton->px(), parton->py(), parton->pz(), parton->energy());
    }
    TLorentzVector* p_parton_p4 = &parton_p4;

    TBranch* branch = genTree->GetBranch("parton_p4");
    if (branch == NULL) {
      branch = genTree->Branch("parton_p4", &p_parton_p4);
    } else {
      branch->SetAddress(&p_parton_p4);
    }

    int flavour = (jet) ? jet->partonFlavour() : 0;
    updateBranch(genTree, &flavour, "parton_flavour", "I");

    genTree->Fill();
  }
}

void GammaJetFilter::metsToTree(const pat::MET& met, const pat::MET& rawMet, const std::vector<TTree*>& trees) {
  metToTree(&met, trees[0], trees[2]);
  metToTree(&rawMet, trees[1], NULL);
}

void GammaJetFilter::metToTree(const pat::MET* met, TTree* tree, TTree* genTree) {
  std::vector<boost::shared_ptr<void> > addresses;
  particleToTree(met, tree, addresses);

  tree->Fill();

  if (genTree) {
    particleToTree((met) ? met->genMET() : NULL, genTree, addresses);
    genTree->Fill();
  }
}

void GammaJetFilter::electronsToTree(const edm::Handle<pat::ElectronCollection>& electrons, const reco::Vertex& pv) {

  int n = electrons->size();
  static int   id[30];
  static float isolation[30];
  static float pt[30];
  static float px[30];
  static float py[30];
  static float pz[30];
  static float eta[30];
  static float phi[30];
  static int   charge[30];

  int i = 0;
  for (pat::ElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); ++it, i++) {
    const pat::Electron& electron = *it;

    if (i >= 30)
      break;

    // See https://twiki.cern.ch/twiki/bin/view/CMS/TopLeptonPlusJetsRefSel_el
    bool elecID = fabs(pv.z() - it->vertex().z()) < 1.;
    elecID     &= it->et() > 30.;
    elecID     &= fabs(it->eta()) < 2.5 && (it->superCluster()->eta() > 1.4442 && it->superCluster()->eta() < 1.5660);
    elecID     &= it->dB() < 0.02;
    elecID     &= ((int) it->electronID("eidLoose") & 0x1);

    float iso     = (it->dr03TkSumPt() + it->dr03EcalRecHitSumEt() + it->dr03HcalTowerSumEt()) / it->et();

    id[i]         = elecID;
    isolation[i]  = iso;
    pt[i]         = electron.pt();
    px[i]         = electron.px();
    py[i]         = electron.py();
    pz[i]         = electron.pz();
    eta[i]        = electron.eta();
    phi[i]        = electron.phi();
    charge[i]     = electron.charge();
  }

  updateBranch(mElectronsTree, &n, "n", "I");
  updateBranchArray(mElectronsTree, id, "id", "n", "I");
  updateBranchArray(mElectronsTree, isolation, "isolation", "n");
  updateBranchArray(mElectronsTree, pt, "pt", "n");
  updateBranchArray(mElectronsTree, px, "px", "n");
  updateBranchArray(mElectronsTree, py, "py", "n");
  updateBranchArray(mElectronsTree, pz, "pz", "n");
  updateBranchArray(mElectronsTree, eta, "eta", "n");
  updateBranchArray(mElectronsTree, phi, "phi", "n");
  updateBranchArray(mElectronsTree, charge, "charge", "n", "I");

  mElectronsTree->Fill();
}

void GammaJetFilter::muonsToTree(const edm::Handle<pat::MuonCollection>& muons, const reco::Vertex& pv) {

  int n = muons->size();
  static int   id[30];
  static float isolation[30];
  static float delta_beta_isolation[30];
  static float pt[30];
  static float px[30];
  static float py[30];
  static float pz[30];
  static float eta[30];
  static float phi[30];
  static int   charge[30];

  int i = 0;
  for (pat::MuonCollection::const_iterator it = muons->begin(); it != muons->end(); ++it, i++) {
    const pat::Muon& muon = *it;

    if (i >= 30)
      break;

    // See https://twiki.cern.ch/twiki/bin/view/CMS/TopLeptonPlusJetsRefSel_mu
    bool muonID = it->isGlobalMuon();
    //FIXME: reco::Tracks need to be keept in PF2PAT.
    //It's not the case right now, so muon ID will be incorrect
    if (it->globalTrack().isNull() || it->innerTrack().isNull() || it->muonBestTrack().isNull() || it->track().isNull()) {
      muonID = false;
    } else {
      muonID     &= it->globalTrack()->normalizedChi2() < 10.;
      muonID     &= it->globalTrack()->hitPattern().numberOfValidMuonHits() > 0;
      muonID     &= it->numberOfMatchedStations() > 1;
      muonID     &= it->dB() < 0.2;
      muonID     &= fabs(it->muonBestTrack()->dz(pv.position())) < 0.5;
      muonID     &= it->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
      muonID     &= it->track()->hitPattern().trackerLayersWithMeasurement() > 5;
    }

    float relIso = (it->chargedHadronIso() + it->neutralHadronIso() + it->photonIso()) / it->pt();
    float deltaBetaRelIso = (it->chargedHadronIso() + std::max((it->neutralHadronIso() + it->photonIso()) - 0.5 * it->puChargedHadronIso(), 0.0)) / it->pt();

    id[i]          = muonID;
    isolation[i]   = relIso;
    delta_beta_isolation[i] = deltaBetaRelIso;
    pt[i]          = muon.pt();
    px[i]          = muon.px();
    py[i]          = muon.py();
    pz[i]          = muon.pz();
    eta[i]         = muon.eta();
    phi[i]         = muon.phi();
    charge[i]      = muon.charge();
  }

  updateBranch(mMuonsTree, &n, "n", "I");
  updateBranchArray(mMuonsTree, id, "id", "n", "I");
  updateBranchArray(mMuonsTree, isolation, "relative_isolation", "n");
  updateBranchArray(mMuonsTree, delta_beta_isolation, "delta_beta_relative_isolation", "n");
  updateBranchArray(mMuonsTree, pt, "pt", "n");
  updateBranchArray(mMuonsTree, px, "px", "n");
  updateBranchArray(mMuonsTree, py, "py", "n");
  updateBranchArray(mMuonsTree, pz, "pz", "n");
  updateBranchArray(mMuonsTree, eta, "eta", "n");
  updateBranchArray(mMuonsTree, phi, "phi", "n");
  updateBranchArray(mMuonsTree, charge, "charge", "n", "I");

  mMuonsTree->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(GammaJetFilter);
