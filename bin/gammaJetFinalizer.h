#pragma once

#include <TRandom3.h>

#include "Tree/AnalysisTree.h"
#include "Tree/PhotonTree.h"
#include "Tree/JetTree.h"
#include "Tree/GenJetTree.h"
#include "Tree/METTree.h"
#include "Tree/MiscTree.h"
#include "Tree/ElectronTree.h"
#include "Tree/MuonTree.h"

#include "etaBinning.h"
#include "ptBinning.h"
#include "vertexBinning.h"
#include "extrapBinning.h"
#include "newExtrapBinning.h"
#include "triggers.h"
#include "GaussianProfile.h"

#include <vector>
#include <memory>
#include <unordered_map>

namespace fwlite {
  class TFileService;
}

namespace pat {
  class Photon;
  class Jet;
}

namespace edm {
  class EventBase;
}

class TTree;
class TChain;
class TFileDirectory;


//edit giulia---------- Add CA8
enum JetAlgo {
  AK5,
  AK7,
  CA8

};

enum JetType {
  PF,
  CALO
};

template <typename T>
struct ExtrapolationVectors {
  typedef std::vector<std::vector<std::vector<T*> > > type;
};


class PUReweighter;

class GammaJetFinalizer
{
  public:
    GammaJetFinalizer();
    ~GammaJetFinalizer();

    void setInputFiles(const std::vector<std::string>& files) {
      mInputFiles = files;
      checkInputFiles();
    }

    void setDatasetName(const std::string& name) {
      mDatasetName = name;
    }

    void setJetAlgo(const std::string& jetType, const std::string& jetAlgo) {
      if (jetType == "pf") {
        mJetType = PF;
      } else {
        mJetType = CALO;
      }


      //edit giulia
      if (jetAlgo == "ak5") {
        mJetAlgo = AK5;
      } else if (jetAlgo == "ak7") {
        mJetAlgo = AK7;
      } else if (jetAlgo == "ca8"){
	mJetAlgo = CA8;

      }
    }

    void setMC(bool isMC) {
      mIsMC = isMC;
    }

    void setMCComparison(bool mcComparison) {
      mDoMCComparison = mcComparison;
    }

    void setUseExternalJEC(bool useExternalJEC) {
      mUseExternalJECCorrecion = useExternalJEC;
    }

    void setBatchJob(int currentJob, int totalJobs) {
      mIsBatchJob = true;
      mCurrentJob = currentJob;
      mTotalJobs = totalJobs;
    }

    void setAlphaCut(float alphaCut) {
      mAlphaCut = alphaCut;
    }

    void setCHS(bool chs) {
      mUseCHS = chs;
    }

    void setVerbose(bool verbose) {
      mVerbose = verbose;
    }

    void setUncutTrees(bool uncutTrees) {
      mUncutTrees = uncutTrees;
    }

    void runAnalysis();

  private:
    void checkInputFiles();
    void loadFiles(TChain& chain);

    //bool passTrigger(const TRegexp& regexp) const;
    int checkTrigger(std::string& passedTrigger, float& weight);

    void cleanTriggerName(std::string& trigger);
    void computePUWeight(const std::string& passedTrigger);

    template<typename T>
      std::vector<std::vector<T*> > buildEtaPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<T*> buildPtVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<T*> buildPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<std::vector<T*> > buildEtaVertexVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<T*> buildVertexVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<std::vector<std::vector<T*> > > buildExtrapolationEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<std::vector<T*> > buildExtrapolationVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);

    std::shared_ptr<GaussianProfile> buildNewExtrapolationVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    std::vector<std::shared_ptr<GaussianProfile>> buildNewExtrapolationEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);

    void cloneTree(TTree* from, TTree*& to);

    std::string buildPostfix();

    // Datas from step 2
    AnalysisTree analysis;
    PhotonTree photon;
    GenTree genPhoton;
    MuonTree muons;
    ElectronTree electrons;

    JetTree firstJet;
    JetTree firstRawJet;
    GenJetTree firstGenJet;

    JetTree secondJet;
    JetTree secondRawJet;
    GenJetTree secondGenJet;

    METTree MET;
    GenTree genMET;
    METTree rawMET;

    MiscTree misc;

    EtaBinning mEtaBinning;
    PtBinning mPtBinning;
    VertexBinning mVertexBinning;
    ExtrapBinning mExtrapBinning;
    NewExtrapBinning mNewExtrapBinning;

    std::vector<std::string> mInputFiles;
    std::string mDatasetName;
    JetType mJetType;
    JetAlgo mJetAlgo;
    bool mNoPUReweighting;
    bool mIsMC;
    bool mIsBatchJob;
    int mTotalJobs;
    int mCurrentJob;
    bool mUseExternalJECCorrecion;

    float  mAlphaCut;
    bool   mDoMCComparison;
    bool   mUseCHS;
    bool   mVerbose;
    bool   mUncutTrees;

    std::unordered_map<std::string, boost::shared_ptr<PUReweighter>> mLumiReweighting;
    float mPUWeight;

    float respMPF;

    // Balancing
    float respBalancing;
    float respBalancingGen;
    float respBalancingRaw;
    float respBalancingRawGen;

    // For DATA/MC comparison
    float respGenPhoton;
    float respGenGamma;
    float respPhotonGamma;

    // Triggers on data
    Triggers* mTriggers;
    MCTriggers* mMCTriggers;
    TRandom3 mRandomGenerator;
};
