import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")


# Do some CHS stuff
process.ak5PFchsL1Fastjet  = process.ak5PFL1Fastjet.clone(algorithm = 'AK5PFchs')
process.ak5PFchsL2Relative = process.ak5PFL2Relative.clone(algorithm = 'AK5PFchs')
process.ak5PFchsL3Absolute = process.ak5PFL3Absolute.clone(algorithm = 'AK5PFchs')
process.ak5PFchsResidual   = process.ak5PFResidual.clone(algorithm = 'AK5PFchs')
process.ak5PFchsL1FastL2L3 = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak5PFchsL1Fastjet', 'ak5PFchsL2Relative','ak5PFchsL3Absolute')
    )
process.ak5PFchsL1FastL2L3Residual = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak5PFchsL1Fastjet', 'ak5PFchsL2Relative','ak5PFchsL3Absolute','ak5PFchsResidual')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      #'/store/user/sbrochet/Photon/Photon_Run2012A-recover-06Aug2012_18Feb13-v1/298629d7efe53cfa022ca63c838ed612/patTuple_PF2PAT_16_1_Z3z.root'
      #'file:patTuple_PF2PAT.root'
      #'file:/afs/cern.ch/work/g/gdimperi/GammaJet/CMSSW_5_3_14/src/JetMETCorrections/GammaJetFilter/analysis/2ndLevel/patTuple_PF2PAT_182_1_cUM.root'
      'file:/afs/cern.ch/work/g/gdimperi/GammaJet/CMSSW_5_3_14/src/JetMETCorrections/GammaJetFilter/crab/SiglePhotonParked_runD.root'
    )
)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('datasetName',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "The dataset currently processed. A folder named 'datasetName' must exists")

options.register ('globalTag',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "The globaltag to use")

options.parseArguments()
if len(options.globalTag) == 0:
  raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")

process.GlobalTag.globaltag = cms.string("%s::All" % options.globalTag) ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#process.GlobalTag.globaltag = cms.string("FT_53_V21_AN5::All") ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)

fullPath = os.path.join(os.getcwd(), options.datasetName)



##############################################
#                                            #
#      connect to a local sqlite file        #
#                                            # 
##############################################

process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
       messageLevel = cms.untracked.int32(0)
       ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
          record = cms.string('JetCorrectionsRecord'),
          tag    = cms.string('JetCorrectorParametersCollection_Winter14_V1_DATA_AK5PFchs'),
          #tag    = cms.string('JetCorrectorParametersCollection_Winter14_V1_MC_AK5PFchs'),
          label  = cms.untracked.string('AK5PFchs')
          ),

      ##..................................................
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file
      ),
      #connect = cms.string('sqlite:Winter14_V1_DATA.db')
      # uncomment above tag lines and this comment to use MC JEC
      connect = cms.string('sqlite:Winter14_V1_MC.db')
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')




process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(False),
    photons = cms.untracked.InputTag("selectedPatPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    #json = cms.string(os.path.join(fullPath, "lumiSummary.json")),
    #csv = cms.string(os.path.join(fullPath, "lumibyls.csv")),
    json = cms.string("lumiSummary.json"),
    csv = cms.string("lumibyls.csv"),     

    filterData = cms.untracked.bool(True),

    runOnNonCHS   = cms.untracked.bool(True),
    runOnCHS      = cms.untracked.bool(True),

    runOnPFAK5    = cms.untracked.bool(True),
    runOnPFAK7    = cms.untracked.bool(False),
    runOnPFCA8    = cms.untracked.bool(False),

    runOnCaloAK5  = cms.untracked.bool(False),
    runOnCaloAK7  = cms.untracked.bool(False),

    # JEC
    doJetCorrection = cms.untracked.bool(True),
    correctJecFromRaw = cms.untracked.bool(True),
    correctorLabel = cms.untracked.string("ak5PFchsL1FastL2L3"),
    #correctorLabel = cms.untracked.string("ak5PFchsL1FastL2L3Residual"),
    #correctorLabel = cms.untracked.string("ak7PFchsL1FastL2L3"),                                                            

    # MET
    redoTypeIMETCorrection = cms.untracked.bool(True)
    )

process.p = cms.Path(process.gammaJet)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string("output.root")
     #fileName = cms.string("output_data_SiglePhotonParked_runD.root")
     )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
