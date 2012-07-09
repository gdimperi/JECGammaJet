import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string("GR_R_52_V9D::All") ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#process.GlobalTag.globaltag = cms.string("GR_R_52_V7::All") ## 2011 JEC
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:patTuple_PF2PAT.root'
        #'/store/user/sbrochet/Photon/JetMet_PF2PAT_Run2012A_PromptReco_09June/292d73a6795d0a493723b2d3a624156c/patTuple_PF2PAT_84_2_tpP.root'        
        '/store/user/sbrochet/SinglePhoton/JetMet_PF2PAT_Run2012B_PromptReco_23June/e316a04ff768fa1f713908f8996bbb95/patTuple_PF2PAT_216_1_xRU.root'
    )
)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('datasetName',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "The dataset currently processed. A folder named 'datasetName' must exists")

options.parseArguments()

fullPath = os.path.join(os.getcwd(), options.datasetName)

process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(False),
    photons = cms.untracked.InputTag("selectedPatPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    json = cms.string(os.path.join(fullPath, "lumiSummary.json")),
    csv = cms.string(os.path.join(fullPath, "lumibyls.csv")),
    filterData = cms.untracked.bool(True),

    runOnNonCHS   = cms.untracked.bool(False),
    runOnCHS      = cms.untracked.bool(True),

    runOnPFAK5    = cms.untracked.bool(True),
    runOnPFAK7    = cms.untracked.bool(False),

    runOnCaloAK5  = cms.untracked.bool(False),
    runOnCaloAK7  = cms.untracked.bool(False),

    # JEC
    doJetCorrection = cms.untracked.bool(True),
    correctJecFromRaw = cms.untracked.bool(True),
    correctorLabel = cms.untracked.string("ak5PFchsL1FastL2L3")
    #correctorLabel = cms.untracked.string("ak5PFResidual")
    )

process.p = cms.Path(process.gammaJet)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("output.root"),
#    SelectEvents = cms.untracked.PSet(
#      SelectEvents = cms.vstring('p')
#      )
#    )

#process.out.outputCommands = cms.untracked.vstring('keep *',
#    'drop *_selectedPatJets*_*_*',
#    'drop *_selectedPatPhotons*_*_*',
#    #'keep *_selectedPatJets*_genJets_*',
#    'keep *_selectedPatJets*_caloTowers_*',
#    # Drop CHS
#    'drop *_*chs*_*_*'
#)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string("output.root")
     )

#process.out.fileName = 'patTuple_cleaned.root'
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
