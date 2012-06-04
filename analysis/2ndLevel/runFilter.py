import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string("GR_R_52_V7::All") ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:patTuple_PF2PAT.root'
        '/store/user/sbrochet/Photon/JetMet_PF2PAT_Run2012A_PromptReco_22May//292d73a6795d0a493723b2d3a624156c/patTuple_PF2PAT_102_3_7Kw.root'
        
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

process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")

process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(False),
    photons = cms.untracked.InputTag("selectedPatPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    json = cms.string(os.path.join(fullPath, "lumiSummary.json")),
    csv = cms.string(os.path.join(fullPath, "lumibyls.csv")),
    filterData = cms.untracked.bool(True),

    # JEC
    doJetCorrection = cms.untracked.bool(False),
    correctJecFromRaw = cms.untracked.bool(False),
    #correctorLabel = cms.untracked.string("ak5PFL1FastL2L3")
    correctorLabel = cms.untracked.string("ak5PFResidual")
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
