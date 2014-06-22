import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")

process.load('Calibration.EleNewEnergiesProducer.elenewenergiesproducer_cfi')

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
      'file:/cmshome/gdimperi/GammaJet/test2/CMSSW_5_3_16_patch1/src/JetMETCorrections/GammaJetFilter/crab/patTuple_PF2PAT.root'
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
                  'FT53_V21A_AN6', #Winter13 2012 A, B, C, D datasets re-reco with CMSSW_5_3_7_patch6 
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "The globaltag to use")

options.parseArguments()
if len(options.globalTag) == 0:
    raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")

process.GlobalTag.globaltag = cms.string("%s::All" % options.globalTag) ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)


fullPath = os.path.join(os.getcwd(), options.datasetName)


process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(False),
    photons = cms.untracked.InputTag("selectedPatPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    #json = cms.string(os.path.join(fullPath, "lumiSummary.json")),
    #csv = cms.string(os.path.join(fullPath, "lumibyls.csv")),
    json = cms.string("Photon_Run2012A-22Jan2013/lumiSummary.json"),
    csv = cms.string("Photon_Run2012A-22Jan2013/lumibyls.csv"),     

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
    #correctorLabel = cms.untracked.string("ak5PFchsL1FastL2L3"),
    correctorLabel = cms.untracked.string("ak5PFchsL1FastL2L3Residual"),
    #correctorLabel = cms.untracked.string("ak7PFchsL1FastL2L3"),                                                            

    # MET
    redoTypeIMETCorrection = cms.untracked.bool(True)
    )

mySequence = cms.Sequence(process.eleNewEnergiesProducer * process.gammaJet)
process.p = cms.Path(mySequence)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string("output.root")
     #fileName = cms.string("output_data_SiglePhotonParked_runD.root")
     )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
