import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = cms.string("START53_V23::All")
process.GlobalTag.globaltag = cms.string("PHYS14_25_V2::All")

process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")

process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))


# Do some CHS stuff
process.ak4PFchsL1Fastjet  = process.ak4PFL1Fastjet.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL2Relative = process.ak4PFL2Relative.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL3Absolute = process.ak4PFL3Absolute.clone(algorithm = 'AK4PFchs')
process.ak4PFchsResidual   = process.ak4PFResidual.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL1FastL2L3 = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak4PFchsL1Fastjet', 'ak4PFchsL2Relative','ak4PFchsL3Absolute')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

from FWCore.ParameterSet.VarParsing import VarParsing
#readFiles = cms.untracked.vstring(
#    )

#readFiles.extend( [
#  'file:/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_3_2/test/test_file_MINIAOD_for_JEC2015.root',
#  ])

process.source = cms.Source (
    "PoolSource", 
    fileNames = cms.untracked.vstring('file:/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_3_2/test/test_file_MINIAOD_for_JEC2015.root')
    )

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('processedEvents',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "The number of processed events")

options.register ('crossSection',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.float,
				  "Dataset cross section")

options.register ('lowPtHat',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.float,
				  "Min. generated pt (-1 if the sample is unbinned)")

options.register ('highPtHat',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.float,
				  "Max. generated pt (-1 if the sample is unbinned)")

options.parseArguments()

processedEvents = int(options.processedEvents) if isinstance(options.processedEvents, int) and int(options.processedEvents) != 0 else 1
crossSection = float(options.crossSection) if isinstance(options.crossSection, float) and float(options.crossSection) != 0 else 1
ptHatMin = options.lowPtHat if isinstance(options.lowPtHat, float) else -1
ptHatMax = options.highPtHat if isinstance(options.highPtHat, float) else -1

print("Running on sample with:")
print("\tNumber of processed events: %d" % processedEvents)
print("\tCross-section: %f" % crossSection)
print("\tPt hat min: %f" % ptHatMin)
print("\tPt hat max: %f" % ptHatMax)


## Add our PhotonIsolationProducer to the analysisSequence. This producer compute pf isolations  for our photons
#process.photonPFIsolation = cms.EDProducer("PhotonIsolationProducer",
#    src = cms.InputTag("slimmedPhotons")
#    )


process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(True),
    photons = cms.untracked.InputTag("slimmedPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    crossSection = cms.double(crossSection),
    generatedEvents = cms.uint64(processedEvents),
    ptHatMin = cms.untracked.double(ptHatMin),
    ptHatMax = cms.untracked.double(ptHatMax),

    dumpAllGenParticles = cms.untracked.bool(False),

    runOnNonCHS   = cms.untracked.bool(False),
    runOnCHS      = cms.untracked.bool(True),

    runOnPFAK4    = cms.untracked.bool(True),
    runOnPFAK8    = cms.untracked.bool(False),

    runOnCaloAK4  = cms.untracked.bool(False),
    runOnCaloAK8  = cms.untracked.bool(False),

    # JEC
    doJetCorrection = cms.untracked.bool(True),
    correctJecFromRaw = cms.untracked.bool(True),
    correctorLabel = cms.untracked.string("ak4PFchsL1FastL2L3"),
    #correctorLabel = cms.untracked.string("ak4PFResidual")

    # MET
    redoTypeIMETCorrection = cms.untracked.bool(False)
    )

process.p = cms.Path(
    process.chs *
    #process.photonPFIsolation*
    process.gammaJet)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("delete_me.root"),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
      )
    )

process.out.outputCommands = cms.untracked.vstring('keep *',
#    'drop *_selectedPatJets*_*_*',
#    'drop *_selectedPatPhotons*_*_*',
#    'keep *_selectedPatJets*_genJets_*',
#    'keep *_selectedPatJets*_caloTowers_*',
#    # Drop CHS
#    'drop *_*chs*_*_*'
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_mc.root")
    )

#process.out.fileName = 'patTuple_cleaned.root'
# set True if you want long output
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.outpath = cms.EndPath(process.out)
