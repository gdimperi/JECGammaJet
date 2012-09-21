import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = cms.string("START53_V11::All")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from FWCore.ParameterSet.VarParsing import VarParsing
maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring(
    #"/store/user/sbrochet/G_Pt-30to50_TuneZ2_7TeV_pythia6/JetMet_PF2PAT_01apr_Fall11/5105f901d1c8d136e9eaa25f22cecd3a/patTuple_PF2PAT_MC_54_3_zYV.root"
    #"file:patTuple_PF2PAT_MC.root"
    #"/store/user/sbrochet/G_Pt-170to300_TuneZ2_7TeV_pythia6/JetMet_PF2PAT_01apr_Fall11/5105f901d1c8d136e9eaa25f22cecd3a/patTuple_PF2PAT_MC_55_1_S8q.root"
    #"/store/user/sbrochet/G_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/JetMet_PF2PAT_2012_27May_Summer12/fbe304f4e69d87969ac8d25e6b0621a8/patTuple_PF2PAT_MC_46_1_l65.root"
    "/store/user/sbrochet/G_Pt-15to30_TuneZ2star_8TeV_pythia6/JetMet_PF2PAT_2012_11July2012_Summer12/a569afd3137dcbe87038077d06761aba/patTuple_PF2PAT_MC_51_1_BkA.root"
    )
process.source = cms.Source ("PoolSource", fileNames = readFiles)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('datasetName',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "The dataset currently processed. A folder named 'datasetName' must exists")

options.parseArguments()

# Extract number of generated events, cross section, and pt hat
if os.path.exists(options.datasetName): 
  file = os.path.join(options.datasetName, "xsection");
  f = open(file, "r");
  (generatedEvents, crossSection, ptHatMax, ptHatMin) = f.read().split()
  f.close()
else:
  raise IOError("A folder named '%s' must exists", options.datasetName)

generatedEvents = int(generatedEvents)
crossSection = float(crossSection)
ptHatMax = float(ptHatMax)
ptHatMin = float(ptHatMin)

print("Running on sample with:")
print("\tNumber of generated events: %d" % generatedEvents)
print("\tCross-section: %f" % crossSection)
print("\tPt hat min: %f" % ptHatMin)
print("\tPt hat max: %f" % ptHatMax)

process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(True),
    photons = cms.untracked.InputTag("selectedPatPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    crossSection = cms.double(crossSection),
    generatedEvents = cms.uint64(int(generatedEvents)),
    ptHatMin = cms.untracked.double(ptHatMin),
    ptHatMax = cms.untracked.double(ptHatMax),

    runOnNonCHS   = cms.untracked.bool(True),
    runOnCHS      = cms.untracked.bool(True),

    runOnPFAK5    = cms.untracked.bool(True),
    runOnPFAK7    = cms.untracked.bool(False),

    runOnCaloAK5  = cms.untracked.bool(False),
    runOnCaloAK7  = cms.untracked.bool(False),

    # JEC
    doJetCorrection = cms.untracked.bool(False),
    correctJecFromRaw = cms.untracked.bool(True),
    correctorLabel = cms.untracked.string("ak5PFchsL1FastL2L3")
    #correctorLabel = cms.untracked.string("ak5PFResidual")
    )

process.p = cms.Path(process.gammaJet)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("delete_me.root"),
#    SelectEvents = cms.untracked.PSet(
#      SelectEvents = cms.vstring('p')
#      )
#    )

#process.out.outputCommands = cms.untracked.vstring('keep *',
#    'drop *_selectedPatJets*_*_*',
#    'drop *_selectedPatPhotons*_*_*',
#    'keep *_selectedPatJets*_genJets_*',
#    'keep *_selectedPatJets*_caloTowers_*',
#    # Drop CHS
#    'drop *_*chs*_*_*'
#)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_mc.root")
    )

#process.out.fileName = 'patTuple_cleaned.root'
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#process.outpath = cms.EndPath(process.out)
