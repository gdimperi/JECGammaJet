import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = cms.string("START53_V7G::All")
process.GlobalTag.globaltag = cms.string("START53_V27::All")

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
readFiles = cms.untracked.vstring(
    )

#readFiles.extend( [                                                                                                                                                        
#       '/store/user/sbrochet/G_Pt-170to300_TuneZ2star_8TeV_pythia6/G_Pt-170to300_START53_V7A_22Feb13-v1/31346b79deb97ac1b786d692cd650a21/patTuple_PF2PAT_MC_10_3_15V.root',
readFiles.extend( [                                                                                                                                                        
       #'file:/afs/cern.ch/work/g/gdimperi/CMSSW_5_3_9_patch2/src/JetMETCorrections/GammaJetFilter/crab/patTuple_PF2PAT_MC_G_Pt-30to50.root',

       #'file:/afs/cern.ch/work/g/gdimperi/CMSSW_5_3_9_patch2/src/JetMETCorrections/GammaJetFilter/crab/patTuple_PF2PAT_MC_G_Pt-50to80__step1.root',
       #'file:/afs/cern.ch/work/g/gdimperi/CMSSW_5_3_9_gammajet/src/JetMETCorrections/GammaJetFilter/crab/patTuple_PF2PAT_MC_G_Pt-50to80__step0-1_START53_V21.root',
       'file:../../crab/patTuple_PF2PAT_MC_G-Pt120to170.root'
       #'file:/afs/cern.ch/work/g/gdimperi/CMSSW_5_3_9_patch2/src/JetMETCorrections/GammaJetFilter/crab/patTuple_PF2PAT_MC_G_Pt-50to80__step1.root',

       ])

process.source = cms.Source ("PoolSource", fileNames = readFiles)

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
          #tag    = cms.string('JetCorrectorParametersCollection_Winter14_V1_DATA_AK5PF'),
          tag    = cms.string('JetCorrectorParametersCollection_Winter14_V1_MC_AK5PFchs'),
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



processedEvents = int(options.processedEvents) if isinstance(options.processedEvents, int) and int(options.processedEvents) != 0 else 1
crossSection = float(options.crossSection) if isinstance(options.crossSection, float) and float(options.crossSection) != 0 else 1
ptHatMin = options.lowPtHat if isinstance(options.lowPtHat, float) else -1
ptHatMax = options.highPtHat if isinstance(options.highPtHat, float) else -1

print("Running on sample with:")
print("\tNumber of processed events: %d" % processedEvents)
print("\tCross-section: %f" % crossSection)
print("\tPt hat min: %f" % ptHatMin)
print("\tPt hat max: %f" % ptHatMax)

process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(True),
    photons = cms.untracked.InputTag("selectedPatPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    crossSection = cms.double(crossSection),
    generatedEvents = cms.uint64(processedEvents),
    ptHatMin = cms.untracked.double(ptHatMin),
    ptHatMax = cms.untracked.double(ptHatMax),

    dumpAllGenParticles = cms.untracked.bool(False),

    runOnNonCHS   = cms.untracked.bool(True),
    runOnCHS      = cms.untracked.bool(True),


    runOnPFAK5    = cms.untracked.bool(True),
    runOnPFAK7    = cms.untracked.bool(False),
    runOnCA8      = cms.untracked.bool(False),

    runOnCaloAK5  = cms.untracked.bool(False),
    runOnCaloAK7  = cms.untracked.bool(False),

    # JEC
    doJetCorrection = cms.untracked.bool(True),
    correctJecFromRaw = cms.untracked.bool(True),

    correctorLabel = cms.untracked.string("ak5PFchsL1FastL2L3"),
    #correctorLabel = cms.untracked.string("ak7PFchsL1FastL2L3"),                                                            

    #correctorLabel = cms.untracked.string("ak5PFResidual")

    # MET
    redoTypeIMETCorrection = cms.untracked.bool(True)
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
    #fileName = cms.string("output_mc_Gj-Pt120to170.root")
    )

#process.out.fileName = 'patTuple_cleaned.root'
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#giulia
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
#process.outpath = cms.EndPath(process.out)
