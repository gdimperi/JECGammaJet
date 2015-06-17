import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

############################################################################
###### Noise Filters -- load here and apply in process path or before ######
############################################################################
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters

#process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi') #DOES NOT WORK
# Error message:
#Principal::getByToken: Found zero products matching all criteria
#Looking for type: reco::BeamHaloSummary
#Looking for module label: BeamHaloSummary

#process.load('RecoMET.METFilters.eeBadScFilter_cfi') #DOES NOT WORK
# Error message:
#   [2] Calling event method for module EEBadScFilter/'eeBadScFilter'
#Exception Message:
#Principal::getByToken: Found zero products matching all criteria
#Looking for type: edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >
#Looking for module label: reducedEcalRecHitsEE

#process.load('RecoMET.METFilters.trackingFailureFilter_cfi') #DOES NOT WORK
#process.goodVertices = cms.EDFilter(
#  "VertexSelector",
#  filter = cms.bool(False),
#  src = cms.InputTag("offlinePrimaryVertices"),
#  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
#)
# Error message:
#   [2] Calling event method for module VertexSelector/'goodVertices'
#Exception Message:
#Principal::getByToken: Found zero products matching all criteria
#Looking for type: std::vector<reco::Vertex>
#Looking for module label: offlinePrimaryVertice

process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = 'START53_V27::All'
#process.GlobalTag.globaltag = 'START70_V6::All'
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'
#process.GlobalTag.globaltag = 'POSTLS170_V7::All'
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All'
process.GlobalTag.globaltag = THISGLOBALTAG


#--------------------- Report and output ---------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.TFileService=cms.Service("TFileService",
                                 #fileName=cms.string('dijetTree_signal_M1000.root'),
                                 #fileName=cms.string('dijetTree_signal_M8000.root'),
                                 #fileName=cms.string('dijetTree_QstarToJJ_M_3000_PHYS14.root'),
                                 #fileName=cms.string('dijetTree_dataTest.root'),
                                 fileName=cms.string(THISROOTFILE),
                                 closeFileFast = cms.untracked.bool(True)
                                 )

## --- suppress long output ---> wantSummary = cms.untracked.bool(False) 

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)

############## output  edm format ###############
process.out = cms.OutputModule('PoolOutputModule',                                                                                                                  
                               fileName = cms.untracked.string('jettoolbox.root'),                                                                              
                               outputCommands = cms.untracked.vstring([
                                                                       # 'keep *_ak4PFJetsCHS_*_*',                                                                    
                                                                       # 'keep *_patJetsAK4PFCHS_*_*',                                                                  
                                                                       # 'keep *_ca8PFJetsCHS_*_*',                                                                     
                                                                       # 'keep *_patJetsCA8PFCHS_*_*',                                                                  
                                                                       # 'keep *_ak8PFJetsCHS_*_*',                                                                     
                                                                       # 'keep *_patJetsAK8PFCHS_*_*',                                                                  
                                                                      'keep *_slimmedJets_*_*',                                                                  
                                                                      'keep *_slimmedJetsAK8_*_*',                                                                  
                                                                       ])                                                                                           
                               )

#### NOT RUNNING OUTPUT MODULE ######                                                                                                                              
# process.endpath = cms.EndPath(process.out)    


# ----------------------- Jet Tool Box  -----------------
# ----- giulia test: do not recluster ak4 and ca8 jets to save time --------


# ##Load the toolBoxMiniHelper
# #from RecoJets.JetProducers.jettoolboxMiniHelper_cff import * 
# ##for some reason doesn't work  ?__?
# ## just copy&paste the cff here

process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.slimmedGenJetsAK8 = ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)

# process.ak4PFJets.src = 'packedPFCandidates'
# process.ak4PFJets.doAreaFastjet = True

# process.ak4PFJetsCHS = process.ak4PFJetsCHS.clone(src = 'chs', doAreaFastjet = True) #boh ?
# process.ak8PFJetsCHS = process.ak4PFJetsCHS.clone(src = 'chs', doAreaFastjet = True, rParam = 0.8)
# process.ca8PFJetsCHS = process.ca4PFJets.clone(src = 'chs', doAreaFastjet = True, rParam = 0.8)

# process.ak4GenJets.src = 'packedGenParticles'
# process.ak8GenJets = process.ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)
# process.ca8GenJets = process.ca4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)

# process.fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'

# process.ak4PFJetsCHSPruned = ak5PFJetsCHSPruned.clone(
#     src='chs'
#     rParam = 0.4, 
#     jetPtMin = 15.0
#     )
# process.ak4PFJetsCHSFiltered = ak5PFJetsCHSFiltered.clone(
#     src='chs'
#     rParam = 0.4, 
#     jetPtMin = 15.0 
#     )                                                                             
# process.ak4PFJetsCHSTrimmed = ak5PFJetsCHSTrimmed.clone(
#     src='chs'
#     rParam = 0.4,
#     jetPtMin = 15.0  
#     )

# process.ak8PFJetsCHSPruned.src = 'chs'
# process.ak8PFJetsCHSTrimmed.src = 'chs'
# process.ak8PFJetsCHSFiltered.src = 'chs'

# process.ca8PFJetsCHSPruned.src = 'chs'
# process.ca8PFJetsCHSTrimmed.src = 'chs'
# process.ca8PFJetsCHSFiltered.src = 'chs'

# process.cmsTopTagPFJetsCHS.src = 'chs'

# from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
# from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

# addJetCollection(
#     process,
#     labelName = 'AK4PFCHS',
#     jetSource = cms.InputTag('ak4PFJetsCHS'),
#     algo = 'ak4',
#     rParam = 0.4,
#     jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
#     ) 

# addJetCollection(
#     process,
#     labelName = 'CA8PFCHS',
#     jetSource = cms.InputTag('ca8PFJetsCHS'),
#     algo = 'ca8',
#     rParam = 0.8,
#     jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     )                                                                                                                                                                   

# addJetCollection(
#     process,
#     labelName = 'AK8PFCHS',
#     jetSource = cms.InputTag('ak8PFJetsCHS'),
#     algo = 'ak8',
#     rParam = 0.8,
#     jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     ) 

# """
# switchJetCollection(
#     process,
#     jetSource = cms.InputTag('ak4PFJets'),
#     algo = 'ak4',
#     rParam = 0.4,
#     jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-1'),
#     # btagDiscriminators = ['jetBProbabilityBJetTags',
#     #                       'jetProbabilityBJetTags',
#     #                       'trackCountingHighPurBJetTags',
#     #                       'trackCountingHighEffBJetTags',
#     #                       'simpleSecondaryVertexHighEffBJetTags',
#     #                       'simpleSecondaryVertexHighPurBJetTags',
#     #                       'combinedSecondaryVertexBJetTags'
#     #                       ],
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     )
# """

# process.patJetsAK4PFCHS.addJetCharge   = False
# process.patJetsAK4PFCHS.addBTagInfo    = True
# process.patJetsAK4PFCHS.getJetMCFlavour = False
# process.patJetsAK4PFCHS.addAssociatedTracks = False
# process.patJetPartonMatchAK4PFCHS.matched='prunedGenParticles'
# process.patJetCorrFactorsAK4PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

# process.patJetsCA8PFCHS.addJetCharge   = False
# process.patJetsCA8PFCHS.addBTagInfo    = False   #For some reason this has to be False
# process.patJetsCA8PFCHS.getJetMCFlavour = False
# process.patJetsCA8PFCHS.addAssociatedTracks = False
# process.patJetPartonMatchCA8PFCHS.matched='prunedGenParticles'
# process.patJetCorrFactorsCA8PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

# process.patJetsAK8PFCHS.addJetCharge   = False
# process.patJetsAK8PFCHS.addBTagInfo    = False    #For some reason this has to be False
# process.patJetsAK8PFCHS.getJetMCFlavour = False
# process.patJetsAK8PFCHS.addAssociatedTracks = False
# process.patJetPartonMatchAK8PFCHS.matched='prunedGenParticles'
# process.patJetCorrFactorsAK8PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

# process.patJets.addJetCharge   = False
# process.patJets.addBTagInfo    = True
# process.patJets.getJetMCFlavour = False
# process.patJets.addAssociatedTracks = False
# process.patJetPartonMatch.matched = 'prunedGenParticles'
# process.patJetCorrFactors.primaryVertices = 'offlineSlimmedPrimaryVertices'


# process.load('RecoBTag.Configuration.RecoBTag_cff')
# process.load('RecoJets.Configuration.RecoJetAssociations_cff')
# process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

# process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag('ak4PFJetsCHS')
# process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag('ak4PFJetsCHS')
# process.ak4JetTracksAssociatorAtVertexPF.tracks = cms.InputTag('unpackedTracksAndVertices')
# process.ak8JetTracksAssociatorAtVertexPF=process.ak4JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak8PFJetsCHS'),
#                                                                                         coneSize = 0.8)
# process.ca8JetTracksAssociatorAtVertexPF=process.ak4JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8PFJetsCHS'),
#                                                                                         coneSize = 0.8)

# process.impactParameterTagInfos.primaryVertex = cms.InputTag('unpackedTracksAndVertices')
# process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag('unpackedTracksAndVertices','secondary','')
# process.combinedSecondaryVertex.trackMultiplicityMin = 1 #silly sv, uses un filtered tracks.. i.e. any pt

# process.load('FWCore.MessageLogger.MessageLogger_cfi')
# #process.MessageLogger.cerr.FwkReport.reportEvery = 10
# process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
# #process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# #process.options.allowUnscheduled = cms.untracked.bool(True)


# #Load the toolbox

# process.load('CMSDIJET.DijetRootTreeMaker.jettoolbox_cff')


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Test with Edmund Code... not working
#Load the GenEvent and GenParticles cff files
#process.load('CMSDIJET.DijetRootTreeMaker.RootTupleMakerV2_GenEventInfo_cfi')
#process.load('CMSDIJET.DijetRootTreeMaker.RootTupleMakerV2_GenParticles_cfi')

#-------------------------------------------------------
# Gen Particles Pruner
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.prunedGenParticlesDijet = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
    "drop  *  ", # by default
    "keep ( status = 3 || (status>=21 && status<=29) )", # keep hard process particles
    )
)


#------------- Recluster Gen Jets to access the constituents -------
#already in toolbox, just add keep statements

#process.out.outputCommands.append("keep *_ak4GenJets_*_*")
#process.out.outputCommands.append("keep *_ak8GenJets_*_*")
#process.out.outputCommands.append("keep *_ca8GenJets_*_*")

process.out.outputCommands.append("keep *_slimmedGenJets_*_*")
process.out.outputCommands.append("keep *_slimmedGenJetsAK8_*_*")

##-------------------- Define the source  ----------------------------



process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:miniAOD_RSGravToJJ_kMpl01_M-8000.root')
    #fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14miniaod__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__MINIAODSIM__PU20bx25_POSTLS170_V5-v1__00000__6AACD832-3707-E411-A167-001E672489D5.root')
    #fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14drAODSIM__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__00000__0622C950-58E4-E311-A595-0025904B130A.root')
    #fileNames = cms.untracked.vstring('file:2CEB70D6-D918-E411-B814-003048F30422.root')    
    #fileNames = cms.untracked.vstring('file:QstarToJJ_M_4000_TuneCUETP8M1_13TeV_pythia8__MINIAODSIM__Asympt50ns_MCRUN2_74_V9A-v1__70000__AA35D1E7-FEFE-E411-B1C5-0025905B858A.root')    
    fileNames = cms.untracked.vstring(THISINPUTFILE)    
)

# #Keep statements for valueMaps (link Reco::Jets to associated quantities)
# #You don't have to keep them to deswizzle below

# #pocess.out.outputCommands += ['keep *_pileupJetIdEvaluator_*_*',
# #                       'keep *_QGTagger_*_*']

# process.out.outputCommands += ['keep *_NjettinessCA8_*_*',
#                                'keep *_NjettinessAK8_*_*',
                               
# #                               'keep *_QJetsAdderCA8_*_*',
#                                'keep *_ca8PFJetsCHSPrunedLinks_*_*',
#                                'keep *_ca8PFJetsCHSTrimmedLinks_*_*',
#                                'keep *_ca8PFJetsCHSFilteredLinks_*_*',
#                                'keep *_cmsTopTagPFJetsCHSLinksCA8_*_*',
# #                               'keep *_QJetsAdderAK8_*_*',
#                                'keep *_ak8PFJetsCHSPrunedLinks_*_*',
#                                'keep *_ak8PFJetsCHSTrimmedLinks_*_*',
#                                'keep *_ak8PFJetsCHSFilteredLinks_*_*',
#                                'keep *_cmsTopTagPFJetsCHSLinksAK8_*_*',
#                                ]


# #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# #Deswizzle valueMaps and attach to PAT::Jets as userFloats

# #process.patJetsAK4PFCHS.userData.userFloats.src += ['pileupJetIdEvaluator:fullDiscriminant','QGTagger:qgLikelihood']
# #process.patJetsAK4PFCHS.userData.userInts.src   += ['pileupJetIdEvaluator:cutbasedId','pileupJetIdEvaluator:fullId']

# #run these modules before pat in the sequence

# process.patJetsCA8PFCHS.userData.userFloats.src += ['NjettinessCA8:tau1',
#                                                     'NjettinessCA8:tau2',
#                                                     'NjettinessCA8:tau3',
# #                                                    'QJetsAdderCA8:QjetsVolatility',
#                                                     'ca8PFJetsCHSPrunedLinks',
#                                                     'ca8PFJetsCHSTrimmedLinks',
#                                                     'ca8PFJetsCHSFilteredLinks',
#                                                     'cmsTopTagPFJetsCHSLinksCA8'
#                                                     ]



# process.patJetsAK8PFCHS.userData.userFloats.src += ['NjettinessAK8:tau1',
#                                                     'NjettinessAK8:tau2',
#                                                     'NjettinessAK8:tau3',
#  #                                                   'QJetsAdderAK8:QjetsVolatility',
#                                                     'ak8PFJetsCHSPrunedLinks',
#                                                     'ak8PFJetsCHSTrimmedLinks',
#                                                     'ak8PFJetsCHSFilteredLinks',
#                                                     'cmsTopTagPFJetsCHSLinksAK8'
#                                                     ]





##-------------------- User analyzer  --------------------------------


process.dijets     = cms.EDAnalyzer('DijetTreeProducer',
  ## JETS/MET ########################################
  # jetsAK4             = cms.InputTag('patJetsAK4PFCHS'), 
  # jetsAK8         = cms.InputTag('patJetsAK8PFCHS'),     
  # jetsCA8         = cms.InputTag('patJetsCA8PFCHS'),
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  #ptMinCA8         = cms.double(10),
  #mjjMin           = cms.double(700),
  #dEtaMax          = cms.double(1.3),

  #forse non serve perche' gia' aggiunti al pat::jet -giulia-
  #Cutbasedid        = cms.InputTag("pileupJetIdEvaluator:cutbasedId"),
  #fullDiscriminant  = cms.InputTag("pileupJetIdEvaluator:fullDiscriminant")
  #fullId            = cms.InputTag("pileupJetIdEvaluator:fullId"),
  ## MC ########################################
  pu               = cms.untracked.InputTag('addPileupInfo'),
  ptHat            = cms.untracked.InputTag('generator'),
  genParticles     = cms.InputTag('prunedGenParticlesDijet'),
  # genJetsAK4             = cms.InputTag('ak4GenJets'), 
  # genJetsAK8         = cms.InputTag('ak8GenJets'),     
  # genJetsCA8         = cms.InputTag('ca8GenJets'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
  genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),     


  ## trigger ###################################
  #triggerAlias     = cms.vstring('Fat','PFHT650','PFNoPUHT650','HT750','HT550'),
  triggerAlias     = cms.vstring('PFHT900'),
  triggerSelection = cms.vstring(
     #'HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v*',
     #'HLT_PFHT650_v*', #giulia : commented because not found in new entuples
     ### giulia
     #'HLT_HT650_v*',
     ### end giulia
     #'HLT_PFNoPUHT650_v*',
     #'HLT_HT750_v*',  
     #'HLT_HT550_v*'
     #'HLT_PFHT900_v*'
     '*'
     ),
  triggerConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMask         = cms.bool(False),
    l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  ),
  ## JECs ######################################
  redoJECs  = cms.bool(True),
  L1corrAK4 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt'),
  L2corrAK4 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L2Relative_AK4PFchs.txt'),
  L3corrAK4 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt'),
  L1corrAK8 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L1FastJet_AK8PFchs.txt'),
  L2corrAK8 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L2Relative_AK8PFchs.txt'),
  L3corrAK8 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L3Absolute_AK8PFchs.txt')
)


# ------------------ path --------------------------
process.filter_step = cms.Path(
              process.EcalDeadCellTriggerPrimitiveFilter*
              process.hcalLaserEventFilter)
              
# Noise filters added first in path as recommended in Twiki
process.p = cms.Path(
                     #process.CSCTightHaloFilter* does not work
                     #process.eeBadScFilter* does not work
                     #process.goodVertices*process.trackingFailureFilter* does not work
                     process.HBHENoiseFilter*
                     
                     
                     #process.prunedGenParticlesDijet* #GENPAR REMOVED
                     process.chs * 

                     #process.slimmedGenJetsAK8 * #GENPAR REMOVED
                     
                     #process.ak4PFJetsCHS *
                     #process.ak4GenJets *
                     #process.patJetsAK4PFCHS *
                     ##process.patJetPartonMatchAK4PFCHS *
                     ##process.patJetCorrFactorsAK4PFCHS* 
                     #process.ak8PFJetsCHS *
                     #process.ak8GenJets 
                     #process.NjettinessAK8 *
                     ##process.QJetsAdderAK8 * #questo rallenta - da capire
                     # process.ak8PFJetsCHSPrunedLinks *
                     # process.ak8PFJetsCHSFilteredLinks *
                     # process.ak8PFJetsCHSTrimmedLinks *
                     # process.patJetsAK8PFCHS *
                     # process.ca8PFJetsCHS *
                     # process.ca8GenJets * 
                     # process.NjettinessCA8 *
                     # #process.QJetsAdderCA8 * #questo rallenta - da capire
                     # process.ca8PFJetsCHSPrunedLinks *
                     # process.ca8PFJetsCHSFilteredLinks *
                     # process.ca8PFJetsCHSTrimmedLinks *
                     # # process.patJetsCA8PFCHS *    

                     # #process.pileupJetId        ##recipe not working for now
                     # #process.puJetIdSequence    ##recipe not working for now
                     # #process.pileupJetIdCalculator *  ##recipe not working for now
                     # #process.pileupJetIdEvaluator     ##recipe not working for now
                     # #process.QGTagger *

                     process.dijets 
                     )
