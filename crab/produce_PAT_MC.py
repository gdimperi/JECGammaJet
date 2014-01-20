import FWCore.ParameterSet.Config as cms

import sys
sys.path.insert(0, '.')


from produce_PAT_COMMON import *
#from produce_PAT_COMMON_test import *

#process = createProcess(True, True, True, False, True, "START53_V21")
process = createProcess(True, False, True, False, "START53_V27")

process.source.fileNames =  cms.untracked.vstring('file:input_mc.root')
process.source.fileNames =  cms.untracked.vstring('/store/relval/CMSSW_5_3_6/RelValZEE/GEN-SIM-RECO/PU_START53_V14-v1/0003/2C92DB85-E82C-E211-B8DE-003048D37560.root')
#process.source.fileNames =  cms.untracked.vstring('file:/afs/cern.ch/work/g/gdimperi/CMSSW_5_3_9_patch2/src/JetMETCorrections/GammaJetFilter/crab/G_Pt-50to80_TuneZ2star_8TeV_pythia6_AODSIM_PU_S10_START53_V7A-v1.root')
#process.source.fileNames =  cms.untracked.vstring('/store/group/phys_jetmet/gdimperi/dataset_gjet/G_Pt-30to50_TuneZ2star_8TeV_pythia6_AODSIM_PU_S10_START53_V7A-v1.root')
#process.out.fileName = 'patTuple_PF2PAT_MC.root'
#process.out.fileName = 'patTuple_PF2PAT_MC_G_Pt-50to80__step0-1.root'
process.out.fileName = 'patTuple_PF2PAT_MC_G_Pt-50to80__step1_test.root'

