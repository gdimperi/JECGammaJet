import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('globalTag',
          'FT53_V21A_AN6' , #Latest analysis GT for data (for CMSSW>=538_patch3)
          VarParsing.multiplicity.singleton,
          VarParsing.varType.string,
          "The globaltag to be used")

options.parseArguments()
if len(options.globalTag) == 0:
  raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")

import sys
sys.path.insert(0, '.')

from produce_PAT_COMMON import *

process = createProcess(False, True, True, False, options.globalTag)
#process = createProcess(False, True, True, False, "FT_53_V21_AN5")

#process.source.fileNames =  cms.untracked.vstring('file:input_data.root')
process.source.fileNames =  cms.untracked.vstring('file:/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_5_3_14/src/JetMETCorrections/GammaJetFilter/crab/007A5A14-1069-E211-8BDD-0025905964B4.root')
process.source.skipEvents = cms.untracked.uint32(601)

#process.source.fileNames =  cms.untracked.vstring('/store/user/sbrochet/../../data/Run2012B/SinglePhoton/AOD/22Jan2013-v1/20000/FC9D17BE-0C72-E211-BCC8-003048678A80.root')
process.out.fileName = 'patTuple_PF2PAT.root'
#process.out.fileName = 'patTuple_PF2PAT_Photon_Run2012A-22Jan2013.root'
