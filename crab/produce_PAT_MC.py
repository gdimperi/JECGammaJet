import FWCore.ParameterSet.Config as cms

import sys
sys.path.insert(0, '.')

from produce_PAT_COMMON import *

process = createProcess(True, True, True, False, "START53_V27")

#process.source.fileNames =  cms.untracked.vstring('file:input_mc.root')
process.source.fileNames =  cms.untracked.vstring(
#just for testing interactively! Z' file!
'/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_3_2/test/test_file_MINIAOD_for_JEC2015.root',
)
process.out.fileName = 'patTuple_PF2PAT_MC.root'
