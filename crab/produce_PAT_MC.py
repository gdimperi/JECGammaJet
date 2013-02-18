import FWCore.ParameterSet.Config as cms

import sys
sys.path.insert(0, '.')

from produce_PAT_COMMON import *

process = createProcess(True, True, True, False, "START53_V7G")

process.source.fileNames =  cms.untracked.vstring('file:input_mc.root')
process.out.fileName = 'patTuple_PF2PAT_MC.root'
