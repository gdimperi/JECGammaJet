import FWCore.ParameterSet.Config as cms

import __builtin__

__builtin__.runOnMC = False
__builtin__.runCHS = True
__builtin__.processCaloJets = False
__builtin__.correctMETWithT1 = True

__builtin__.GlobalTag = cms.string("GR_R_53_V8::All")

import sys
sys.path.insert(0, '.')

from produce_PAT_COMMON import *
process.source.fileNames =  cms.untracked.vstring('file:input_data.root')
process.out.fileName = 'patTuple_PF2PAT.root'
