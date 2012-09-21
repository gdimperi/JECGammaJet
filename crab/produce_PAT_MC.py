import FWCore.ParameterSet.Config as cms
import __builtin__

__builtin__.runOnMC = True
__builtin__.runCHS = True
__builtin__.processCaloJets = False
__builtin__.correctMETWithT1 = True

__builtin__.GlobalTag = cms.string("START52_V11C::All")

from produce_PAT_COMMON import *
process.source.fileNames =  cms.untracked.vstring('file:input_mc.root')
process.out.fileName = 'patTuple_PF2PAT_MC.root'
