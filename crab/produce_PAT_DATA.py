import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('globalTag',
          '',
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

process.source.fileNames =  cms.untracked.vstring('file:input_data.root')
process.out.fileName = 'patTuple_PF2PAT.root'
