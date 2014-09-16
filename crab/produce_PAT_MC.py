import FWCore.ParameterSet.Config as cms

import sys
sys.path.insert(0, '.')

from produce_PAT_COMMON import *

process = createProcess(True, True, True, False, "START53_V27")

#process.source.fileNames =  cms.untracked.vstring('file:input_mc.root')
process.source.fileNames =  cms.untracked.vstring(
#'/store/mc/Summer12_DR53X/G_Pt-1400to1800_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/043C2EB6-371A-E211-AA0F-003048F0E18C.root'
#'/store/mc/Summer12_DR53X/G_Pt-470to800_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/044A52B2-DF16-E211-8F33-0030487F172B.root',
#
'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/028B7F68-2518-E211-9FAF-002481E0D480.root',
#'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/02CD4824-4A18-E211-9E48-0030487FA349.root',
#'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/0434834B-2518-E211-B153-003048F0E1EA.root',
#'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/04A1A0BE-0E18-E211-B34E-0030487F932D.root',
#'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/06505842-2C18-E211-AA52-003048C693FC.root',
#'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/06C87A70-2C18-E211-9885-0030487D5DB5.root',
#'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/08EB3DE6-4618-E211-826C-003048D43982.root',
#'/store/mc/Summer12_DR53X/G_Pt-300to470_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/0AFA57A8-2C18-E211-8902-002590494C98.root'
#
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/00D437FF-0E1B-E211-B452-0030487F92AB.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/00FEFD78-FF1A-E211-8436-0030487F1651.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/0237891C-FF1A-E211-B2C2-0030487F1717.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/045C80D1-FE1A-E211-8205-0030487F1651.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/06134123-011B-E211-BB33-0030487F1651.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/062E914B-B91A-E211-A5C0-0030487E52A1.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/064DF8B8-011B-E211-A988-0030487F92AB.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/06C19486-051B-E211-A0EF-0030487F1A4F.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/06F591E4-021B-E211-B184-0030487F1717.root',
#'/store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/0A70C7CF-031B-E211-8FD8-0030487F1A4F.root'
)
process.out.fileName = 'patTuple_PF2PAT_MC.root'
