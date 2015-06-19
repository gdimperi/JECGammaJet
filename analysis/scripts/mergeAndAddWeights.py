#!/usr/bin/env python

import os
import sys
import argparse
import math
from ROOT import *
from array import array

usage = "usage:  python mergeAndAddWeights.py --inputList list_of_files_to_merge.txt -o ./ --xsec 2022100000 [in pb]  --lumi_tot 1000. [in /pb]"

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="",
        help="input list of files to be merged",
	    )
parser.add_argument("-o", "--outputDir", type=str, dest="outputDir", default="./",
        help="output directory",
	    )
parser.add_argument("--xsec", type=float, dest="xsec", default=1,
        help="cross section",
	    )
parser.add_argument("--lumi_tot", type=float, dest="lumi_tot", default=1,
        help="cross section",
	    )

args = parser.parse_args()
print args 

inputList = args.inputList
outputDir = args.outputDir
xsec = args.xsec
lumi_tot = args.lumi_tot
###################
#read input file
ins = open(args.inputList,"r")
files = " "

for line in ins:
  file = line.strip()
  files += str(" "+file)
  pathT2 = file.split("dcap://cmsrm-se01.roma1.infn.it/")[1]
  fullname = os.path.split(pathT2)[1]
  name = fullname.split("_")
  
#print files  

filename_out = outputDir+"/PhotonJet_2ndLevel_"+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+name[5]+".root" 
os.system("hadd -f "+filename_out+"  "+files )

inputFile = TFile(filename_out,"UPDATE")
h_sumOfWeights = inputFile.Get("gammaJet/h_sumW")
sumOfWeights = h_sumOfWeights.Integral()

##### update total luminosity ######
lumi = inputFile.Get("gammaJet/total_luminosity")
lumi.SetVal(lumi_tot) # in /pb

##### update tree with weight for total normalization #####
analysis_tree = inputFile.Get("gammaJet/analysis")
print analysis_tree

evtWeightTot = array("f", [0.] )
b_evtWeightTot = analysis_tree.Branch("evtWeightTot", evtWeightTot,"evtWeightTot/F")

for event in analysis_tree:
  evtWeightTot[0] = (xsec / sumOfWeights)
  b_evtWeightTot.Fill()

inputFile.cd("gammaJet")
analysis_tree.Write()
lumi.Write()
