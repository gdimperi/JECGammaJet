#! /bin/env python
# -*- coding: UTF-8 -*-

import argparse, os, tempfile, shutil, sys
from subprocess import call, PIPE, STDOUT, Popen

jsonFiles = [
    "Photon_Run2012A-13Jul2012.json",
    "Photon_Run2012A-recover-06Aug2012.json",
    "SinglePhoton_Run2012B-13Jul2012.json",
    "SinglePhoton_Run2012C-24Aug2012.json",
    "SinglePhoton_Run2012C-EcalRecover_11Dec2012.json",
    "SinglePhoton_Run2012C-PromptReco.json",
    "SinglePhoton_Run2012D-PromptReco.json"
 ]

triggers = [
#    "HLT_Photon20_*"
    "HLT_Photon30_*",
    "HLT_Photon50_*",
#    "HLT_Photon75_*",
    "HLT_Photon90_*",
    "HLT_Photon135_*",
    "HLT_Photon150_*",
    "HLT_Photon160_*",
#    "HLT_Photon250_*",
#    "HLT_Photon300_*"
    ]

def parallelize(cmd):

  print("Launching parallel for %r" % cmd)

  p_cmd = ["parallel", "-j", "7", "-u"]
  p_cmd.extend(cmd)

  p = Popen(p_cmd, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
  out = p.communicate(input = "\n".join(jsonFiles))[0]

  print(out)


parser = argparse.ArgumentParser(description='Generate PU profile for each prescaled triggers')
#parser.add_argument('crab_json', nargs=1)
parser.add_argument('lumi_json', nargs=1)
args = parser.parse_args()

lumi_json = args.lumi_json[0]

if not os.path.exists(lumi_json):
  print("Error: %s not found" % lumi_json)
  exit(1)

tempFolder = tempfile.mkdtemp(dir="/scratch")
print(tempFolder)

for trigger in triggers:
  print("Generate PU profile for '%s'" % trigger)

  tempFilename = trigger.rstrip("_*")
  outputCSV = os.path.join(tempFolder, "{.}_%s.csv" % tempFilename)
  outputJSON = os.path.join(tempFolder, "{.}_%s_JSON.txt" % tempFilename)
  outputROOT = "pu_truth_data_photon_2012_true_{.}_%s_75bins.root" % tempFilename

  print("\tRunning lumiCalc2...")

  parallelize(["lumiCalc2.py", "--without-checkforupdate", "lumibyls", "-i", "{}", "--hltpath", trigger, "-o", outputCSV])
  #parallelize(["echo", "lumibyls", "-i", "{}", "--hltpath", trigger, "-o", outputCSV])

  print("\tRunning pileupReCalc_HLTpaths.py...")

  parallelize(["pileupReCalc_HLTpaths.py", "-i", outputCSV, "--inputLumiJSON", lumi_json, "-o", outputJSON])

  print("\tRunning pileupCalc...")
 
  parallelize(["pileupCalc.py", "-i", "{}", "--inputLumiJSON", outputJSON, "--calcMode", "true", "--minBiasXsec", "69300", "--maxPileupBin", "75", "--numPileupBins", "75", "--pileupHistName", "pileup", outputROOT, "--verbose"])

  print

#shutil.rmtree(tempFolder)
