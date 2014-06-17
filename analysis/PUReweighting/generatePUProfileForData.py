#! /bin/env python
# -*- coding: UTF-8 -*-

import argparse, os, tempfile, shutil, sys
from subprocess import call, PIPE, STDOUT, Popen

# JSON files for Run2012 before Winter13 rereco
#jsonFiles = [
    #["Photon_Run2012A-13Jul2012.json", [190645, 193621]],
    #["Photon_Run2012A-recover-06Aug2012.json", [190782, 190949]],
    #["SinglePhoton_Run2012B-13Jul2012.json", [193834, 196531]],
    #["SinglePhoton_Run2012C-24Aug2012.json", [198049, 198522]],
    #["SinglePhoton_Run2012C-PromptReco.json", [198941, 203002]],
    #["SinglePhoton_Run2012C-EcalRecover_11Dec2012.json", [201191, 201191]],
    #["SinglePhoton_Run2012D-PromptReco.json", [203894, 208686]]
 #]

# JSON files for Winter13 rereco. 
jsonFiles = [
    #["Photon_Run2012A-22Jan2013.json", [190645, 193621]],
    #["SinglePhoton_Run2012B-22Jan2013.json", [193834, 196531]],
    #["SinglePhoton_Run2012C-22Jan2013.json", [198049, 203742]],
    #["SinglePhoton_Run2012D-22Jan2013.json", [205515, 208686]]

    #["lumiSummary_Photon_Run2012A-22Jan2013_24Apr13.json", [190645, 193621]],       
    #["lumiSummary_SinglePhoton_Run2012B-22Jan2013_24Apr13.json", [193834, 196531]], 
    #["lumiSummary_SinglePhoton_Run2012C-22Jan2013_24Apr13.json", [198049, 203742]], 
    #["lumiSummary_SinglePhoton_Run2012D-22Jan2013_16May13.json", [205515, 208686]],
    #["lumiSummary_SinglePhoton_Run2012D-22Jan2013_part2_07June13.json", [205515, 208686]]

    ["lumiSummary_Photon_Run2012A-22Jan2013.json", [190645, 193621]],
    ["lumiSummary_SinglePhoton_Run2012B-22Jan2013.json", [193834, 196531]],
    ["lumiSummary_SinglePhoton_Run2012C-22Jan2013.json", [198049, 203742]],
    ["lumiSummary_SinglePhoton_Run2012D-22Jan2013.json", [203777, 208686]]

 ]

# Trigger list for all 2012
#triggers = {
    #tuple([190645, 199608]) : ["HLT_Photon30_CaloIdVL_IsoL_*",
                        #"HLT_Photon50_CaloIdVL_IsoL_*",
                        #"HLT_Photon90_CaloIdVL_IsoL_*",
                        #"HLT_Photon135_*",
                        #"HLT_Photon150_*",
                        #"HLT_Photon160_*"],
    #tuple([199609, 208686]) : ["HLT_Photon30_CaloIdVL_*",
                        #"HLT_Photon50_CaloIdVL_IsoL_*",
                        #"HLT_Photon90_CaloIdVL_*",
                        #"HLT_Photon135_*",
                        #"HLT_Photon150_*",
                        #"HLT_Photon160_*"]
    #}

# Trigger list for 2012, for Run2012ABCD.
triggers = {
    tuple([190645, 199608]) : ["HLT_Photon30_CaloIdVL_IsoL_*",
                        "HLT_Photon50_CaloIdVL_IsoL_*",
                        "HLT_Photon90_CaloIdVL_IsoL_*",
                        "HLT_Photon135_*",
                        "HLT_Photon150_*"],
    tuple([199609, 208686]) : ["HLT_Photon30_CaloIdVL_*",
                        "HLT_Photon50_CaloIdVL_IsoL_*",
                        "HLT_Photon90_CaloIdVL_*",
                        "HLT_Photon135_*",
                        "HLT_Photon150_*"]
    }

def parallelize(cmd):

  #print("Launching parallel for %r" % cmd)

  #p_cmd = ["parallel", "-j", "20", "-u"]
  #p_cmd.extend(cmd)

  #p = Popen(p_cmd, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
  #out = p.communicate(input = "\n".join(validJSONFiles))[0]

  print("To be ran: ")
  for jsonFile in validJSONFiles:
      #print((" ".join(cmd)).replace("{}", jsonFile).replace("{.}", os.path.splitext(jsonFile)[0]))
      cmd_to_execute = (" ".join(cmd)).replace("{}", jsonFile).replace("{.}", os.path.splitext(jsonFile)[0])
      #cmd_to_execute = (" ".join(cmd)).replace("{}", jsonFile).replace("{.}", "/afs/cern.ch/work/g/gdimperi/GammaJet/CMSSW_5_3_14/src/JetMETCorrections/GammaJetFilter/analysis/PUReweighting/")
      print(cmd_to_execute)

      #os.system("bsub -q 8nh  source /afs/cern.ch/work/g/gdimperi/GammaJet/CMSSW_5_3_14/src/JetMETCorrections/GammaJetFilter/analysis/PUReweighting/submit.sh %s " % cmd_to_execute )
      #os.system( cmd_to_execute )

  #print(out)
  print
  print


parser = argparse.ArgumentParser(description='Generate PU profile for each prescaled triggers')
#parser.add_argument('crab_json', nargs=1)
parser.add_argument('lumi_json', nargs=1)
args = parser.parse_args()

lumi_json = args.lumi_json[0]

if not os.path.exists(lumi_json):
  print("Error: %s not found" % lumi_json)
  exit(1)

#tempFolder = tempfile.mkdtemp(dir="/scratch")
tempFolder = tempfile.mkdtemp(dir="/afs/cern.ch/work/g/gdimperi/GammaJet/CMSSW_5_3_14/src/JetMETCorrections/GammaJetFilter/tmp")
print(tempFolder)

for trigger_runs, triggerslist in triggers.items():
  # Get list of JSON for this run range
  validJSONFiles = [];
  for json in jsonFiles:
    json_runs = json[1]

    if (json_runs[0] >= trigger_runs[0] and json_runs[0] <= trigger_runs[1]) or (json_runs[1] >= trigger_runs[0] and json_runs[1] <= trigger_runs[1]):
      validJSONFiles.append(json[0])

  for trigger in triggerslist:
    print("Generate PU profile for '%s'" % trigger)

    tempFilename = "%s_%d_%d" % (trigger.rstrip("_*"), trigger_runs[0], trigger_runs[1])
    outputCSV = os.path.join(tempFolder, "{.}_%s.csv" % tempFilename)
    outputJSON = os.path.join(tempFolder, "{.}_%s_JSON.txt" % tempFilename)
    outputROOT = "pu_truth_data_photon_2012_true_{.}_%s_75bins.root" % tempFilename

    print("\tRunning lumiCalc2...")
    triggerBof = "\"%s\"" % (trigger)
    parallelize(["lumiCalc2.py", "--without-checkforupdate", "--begin", str(trigger_runs[0]), "--end", str(trigger_runs[1]), "lumibyls", "-i", "{}", "--hltpath", triggerBof, "-o", outputCSV])
    #parallelize(["echo", "lumibyls", "-i", "{}", "--hltpath", trigger, "-o", outputCSV])

    print("\tRunning pileupReCalc_HLTpaths.py...")

    parallelize(["pileupReCalc_HLTpaths.py", "-i", outputCSV, "--inputLumiJSON", lumi_json, "-o", outputJSON])

    print("\tRunning pileupCalc...")
   
    parallelize(["pileupCalc.py", "-i", "{}", "--inputLumiJSON", outputJSON, "--calcMode", "true", "--minBiasXsec", "69300", "--maxPileupBin", "75", "--numPileupBins", "75", "--pileupHistName", "pileup", outputROOT, "--verbose"])

    print

 # Merge everything
#toMerge = {}
#for trigger_runs, triggerslist in triggers.items():
  #for json in jsonFiles:
    #json_name = os.path.splitext(json[0])[0]

    #for trigger in triggerslist:
      #root_file = "pu_truth_data_photon_2012_true_%s_%s_%d_%d_75bins.root" % (json_name, trigger.rstrip("_*"), trigger_runs[0], trigger_runs[1])
      #if os.path.exists(root_file):
        #toMerge.setdefault(trigger.rstrip("_*"), []).append(root_file)

#for trigger, files in toMerge.items():
  #output_file = "pu_truth_data_photon_2012_true_%s_75bins.root" % trigger
  #input_files = " ".join(files)
  #os.system("hadd %s %s" % (output_file, input_files))

#shutil.rmtree(tempFolder)
