#! /usr/bin/env python

import os, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = [
    # Format
    #["dataset publish name", "dataset name", "Low pt hat", "High pt hat", "Number of processed events", "Cross section"]
    # Gamma + jets
    ["/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM", "G_HT-100to200", 100, 200, 407798, 1534 ],
    ["/GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM", "G_HT-200to400", 200, 400, 459368, 489.9],
    ["/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM", "G_HT-400to600", 400, 600, 393507, 62.05],
    ["/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM", "G_HT-600toinf", 600, 100000, 483900, 20.87],
    ]
    # QCD
email = "giulia.dimperio@cern.ch"
d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

for dataset in datasets:

  dataset_path = dataset[0]
  dataset_name = dataset[1]

  pt_hat_low = dataset[2]
  pt_hat_high = dataset[3]

  events = dataset[4]
  xsection = dataset[5]

  ui_working_dir = ("crab_MC_%s") % (dataset_name)
  output_file = "crab_MC_%s.cfg" % (dataset_name)
  output_dir = "JetMET/MC/%s/%s" % (d, dataset_name)

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tOutput directory: %s" % output_dir)
  print("")

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@outputdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@processed_events@#%d#g\" -e \"s#@cross_section@#%.10e#g\" -e \"s#@low_pt_hat@#%f#g\" -e \"s#@high_pt_hat@#%f#g\" crab_MC.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, output_dir, email, events, xsection, pt_hat_low, pt_hat_high, output_file))

  cmd = "crab -create -submit -cfg %s" % (output_file)
  if options.run:
    os.system(cmd)
