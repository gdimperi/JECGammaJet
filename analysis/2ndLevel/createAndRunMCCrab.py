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
    ["/G_Pt-30to50_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-30to50", 30, 50, 0, 0],
    ["/G_Pt-50to80_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-50to80", 50, 80, 0, 0],
    ["/G_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-80to120", 80, 120, 0, 0],
    ["/G_Pt-120to170_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-120to170", 120, 170, 0, 0],
    ["/G_Pt-170to300_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-170to300", 170, 300, 0, 0],
    ["/G_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-300to470", 300, 470, 0, 0],
    ["/G_Pt-470to800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-470to800", 470, 800, 0, 0],
    ["/G_Pt-800to1400_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-800to1400", 800, 1400, 0, 0],
    ["/G_Pt-1400to1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-1400to1800", 1400, 1800, 0, 0],
    ["/G_Pt-1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "G_Pt-1800", 1800, 8000, 0, 0],

    # QCD
    ["/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM", "QCD_Pt-15-3000", -1, -1, 0, 0],
    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

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

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@outputdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@processed_events@#%d#g\" -e \"s#@cross_section@#%f#g\" -e \"s#@low_pt_hat@#%f#g\" -e \"s#@high_pt_hat@#%f#g\" crab_MC.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, output_dir, email, events, xsection, pt_hat_low, pt_hat_high, output_file))

  cmd = "crab -create -submit -cfg %s" % (output_file)
  if options.run:
    os.system(cmd)
