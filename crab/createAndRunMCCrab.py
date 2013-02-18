#! /usr/bin/env python

import os, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = {
    # Gamma + jets
    "/G_Pt-30to50_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-30to50",
    "/G_Pt-50to80_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-50to80",
    "/G_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-80to120",
    "/G_Pt-120to170_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-120to170",
    "/G_Pt-170to300_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-170to300",
    "/G_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-300to470",
    "/G_Pt-470to800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-470to800",
    "/G_Pt-800to1400_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-800to1400",
    "/G_Pt-1400to1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-1400to1800",
    "/G_Pt-1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "G_Pt-1800",

    # QCD
    "/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM": "QCD_Pt-15-3000",
    }

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

for dataset_path, dataset_name in datasets.items():

  dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)

  publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
  ui_working_dir = ("crab_MC_%s") % (dataset_name)
  output_file = "crab_MC_%s.cfg" % (dataset_name)

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tPublishing name: %s" % publish_name)
  print("")

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@publish_data_name@#%s#g\" -e \"s#@email@#%s#g\" crab_MC.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, publish_name, email, output_file))

  cmd = "crab -create -submit -cfg %s" % (output_file)
  if options.run:
    os.system(cmd)
