#! /usr/bin/env python

import os, datetime, pwd

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = [

    ["/Photon/sbrochet-Photon_Run2012A-13Jul2012_22Feb13-v1-0a006afc6a827f1ca7844a0cdf309f58/USER", "Photon_Run2012A-13Jul2012", "FT_53_V6C_AN4"],
    ["/Photon/sbrochet-Photon_Run2012A-recover-06Aug2012_22Feb13-v1-d40b2edbf331963e2dd3f62e946de959/USER", "Photon_Run2012A-recover-06Aug2012", "FT_53_V6C_AN4"],
    ["/SinglePhoton/sbrochet-SinglePhoton_Run2012B-13Jul2012_22Feb13-v1-0a006afc6a827f1ca7844a0cdf309f58/USER", "SinglePhoton_Run2012B-13Jul2012", "FT_53_V6C_AN4"],
    ["/SinglePhoton/sbrochet-SinglePhoton_Run2012C-24Aug2012_22Feb13-v1-5e618a4e6f4561b627e6e2f68b03adbe/USER", "SinglePhoton_Run2012C-24Aug2012", "FT53_V10A_AN4" ],
    ["/SinglePhoton/sbrochet-SinglePhoton_Run2012C-PromptReco_22Feb13-v1-8307c25f31c3b8496f9ebdaf71f16d1d/USER", "SinglePhoton_Run2012C-PromptReco", "GR_P_V42_AN4" ],
    ["/SinglePhoton/sbrochet-SinglePhoton_Run2012C-EcalRecover_11Dec2012_22Feb13-v1-8165992360f0663abf30a1416afd8911/USER", "SinglePhoton_Run2012C-EcalRecover_11Dec2012", "FT_P_V42C_AN4"],
    ["/SinglePhoton/sbrochet-SinglePhoton_Run2012D-PromptReco_22Feb13-v1-bd1bc616d5258999184a61398b70058f/USER", "SinglePhoton_Run2012D-PromptReco", "GR_P_V42_AN4"],

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
  dataset_globaltag = dataset[2]

  ui_working_dir = ("crab_data_%s") % (dataset_name)
  output_file = "crab_data_%s.cfg" % (dataset_name)
  output_dir = "JetMET/data/%s/%s" % (d, dataset_name)

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tGlobal tag: %s" % dataset_globaltag)
  print("\tOutput directory: %s" % output_dir)
  print("")

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@globaltag@#%s#g\" -e \"s#@remote_dir@#%s#g\" -e \"s#@dataset_name@#%s#g\" crab_data.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, email, dataset_globaltag, output_dir, dataset_name, output_file))

  cmd = "crab -create -submit -cfg %s" % (output_file)
  if options.run:
    os.system(cmd)
