#! /usr/bin/env python

import os, datetime, pwd

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = [

    ["/Photon/sbrochet-Photon_Run2012A-22Jan2013_24Apr13-v1-37e3bf2409397e623ffd52beab84a202/USER", "Photon_Run2012A-22Jan2013", "FT_53_V21_AN4"],
    ["/SinglePhoton/sbrochet-SinglePhoton_Run2012B-22Jan2013_24Apr13-v1-37e3bf2409397e623ffd52beab84a202/USER", "SinglePhoton_Run2012B-22Jan2013", "FT_53_V21_AN4"],
    ["/SinglePhoton/sbrochet-SinglePhoton_Run2012C-22Jan2013_24Apr13-v1-37e3bf2409397e623ffd52beab84a202/USER", "SinglePhoton_Run2012C-22Jan2013", "FT_53_V21_AN4" ],
    ["/SinglePhotonParked/sbrochet-SinglePhoton_Run2012D-22Jan2013_16May13-v1-d8690a0448f852b4d70656ff31f27990/USER", "SinglePhoton_Run2012D-22Jan2013", "FT_53_V21_AN4" ],
    ["/SinglePhotonParked/sbrochet-SinglePhoton_Run2012D-22Jan2013_part2_07June13-v1-d8690a0448f852b4d70656ff31f27990/USER", "SinglePhoton_Run2012D-22Jan2013_part2", "FT_53_V21_AN4"],

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
