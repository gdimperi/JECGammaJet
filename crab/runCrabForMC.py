#! /bin/env python
# Launch crab for every datasets in mc_signal_datasets.list

import os, shutil
from optparse import OptionParser

isCastor = os.system("uname -n | grep cern &> /dev/null") == 0

parser = OptionParser()
parser.add_option("-p", "--path", dest="path", type="string", help="where to store crab folders")
parser.add_option("--mc", dest="mc", type="choice", choices=['Summer12', 'Fall11'], help="The MC type (Summer12 or Fall11)", default="Summer12")

(options, args) = parser.parse_args()

version = 1

if options.path is None or not os.path.isdir(options.path):
  parser.error("you must specify a valid path")

datasetfile = "%s/datasets.list" % options.path

os.path.exists(datasetfile) or exit("'%s' does not exist. Exiting." % datasetfile)

f = open(datasetfile)
datasets = f.readlines();
f.close();

print "Working ..."
i = 1;

# Copy common_dump_config.py and dump_MC.py
shutil.copy2("produce_PAT_MC.py", "%s/produce_PAT_MC.py" % options.path)

for dataset in datasets:

  print "Processing entry %02d/%02d" % (i, len(datasets))
  i = i + 1
  dataset = dataset.rstrip()
  name = dataset.replace("/", "_")[1:dataset.lower().find('tunez2') - 1]
  publish_name = "JetMet_PF2PAT_2012_27May_%s" % options.mc

  if (isCastor):
    template = "crab_MC.cfg.template.castor"
    remoteOutputDir = "/user/s/sbrochet/JetMet/MC/%s/%s" % (options.mc, name)
  else:
    template = "crab_MC.cfg.template.ipnl"
    remoteOutputDir = "JetMet/MC/%s/%s" % (options.mc, name)

  outputConfigFile = "%s/crab_MC_%s.cfg" % (options.path, name)

  os.system("sed -e \"s/@datasetname@/%s/g\" -e \"s/@uiworkingdir@/crab_%s_v%d/g\" -e \"s:@outputdir@:%s:g\" -e \"s:@publish_data_name@:%s:g\" %s > %s" % (dataset.replace("/", "\\/"), name, version, remoteOutputDir, publish_name, template, outputConfigFile))
  # Be sure to create the directory on dpm, and chmod it to 777
  if (isCastor):
    fullRemoveOutputDir = ("/dpm/in2p3.fr/home/cms/data/store/user/sbrochet/%s") % (remoteOutputDir)
  else:
    fullRemoveOutputDir = ("/castor/cern.ch/%s") % (remoteOutputDir)

  os.system("rfmkdir %s &> /dev/null" % fullRemoveOutputDir)
  
  os.system("cd %s && crab -create -submit -cfg crab_MC_%s.cfg && cd -" % (options.path, name))

os.remove("%s/produce_PAT_MC.py" % options.path)
os.remove("%s/produce_PAT_MC.pyc" % options.path)

print "Done."
