#! /usr/bin/env python

import os, copy, datetime

datasets = {
    "/Photon/sbrochet-JetMet_PF2PAT_Run2012A_13July2012_ReReco_27October2012-edebda1c789bfeceb29ea0e5a54c0be8/USER": "Photon_Run2012A",
    "/SinglePhoton/sbrochet-JetMet_PF2PAT_Run2012B_13July2012_ReReco_14August2012-6eb44b0773a87f5608d57e6dc7821956/USER": "Photon_Run2012B",
    "/SinglePhoton/sbrochet-JetMet_PF2PAT_Run2012C_24Aug2012ReReco_29October2012-edebda1c789bfeceb29ea0e5a54c0be8/USER": "Photon_Run2012C_v1",
    "/SinglePhoton/sbrochet-JetMet_PF2PAT_Run2012C_PromptReco_v2_14August2012-6eb44b0773a87f5608d57e6dc7821956/USER": "Photon_Run2012C_v2",
    "/SinglePhoton/sbrochet-JetMet_PF2PAT_Run2012C_PromptReco_v2_14August2012-edebda1c789bfeceb29ea0e5a54c0be8/USER": "Photon_Run2012C_v2_part2"
    }

d = datetime.datetime.now().strftime("%d%b")

for dataset, ui in datasets.items():
  name = ("%s_%s") % (ui, d)
  ui_working_dir = ("crab_%s") % (name)
  output_file = "crab_%s.cfg" % (name)
  #output_dir_semie = ("MTT/Extracted/MC/Summer12/semie/%s_v%d" % (ui, version)).replace("/", "\\/")
  #output_dir_semimu = ("MTT/Extracted/MC/Summer12/semimu/%s_v%d" % (ui, version)).replace("/", "\\/")
  output_dir = ("JetMet/data/5.3/%s/%s" % (d, ui))

  print "Creating config file for '%s'" % (dataset)
  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@outputdir@#%s#g\" -e\"s:@datasetfolder@:%s:g\" crab_data.cfg.template.ipnl > %s" % (dataset, ui_working_dir, output_dir, ui, output_file))

  cmd = "crab -create -submit -cfg %s" % (output_file)
  os.system(cmd);
