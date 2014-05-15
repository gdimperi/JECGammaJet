#! /usr/bin/env python

import os, datetime, pwd

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

#global_json = "Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"
#global_json = "crab_data_Photon_Run2012A-22Jan2013/res//missingLumiSummary.json"
#global_json = "crab_data_SinglePhoton_Run2012B-22Jan2013/res//missingLumiSummary.json"
#global_json = "crab_data_SinglePhoton_Run2012C-22Jan2013/res//missingLumiSummary.json"
#global_json = "crab_data_SinglePhoton_Run2012D-22Jan2013/res//missingLumiSummary.json"
global_json = "crab_data_SinglePhoton_Run2012D-22Jan2013_splitted//res//missingLumiSummary_split3.json"

datasets = [

#    ["/Photon/Run2012A-22Jan2013-v1/AOD", "Photon_Run2012A-22Jan2013", "", "FT_53_V21_AN3", [190456, 193621]],                  
#    ["/SinglePhoton/Run2012B-22Jan2013-v1/AOD", "SinglePhoton_Run2012B-22Jan2013", "", "FT_53_V21_AN3", [193833, 196531]],      
#    ["/SinglePhoton/Run2012C-22Jan2013-v1/AOD", "SinglePhoton_Run2012C-22Jan2013", "", "FT_53_V21_AN3", [198022, 203742]],      
#    ["/SinglePhotonParked/Run2012D-22Jan2013-v1/AOD", "SinglePhoton_Run2012D-22Jan2013", "", "FT_53_V21_AN3", [203768, 208686]],

#  ["/Photon/Run2012A-22Jan2013-v1/AOD", "Photon_Run2012A-22Jan2013", "", "FT_53_V21_AN5", [190456, 193621]],                  
#  ["/SinglePhoton/Run2012B-22Jan2013-v1/AOD", "SinglePhoton_Run2012B-22Jan2013", "", "FT_53_V21_AN5", [193833, 196531]],      
#  ["/SinglePhoton/Run2012C-22Jan2013-v1/AOD", "SinglePhoton_Run2012C-22Jan2013", "", "FT_53_V21_AN5", [198022, 203742]],      
  ["/SinglePhotonParked/Run2012D-22Jan2013-v1/AOD", "SinglePhoton_Run2012D-22Jan2013", "", "FT_53_V21_AN5", [203768, 208686]],e

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
  dataset_json = dataset[2]
  if len(dataset_json) == 0:
    dataset_json = global_json

  dataset_globaltag = dataset[3]

  runselection = ""
  if len(dataset) > 4 and len(dataset[4]) == 2:
    runselection = "runselection = %d-%d" % (dataset[4][0], dataset[4][1])

  date = "07Feb14"
  publish_name = "%s_%s-v%d" % (dataset_name, date, version)
  ui_working_dir = ("crab_data_%s_splittedNew_3") % (dataset_name)
  #ui_working_dir2 = ("crab_data_%s_splitted2_2") % (dataset_name)
  output_file1 = "crab_data_%s_splittedNew_3.cfg" % (dataset_name)
  #output_file2 = "crab_data_%s_splitted2_2.cfg" % (dataset_name)
  output_dir = "/store/group/phys_jetmet/gdimperi/step1/data/%s/%s" % (d, dataset_name)

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tJSON: %s" % dataset_json)
  print("\tGlobal tag: %s" % dataset_globaltag)
  print("\tRun selection: %s" % runselection)
  print("\tPublishing name: %s" % publish_name)
  print("")

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@lumi_mask@#%s#g\" -e \"s#@runselection@#%s#g\" -e \"s#@publish_data_name@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@remote_dir@#%s#g\" -e \"s#@globaltag@#%s#g\" crab_data.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, dataset_json, runselection, publish_name, email, output_dir, dataset_globaltag, output_file1))

  #os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@lumi_mask@#%s#g\" -e \"s#@runselection@#%s#g\" -e \"s#@publish_data_name@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@remote_dir@#%s#g\" -e \"s#@globaltag@#%s#g\" crab_data.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir2, dataset_json, runselection, publish_name, email, output_dir, dataset_globaltag, output_file2))
  

  if options.run:
    cmd = "crab -create 1-5000 -cfg %s" % (output_file1)
    os.system(cmd)

    #cmd = "crab -create 2500-4499 -cfg %s" % (output_file2)
    os.system(cmd)
    
    cmd = "crab -submit 1-499 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)
    
    cmd = "crab -submit 500-999 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)
    
    cmd = "crab -submit 1000-1499 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)

    cmd = "crab -submit 1500-1999 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)

    cmd = "crab -submit 2000-2499 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)

    cmd = "crab -submit 2500-2999 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)

    cmd = "crab -submit 3000-3499 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)

    cmd = "crab -submit 3500-3999 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)

    cmd = "crab -submit 4000-4499 -c %s" % (ui_working_dir)
    print("%s" % cmd)
    os.system(cmd)
