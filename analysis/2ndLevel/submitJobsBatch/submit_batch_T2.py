#! /usr/bin/env python

import os
import sys
import optparse
import re

usage = "usage: To be run from DijetRootTreeAnalyzer/ :  python scripts/submit_batch_T2.py -i directory_containing_lists_to_run -o output_directory"

parser = optparse.OptionParser("submitAllGJetsID.py")
parser.add_option('-q', '--queue',       action='store',     dest='queue',       
                  help='run in batch in queue specified as option (default -q cmslong)', 
                  default='cmsan',
                  metavar="QUEUE")

parser.add_option("-i", "--inputList", dest="inputlist",
                  help="list of samples to be analyzed",
                  )

parser.add_option("-o", "--output", dest="output",
                  help="the directory OUTDIR contains the output of the program",
                  metavar="OUTDIR")

#parser.add_option("-m", "--match", dest="match",
#                  help="run only the samples containing this string in the name",
#                  default="")
parser.add_option("-c",  type=str, dest="template_cmssw", default="",
                  help="template of cmssw file to execute")

parser.add_option('-I', '--interactive',      
                  action='store_true',
                  dest='interactive',      
                  help='run the jobs interactively, 2 jobs at a time',
                  default=False)


(opt, args) = parser.parse_args()
################################################



os.system("mkdir -p batch")
pwd = os.environ['PWD']



ins = open(opt.inputlist, "r") 
j=0
for line in  ins:
  fullpath = line.split("dcap://cmsrm-se01.roma1.infn.it/")[1]
  #print fullpath
  file = os.path.basename(fullpath)
  #print namefile
  filename = os.path.splitext(file)[0]
  print ("process %s" % filename)
  line = line.rstrip('\n')
  filename = filename.rstrip('\n')


  dict = {
      "THISROOTFILE":"\""+opt.output+"/rootfile_"+filename+".root\"",
      #"THISGLOBALTAG":"\"MCRUN2_74_V8\"",
      "THISGLOBALTAG":"\"GR_P_V56\"",
      "THISINPUTFILE":"\""+line+"\""
      }

  ##create cmssw configuration file
  cmssw_cfgfile = "cfg/"+filename+"_cmssw.py"
  with open(cmssw_cfgfile, "wt") as fout:
    with open(opt.template_cmssw, "rt") as fin:
      for line_ in fin:
	#fout.write(line_)
      	line_=line_.rstrip()
	for k,v in dict.items():
	  line_ = re.sub(k,v,line_)
	fout.write(line_+"\n")




  command = "cmsRun cfg/"+filename+"_cmssw.py"
  print "submit "+command
  print ""
  
  logfile = "logfile_"+filename+".log"
  outputname = "batch/submit_"+filename+".src"
  outputfile = open(outputname,'w')
  outputfile.write('#!/bin/bash\n')
  outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc491\n')
  outputfile.write('cd '+pwd+' \n')
  outputfile.write('eval `scramv1 runtime -sh`\n')
  outputfile.write(command+"\n")
  print outputname 
  if opt.interactive==False:
      os.system("bsub -q "+opt.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)
  else:
      print logfile
      if imc==0: os.system(command+" >&! "+logfile+"&")
      else: os.system(command+" >&! "+logfile)
               
          
