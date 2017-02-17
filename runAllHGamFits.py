import ROOT
import os
import sys
import subprocess

inFile = os.getcwd()+"/samples/diphoton/data1516_13TeV.xsec_hists.root"

# To use
command = "python runIndividualHGamFit.py -b --hist {0}"
batchScript = "scripts/batchScript_Template.sh"
scriptLocation = "submit_scripts"

def GetKeyNames(self,dir=""):
  self.cd(dir)
  return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames

def batchSubmit(command,signame = "") :

  # Open batch script as fbatchin
  fbatchin = open(batchScript, 'r')
  fbatchindata = fbatchin.read()
  fbatchin.close()

  # open modified batch script (fbatchout) for writing
  batchtempname = '{0}/HGamFits_batchScript_{1}.sh'.format(scriptLocation,signame)
  fbatchout = open(batchtempname,'w')
  fbatchoutdata = fbatchindata.replace("ZZZ",command) # In batch script replace ZZZ for submit command
  fbatchout.write(fbatchoutdata)
  modcommand = 'chmod 744 {0}'.format(batchtempname)
  subprocess.call(modcommand, shell=True)
  fbatchout.close()
  submitcommand = "qsub {0}".format(batchtempname)
  print submitcommand
  subprocess.call(submitcommand, shell=True)


fileList = [f for f in os.listdir("results/diphoton/") if os.path.isfile(os.path.join("results/diphoton/", f))]
fileList.sort()

runningList = [l for l in subprocess.check_output("qstat -f",shell=True).split("\n") if "Job_Name" in l]
runningList.sort()

readFile = ROOT.TFile.Open(inFile)
keys = readFile.GetKeyNames()
for key in sorted(keys) :
  
  if "fit" in key :
    continue

  #print "Checking item",key
  done = False

  # Did this output file already get made?
  for file in fileList :
    if key in file :
      #print "Completed: removing."
      done=True
      continue

  # Is this file currently getting made?
  for item in runningList :
    if key in item :
      #print "Completed: removing."
      done=True
      continue

  if done : continue

  print "Submitting job for hist",key

  thiscommand = command.format(key)
  batchSubmit(thiscommand,key)



