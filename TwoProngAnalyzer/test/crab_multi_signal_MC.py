import CRABClient
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
from multiprocessing import Process
import sys
import subprocess

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--asym", action="store_true", dest="asym", default=False, help="")
parser.add_option("--loose", action="store_true", dest="loose", default=False, help="")
parser.add_option("-f","--filter", action="store", dest="filt", default="none", help="")

parser.add_option("-t","--tag", action="store", dest="tag", default="", help="tag for crab directory and output directory")
parser.add_option("-i","--input", action="store", dest="input", default="", help="input .dat file of miniaod datasets, one on each line")
parser.add_option("-v","--version", action="store", dest="version", default="v2", help="v2 (default) or v3")
parser.add_option("-p", "--print", action="store_true", dest="printOnly", default=False, help="will not run, only print summary of tasks to be created")
parser.add_option("-c", "--command", action="store_true", dest="command", default=False, help="print command to test crab config, then exit")
parser.add_option("-e", "--exec", action="store_true", dest="execute", default=False, help="print command to test crab config, run the command, then exit")
(options, args) = parser.parse_args()
if options.tag == "":
  dir_tag = ""
else:
  dir_tag = "_" + options.tag
out_tag = options.tag

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_signal_MC' + dir_tag
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['includeMCInfoBranches=True', 'includeSignalMCBranches=True']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.outLFNDirBase = '/store/user/bchiari1/cms_area/twoprong/trees/signal/' + out_tag
config.Data.publication = False
config.Data.unitsPerJob = 100000
config.Data.totalUnits = -1
config.Data.splitting = 'EventAwareLumiBased'
config.Data.lumiMask = ''
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_Rutgers'

if options.loose:
  config.JobType.pyCfgParams.extend(['includeLooseTwoProngs=True'])

if options.asym:
  config.JobType.pyCfgParams.extend(['includeAsymTwoProngs=True'])

if options.filt == "photon":
  config.JobType.pyCfgParams.extend(['filterOnPhoton=True'])

if options.filt == "muon":
  config.JobType.pyCfgParams.extend(['filterOnLepton=True'])    

if options.version == "v2":
  globalTagOptions = ['globalTag=mc2016', 'miniAODv=2']
  testfile = "/cms/chiarito/datafiles/signal/miniaod/MINIAOD_v2.root"
if options.version == "v3":
  globalTagOptions = ['globalTag=mc2016', 'miniAODv=3']
  testfile = "/cms/chiarito/datafiles/signal/miniaod/MINIAOD_v3.root"

config.JobType.pyCfgParams.extend(globalTagOptions)

if options.command or options.execute:
  command = "cmsRun " + config.JobType.psetName + " "
  for param in config.JobType.pyCfgParams:
    command += param + " "
  command += "sample=file:"+testfile+" maxEvents=10"
  if not options.execute:
    print command
    sys.exit()
  else:
    print command, '\n'
    subprocess.call(command.split())
    sys.exit()

if options.input == "":
  print "Must supply input file"
  sys.exit()
datasets = []
with open(options.input) as fi:
  for line in fi:
    datasets.append(line.strip())

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    if options.printOnly:
      print "creating mutlitask directory:", config.General.workArea
    for dataset in datasets:
      rN = ((dataset.replace('/','_')).split('_TuneCUETP8M1')[0])
      config.General.requestName = rN[1:len(rN)]
      config.Data.inputDataset = dataset
      if options.printOnly:
        print "  requestName =", config.General.requestName
        print "  inputDataset =", config.Data.inputDataset
      else:
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
