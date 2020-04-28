import sys
import subprocess

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-c", "--command", action="store_true", dest="command", default=False, help="")
parser.add_option("-e", "--exec", action="store_true", dest="execute", default=False, help="")
parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="")
parser.add_option("--fnal", action="store_true", dest="fnal", default=False, help="")
parser.add_option("--asym", action="store_true", dest="asym", default=False, help="")
parser.add_option("--loose", action="store_true", dest="loose", default=False, help="")
parser.add_option("-f","--filter", action="store", dest="filt", default="none", help="")
(options, args) = parser.parse_args()
if len(args) == 0:
  tag = ""
  path = ""
elif len(args) == 1:
  tag = "_"+args[0]
  path = args[0]+'/'
else:
  print "bad arguments", len(args)
  sys.exit()

def modify_config(config, testfile):
  global tag
  config.General.workArea = config.General.workArea + tag
  config.Data.outLFNDirBase = config.Data.outLFNDirBase + path
  
  if options.test:
    config.Data.totalUnits = 1

  if options.fnal:
    config.Site.storageSite = 'T3_US_FNALLPC'
    config.Data.outLFNDirBase = config.Data.outLFNDirBase.replace('/cms_area/','/')
    tag+='_FNAL'

  if options.loose:
    config.JobType.pyCfgParams.extend(['includeLooseTwoProngs=True'])

  if options.asym:
    config.JobType.pyCfgParams.extend(['includeAsymTwoProngs=True'])

  if options.filt == "photon":
    config.JobType.pyCfgParams.extend(['filterOnPhoton=True'])

  if options.filt == "muon":
    config.JobType.pyCfgParams.extend(['filterOnLepton=True'])    

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
