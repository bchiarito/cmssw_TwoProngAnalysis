import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-c", "--command", action="store_true", dest="command", default=False, help="")
parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="")
parser.add_option("--fnal", action="store_true", dest="fnal", default=False, help="")
parser.add_option("--asym", action="store_true", dest="asym", default=False, help="")
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

def print_command(pset, params, filename):
  command = "cmsRun " + pset + " "
  for param in params:
    command += param + " "
  command += "sample=file:"+filename+" maxEvents=10"
  print command

def modify_config(config):
  if options.test:
    config.Data.totalUnits = 1

  if options.fnal:
    config.Site.storageSite = 'T3_US_FNALLPC'
    config.Data.outLFNDirBase = config.Data.outLFNDirBase.replace('/cms_area/','/')
    global tag
    tag+='_FNAL'
  config.General.workArea = config.General.workArea + tag
  config.Data.outLFNDirBase = config.Data.outLFNDirBase + path

  if options.asym:
    config.JobType.pyCfgParams.extend(['flipAsymReq=True'])

  if options.filt == "photon":
    config.JobType.pyCfgParams.extend(['filterOnPhoton=True'])
    
