import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-c", "--command", action="store_true", dest="command", default=False, help="")
parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="")
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
