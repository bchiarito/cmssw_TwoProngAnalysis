import classad
import htcondor
import sys
import os

from optparse import OptionParser
usage = "Usage: %prog [options] <directory with MINIAOD signal mc files>"
parser = OptionParser(usage=usage)
parser.add_option("-t","--tag", action="store", dest="tag", default="", help="tag for run, used for condor log file and optional .out files")
parser.add_option("--v2", action="store_true", dest="old", default=False, help="use when running on Miniaodv2 (default v3)")
parser.add_option("-o", "--out", action="store_true", dest="out", default=False, help="save .out files with job stdout and stderr")
parser.add_option("--asym", action="store_true", dest="asym", default=False, help="inactive")
parser.add_option("--loose", action="store_true", dest="loose", default=False, help="inactive")

(options, args) = parser.parse_args()

# setup
if options.old: old_data_flag = 'True'
else: old_data_flag = 'False'
directory = args[0]
job_dir = "condor_signal_multirun_"+options.tag
if job_dir in os.listdir('.'):
  print "Job directory " + job_dir + "already exists!"
  print "  please delete or use a different tag"
  sys.exit()
os.system("mkdir "+job_dir)

# schedd
coll = htcondor.Collector()
schedd_deamon_ad = coll.locate(htcondor.DaemonTypes.Schedd)
schedd = htcondor.Schedd(schedd_deamon_ad)

# submit job
with schedd.transaction() as txn:
  sub = htcondor.Submit()
  sub['executable'] = "condor_ntuplizer.sh"
  sub['log'] = job_dir+"/condor.log"
  for fi in os.listdir(directory):
    if os.path.isdir(fi): continue
    fi_path = os.path.join(directory, fi)
    fi_tag = fi.replace('.root','')
    if options.out:
      sub['output'] = job_dir + "/" + fi_tag+"_condor.out"
      sub['error'] = job_dir + "/" + fi_tag+"_condor.out"
    sub['arguments'] = fi_path + " " + old_data_flag + " " + fi_tag + " " + job_dir
    print sub['arguments']
    sub.queue(txn)

# query job
for job_ad in schedd.xquery(requirements='Owner=="chiarito"',projection=['Owner', 'ClusterId', 'ProcId', 'JobStatus']):
  print "  ", job_ad.__repr__()
