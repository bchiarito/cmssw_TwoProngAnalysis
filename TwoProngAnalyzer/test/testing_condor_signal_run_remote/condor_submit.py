import classad
import htcondor
import sys
import os

# schedd
coll = htcondor.Collector()
schedd_deamon_ad = coll.locate(htcondor.DaemonTypes.Schedd)
schedd = htcondor.Schedd(schedd_deamon_ad)

#print "all schedds:"
#schedd_deamons = coll.locateAll(htcondor.DaemonTypes.Schedd)
#for s in schedd_deamons:
#  print ' ', s["Name"], s['MyAddress']
#print "default schedd:"
#print ' ', schedd_deamon_ad["Name"], schedd_deamon_ad["MyAddress"]

# submit job
old_data_flag = 'True'
with schedd.transaction() as txn:
  sub = htcondor.Submit()
  sub['executable'] = "ntuplizer.sh"
  sub['log'] = "ntuplizer.log"
  directory = sys.argv[1]
  for fi in os.listdir(directory):
    if os.path.isdir(fi): continue
    fi_path = os.path.join(directory, fi)
    fi_tag = fi.replace('.root','')
    sub['output'] = fi_tag+"ntuplizer.out"
    sub['error'] = fi_tag+"ntuplizer.out"
    sub['transfer_output_files'] = "CMSSW_9_4_0/src/TwoProngAnalysis/TwoProngAnalyzer/test/TwoProngNtuplizer.root"
    sub['transfer_output_remaps'] = '"TwoProngNtuplizer.root = TwoProngNtuplizer_"+fi_tag+".root"'
    sub['arguments'] = fi_path + " " + old_data_flag
    print sub
    print sub.queue(txn)

# query job
for job_ad in schedd.xquery(requirements='Owner=="chiarito"',projection=['Owner', 'ClusterId', 'ProcId', 'JobStatus']):
  print "  ", job_ad.__repr__()
