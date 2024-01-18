#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "Arguments: $1 $2 $3 $4"
echo ''

echo "I am here:"
pwd
ls -l
echo ''

echo 'cmsenv...'
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
echo ''

echo 'now do cmsRun...'
cmsRun cmssw_twoprongntuplizer_cfg.py globalTag=mc2016 sample=file:$1 includeSignalMCBranches=True oldData=$2 out=$3 maxEvents=-1
echo ''

echo 'done job, contents now:'
ls -l
mv TwoProngNtuplizer_$3.root $4/
