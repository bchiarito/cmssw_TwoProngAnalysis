#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "Arguments: $1 $2"
echo ''

cd ${_CONDOR_SCRATCH_DIR}
echo "I am here:"
pwd
echo ''

echo "This is the scratch dir:"
echo ${_CONDOR_SCRATCH_DIR}
echo ''

source /cvmfs/cms.cern.ch/cmsset_default.sh
echo 'making CMSSW project directory'
eval "scramv1 project CMSSW CMSSW_9_4_0"
echo ''

cd CMSSW_9_4_0/src
echo 'cloning git repository'
git clone https://github.com/bchiarito/cmssw_TwoProngAnalysis.git TwoProngAnalysis
echo ''

echo 'cmsenv...'
eval `scramv1 runtime -sh`
echo 'scram b...'
scram b -j 4
echo ''

cd TwoProngAnalysis/TwoProngAnalyzer/test/
echo "moved to test dir, here is location and contents:"
pwd
ls -l
echo ''

echo 'now do cmsRun...'
cmsRun cmssw_twoprongntuplizer_cfg.py globalTag=mc2016 sample=file:$1 includeSignalMCBranches=True oldData=$2 
echo ''

echo 'done job, contents now:'
ls -l
