#!/bin/bash

SAMPLE=$1;

mkdir $SAMPLE;
cp exodiphotonsignalmcanalyzer_cfg.py $SAMPLE;

cat > ${SAMPLE}/crab.cfg <<EOF


[CRAB]

jobtype = cmssw
scheduler = glite
use_server = 1
server_name=slc5cern

[CMSSW]

datasetpath=/${SAMPLE}_7TeV-pythia6/Spring10-START3X_V26-v1/GEN-SIM-RECO

pset=exodiphotonsignalmcanalyzer_cfg.py

total_number_of_events=-1
#total_number_of_events=1000
events_per_job = 20000
### number_of_jobs = 10

[USER]
return_data = 1

[GRID]

### for RS gravitons
### and for SM diphotons
#ce_white_list = T2_US_Caltech
ce_black_list = T2_US_Caltech, T2_US_Nebraska

#ce_black_list = T2_US_MIT, T3_GR_Ioannina

EOF

echo " setting up crab env"
# setup crab environment
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh;
eval `scramv1 runtime -sh`;
#source /afs/cern.ch/cms/ccs/wm/scripts/Crab/CRAB_2_6_6/crab.sh;
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh;


cd ${SAMPLE}

echo " launching crab jobs"
crab -create;
crab -submit;
crab -status;
cd -;

exit
