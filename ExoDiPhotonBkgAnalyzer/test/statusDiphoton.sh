#!/bin/bash                                                                                                
SAMPLE=$1;

# setup crab environment
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh;
eval `scramv1 runtime -sh`;
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh;

cd ${SAMPLE};

crab -status;

cd -;

