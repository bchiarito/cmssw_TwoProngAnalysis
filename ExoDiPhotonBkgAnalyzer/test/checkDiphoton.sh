#!/bin/bash                                                                                                
SAMPLE=$1;

# setup crab environment
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh;
eval `scramv1 runtime -sh`;
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh;

cd ${SAMPLE};

njobs=`crab -status | grep "Total Jobs" | awk '{print $2}'`
#echo $njobs;

nsubmitted=`crab -status | grep -vi job | grep -c Submitted `;
nscheduled=`crab -status | grep -vi job | grep -c Scheduled `;
nrunning=`crab -status | grep -vi job | grep -c Running `;
ndone=`crab -status | grep -vi job | grep -ic Done `;

if [ ${ndone} == ${njobs} ]
then
    echo $SAMPLE " is done..." $ndone
    crab -get;
else
    echo $SAMPLE "NOT yet done..." "done:" $ndone "run:" $nrunning "scheduled:" $nscheduled "submitted:" $nsubmitted "...TOTAL:" $njobs
fi

#cd -;
