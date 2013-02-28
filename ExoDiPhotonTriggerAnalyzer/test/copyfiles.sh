#!/bin/bash

sourcedir=$1;
destfile=$2;
option=$3;
eoscommand="/afs/cern.ch/project/eos/installation/0.1.0-22d/bin/eos.select";
chainresult="hadd -f "$destfile" ";
for i in $($eoscommand ls $sourcedir);
do
#echo $i;
chainresult=$chainresult"root://eoscms/"$sourcedir/$i" ";
#echo $chainresult;
done
echo $chainresult;
output=$($chainresult);
