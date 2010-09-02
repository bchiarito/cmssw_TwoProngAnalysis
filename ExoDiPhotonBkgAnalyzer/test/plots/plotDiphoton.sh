#!/bin/bash

SAMPLE=$1;
DATA=$2; #"data"; # data or mc ?
MERGED=$3; #0; # (0) summed or individual (1) ?

if [ "X"${SAMPLE} == "X" ]
then
    echo "Please specify sample";
    exit;
fi

if [ "X"${DATA} == "X" ]
then
    echo "Please specify data or mc";
    exit;
fi

if [ "X"${MERGED} == "X" ]
then
    MERGED=0;
fi

echo $SAMPLE $DATA $MERGED;

eval `scramv1 runtime -sh`;

if [ ${DATA} == "data" ]; then
  ISDATA="kTRUE";
  TITLE="Diphoton: ${SAMPLE}"
else
  ISDATA="kFALSE";
  TITLE="Diphoton: ${SAMPLE} Background MC"
fi

if [ ${MERGED} -ne 0 ] ; then
{
if [ ${SAMPLE} == "allMC" ]; then
root -b <<!
.L make_diphoton_plots.C
 mergeAllMC();
 .q
!
elif [  ${SAMPLE} == "DataMC" ]; then
echo "data MC comparison";
root -b <<!
.L make_diphoton_plots.C
 overlayDataMC();
 .q
!
else
root -b <<!
.L make_diphoton_plots.C
 merge();
 .q
!
SAMPLE=${SAMPLE}_all;
fi
}
else
root -b <<!
 .L make_diphoton_plots.C
 make_diphoton_plots("${SAMPLE}",kTRUE,${ISDATA});
 .q
!
fi

rfmkdir /castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/${DATA}/${SAMPLE};
rfcp histograms_${SAMPLE}.root  /castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/${DATA}/${SAMPLE};

# now make webpge

WEBDIR=/afs/cern.ch/user/t/torimoto/www/public_html/physics/diphoton/${DATA}/${SAMPLE};
mkdir $WEBDIR;
mv *${SAMPLE}*.png $WEBDIR;

cat > ${WEBDIR}/index.html <<EOF

<HTML>

<HEAD><TITLE>${TITLE}</TITLE></HEAD>

<h1><FONT color="Blue">${TITLE}</FONT></h1>
 
<BODY link="Red">
<FONT color="Red"> 

<h2> Trigger Info </h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_TrigHLT_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_TrigHLT_${SAMPLE}.png">

<h2> Photon1 Kinematics </h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_pt_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_pt_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_eta_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_eta_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_phi_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_phi_${SAMPLE}.png">

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_r9_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_r9_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_sigmaIetaIeta_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_sigmaIetaIeta_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_sigmaEtaEta_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_sigmaEtaEta_${SAMPLE}.png">

<h2> Photon1 HCAL & ECAL Isolation </h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_hadOverEm_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_hadOverEm_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_hcalIso04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_hcalIso04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_hcalIso03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_hcalIso03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_ecalIso04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_ecalIso04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_ecalIso03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_ecalIso03_${SAMPLE}.png">


<h2> Photon1 Track Isolation</h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtHollow03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtHollow03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtSolid03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtSolid03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtHollow04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtHollow04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtSolid04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoSumPtSolid04_${SAMPLE}.png">

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksHollow03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksHollow03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksSolid03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksSolid03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksHollow04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksHollow04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksSolid04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_trkIsoNtrksSolid04_${SAMPLE}.png">

</A>

<h2> Photon1 Spike Rejection</h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_swisscross_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_swisscross_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_severityLevel_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_severityLevel_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_recHitFlag_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_recHitFlag_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_maxRecHitTime_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon1_maxRecHitTime_${SAMPLE}.png">


<h2> Photon2 Kinematics </h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_pt_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_pt_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_eta_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_eta_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_phi_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_phi_${SAMPLE}.png">

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_r9_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_r9_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_sigmaIetaIeta_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_sigmaIetaIeta_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_sigmaEtaEta_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_sigmaEtaEta_${SAMPLE}.png">

<h2> Photon2 HCAL & ECAL Isolation </h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_hadOverEm_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_hadOverEm_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_hcalIso04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_hcalIso04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_hcalIso03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_hcalIso03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_ecalIso04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_ecalIso04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_ecalIso03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_ecalIso03_${SAMPLE}.png">


<h2> Photon2 Track Isolation</h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtHollow03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtHollow03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtSolid03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtSolid03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtHollow04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtHollow04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtSolid04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoSumPtSolid04_${SAMPLE}.png">

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksHollow03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksHollow03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksSolid03_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksSolid03_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksHollow04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksHollow04_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksSolid04_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_trkIsoNtrksSolid04_${SAMPLE}.png">

</A>

<h2> Photon2 Spike Rejection</h2>

<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_swisscross_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_swisscross_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_severityLevel_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_severityLevel_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_recHitFlag_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_recHitFlag_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_maxRecHitTime_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Photon2_maxRecHitTime_${SAMPLE}.png">

<h2> DiPhoton plots</h2>


<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_Minv_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_Minv_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_qt_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_qt_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_deltaPhi_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_deltaPhi_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_deltaEta_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_deltaEta_${SAMPLE}.png">
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_deltaR_${SAMPLE}.png><img height="300" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/h_Diphoton_deltaR_${SAMPLE}.png">


</FONT>
</BODY>
</HTML>

EOF
