#!/bin/bash

SAMPLE=$1;
DATA=$2; #"data"; # data or mc ?
MERGED=$3; #0; # (0) summed or individual (1) ?
CATEGORY=$4;
LUMI=$5; #pb
DATAMCDIR=$6;
DATADIR=$7;
MCDIR=$8;
SIGNAL=0;
VERSION=BLAH;
if [ "X"${SAMPLE} == "X" ]; then {
    echo "Please specify sample";
    exit;
} fi

if [ "X"${DATA} == "X" ]; then {
    echo "Please specify data or mc";
    exit;
} fi

if [ "X"${MERGED} == "X" ]; then {
    MERGED=0;
} fi

if [ "X"${LUMI} == "X" ]; then {
    echo "Setting LUMI to default of 1.0";
    LUMI=1.0;
} fi

if [ "X"${CATEGORY} == "X" ]; then {
	echo "running on ALL categories..";
    CATEGORY="ALL";
} fi
echo "running on " $CATEGORY;     

if [ ${SIGNAL} -eq 0 ]; then {
    SIGNAL=0;
    ISSIGNAL="kFALSE";
} else {
    SIGNAL=1;    
    ISSIGNAL="kTRUE";
} fi

#if [ "X"${DATAMCDIR} == "X" ]
#then
#    DATAMCDIR="DataMC"
#fi

echo $SAMPLE $DATA $MERGED $LUMI;

#eval `scramv1 runtime -sh`;

if [ ${DATA} == "data" ]; then {
  ISDATA="kTRUE";
  TITLE="Diphoton: ${SAMPLE} DATA, ${VERSION}";
#  VERSION=Data_41X_V3_K1.0_M100_Pt50_scaleMC2DATA;
#  VERSION=BLAH;
#  VERSION=Data_41X_V3_K1.0_M120_Pt60_scaleMC2DATA_PUreweight;
  VERSION=Data_41X_V3_K1.0_M140_Pt70_scaleMC2DATA_PUreweight;
#  VERSION=Data_41X_V3_K1.0_M120_Pt50_scaleMC2DATA_PUreweight;
#  VERSION=Data_41X_V3_K1.0_M100_Pt30_scaleMC2DATA_PUreweight;
#  VERSION=Data_41X_V3_K1.0_M100_Pt30_scaleMC2DATA;
#  VERSION=Data_41X_V3_K1.0_M100_scaleMC2DATA;
} else {
  ISDATA="kFALSE";
  TITLE="Diphoton: ${SAMPLE} MC, ${VERSION}";
} fi

##
#### MERGED FILES
##
if [ ${MERGED} -ne 0 ] ; then
{
if [ ${SAMPLE} == "allMC" ]; then {
root -b <<!
.L make_diphoton_plots.C
 mergeAllMC(${LUMI},"${CATEGORY}"); 
 .q
!
SAMPLE=${SAMPLE}_${LUMI}pb;
} elif [  ${SAMPLE} == "DataMC" ]; then {
root -b <<!
.L make_diphoton_plots.C
 overlayDataMC("${DATADIR}","${MCDIR}","${CATEGORY}",kTRUE);
 .q
!
} elif [  ${SAMPLE} == "allData" ]; then {
root -b <<!
.L make_diphoton_plots.C
mergeAllData();
 .q
!
} else {
root -b <<!
.L make_diphoton_plots.C
 merge("${SAMPLE}","${CATEGORY}");
 .q
!
SAMPLE=${SAMPLE}_all;
} fi
}
else {
##
#### UN-MERGED FILES
##
root -b <<!
 .L make_diphoton_plots.C
 make_diphoton_plots("${SAMPLE}",kTRUE,${ISDATA},${LUMI},"${CATEGORY}",${ISSIGNAL});
 .q
!
} fi

if [ ${DATA} == "data" ]; then
if [[ "${SAMPLE}" != "DataMC" && "${SAMPLE}" != "allData" ]]; then
SAMPLE=data_${SAMPLE}
fi
fi

CASTORDIR=/castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/${DATA}/${VERSION}/${SAMPLE};
AFSWEBDIR=/Users/toyoko/Work/cms/physics/diphoton/2011/web/${DATA}/${VERSION}/${SAMPLE};
WEBDIR=/Users/toyoko/Work/cms/physics/diphoton/2011/web/${DATA}/${VERSION}/${SAMPLE}

if [[ ${SAMPLE} == "DataMC" && "X"${DATAMCDIR} != "X" ]]; then
    CASTORDIR=/castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/${DATA}/${VERSION}/DataMC_${DATAMCDIR};
    AFSWEBDIR=/Users/toyoko/Work/cms/physics/diphoton/2011/web/${DATA}/${VERSION}/DataMC_${DATAMCDIR};
    WEBDIR=/Users/toyoko/Work/cms/physics/diphoton/2011/web/${DATA}/${VERSION}/DataMC_${DATAMCDIR};
fi


if [ ${CATEGORY} == "EBEB" ]; then {
    CASTORDIR=${CASTORDIR}_EBEB;
    AFSWEBDIR=${AFSWEBDIR}_EBEB;
    WEBDIR=${WEBDIR}_EBEB;
} elif [ ${CATEGORY} == "EEEE" ]; then {
    CASTORDIR=${CASTORDIR}_EEEE;
    AFSWEBDIR=${AFSWEBDIR}_EEEE;
    WEBDIR=${WEBDIR}_EEEE;
} elif [ ${CATEGORY} == "EBEE" ]; then {
    CASTORDIR=${CASTORDIR}_EBEE;
    AFSWEBDIR=${AFSWEBDIR}_EBEE;
    WEBDIR=${WEBDIR}_EBEE;
} elif [ ${CATEGORY} == "ALL" ]; then {
    CASTORDIR=${CASTORDIR}_ALL;
    AFSWEBDIR=${AFSWEBDIR}_ALL;
    WEBDIR=${WEBDIR}_ALL;
} elif [ ${CATEGORY} == "NOTEE" ]; then {
    CASTORDIR=${CASTORDIR}_NOTEE;
    AFSWEBDIR=${AFSWEBDIR}_NOTEE;
    WEBDIR=${WEBDIR}_NOTEE;
} else {
	echo "Wrong CATEGORY !!";
	exit;
} fi

#rfmkdir ${CASTORDIR};
#rfcp histograms_${SAMPLE}_${CATEGORY}.root ${CASTORDIR};

#VERSION=TEST;

# now make webpge
mkdir $AFSWEBDIR;
mv *${SAMPLE}*.png $AFSWEBDIR;
mv *${SAMPLE}*.pdf $AFSWEBDIR;
mv *${SAMPLE}*.C $AFSWEBDIR;

cat > ${AFSWEBDIR}/index.html <<EOF

<HTML>

<HEAD><TITLE>${TITLE}</TITLE></HEAD>

<h1><FONT color="Blue">${TITLE}</FONT></h1>
 
<BODY link="Red">
<FONT color="Red"> 

<h2> Primary Vertex Info </h2>

<A HREF=${WEBDIR}/h_Nvtx_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Nvtx_${SAMPLE}.png">

<h2> Trigger Info </h2>

<A HREF=${WEBDIR}/h_TrigHLT_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_TrigHLT_${SAMPLE}.png">

<h2> Photon1 Kinematics </h2>

<A HREF=${WEBDIR}/h_Photon1_pt_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_pt_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_pt_log_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_pt_log_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_pt_zoom_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_pt_zoom_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_eta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_eta_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_phi_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_phi_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_occupancy_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_occupancy_${SAMPLE}.png">

<A HREF=${WEBDIR}/h_Photon1_r9_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_r9_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_sigmaIetaIeta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_sigmaIetaIeta_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_sigmaEtaEta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_sigmaEtaEta_${SAMPLE}.png">

<h2> Photon1 HCAL & ECAL Isolation </h2>

<A HREF=${WEBDIR}/h_Photon1_hadOverEm_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_hadOverEm_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_hcalIso04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_hcalIso04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon1_hcalIso03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_hcalIso03_${SAMPLE}.png">-->
<A HREF=${WEBDIR}/h_Photon1_ecalIso04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_ecalIso04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon1_ecalIso03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_ecalIso03_${SAMPLE}.png">-->


<h2> Photon1 Track Isolation</h2>

<!--
<A HREF=${WEBDIR}/h_Photon1_trkIsoSumPtHollow03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoSumPtHollow03_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_trkIsoSumPtSolid03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoSumPtSolid03_${SAMPLE}.png">
-->
<A HREF=${WEBDIR}/h_Photon1_trkIsoSumPtHollow04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoSumPtHollow04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon1_trkIsoSumPtSolid04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoSumPtSolid04_${SAMPLE}.png">-->

<!--
<A HREF=${WEBDIR}/h_Photon1_trkIsoNtrksHollow03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoNtrksHollow03_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_trkIsoNtrksSolid03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoNtrksSolid03_${SAMPLE}.png">
-->
<A HREF=${WEBDIR}/h_Photon1_trkIsoNtrksHollow04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoNtrksHollow04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon1_trkIsoNtrksSolid04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_trkIsoNtrksSolid04_${SAMPLE}.png">-->

</A>

<h2> Photon1 Spike Rejection</h2>

<A HREF=${WEBDIR}/h_Photon1_swisscross_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_swisscross_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_e2e9_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_e2e9_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_e2x2e4x4_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_e2x2e4x4_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_severityLevel_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_severityLevel_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_recHitFlag_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_recHitFlag_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_maxRecHitTime_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_maxRecHitTime_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_maxRecHitTime_wide_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_maxRecHitTime_wide_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon1_e2e9_v_maxRecHitTime_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon1_e2e9_v_maxRecHitTime_${SAMPLE}.png">


<h2> Photon2 Kinematics </h2>

<A HREF=${WEBDIR}/h_Photon2_pt_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_pt_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_pt_log_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_pt_log_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_pt_zoom_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_pt_zoom_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_eta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_eta_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_phi_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_phi_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_occupancy_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_occupancy_${SAMPLE}.png">

<A HREF=${WEBDIR}/h_Photon2_r9_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_r9_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_sigmaIetaIeta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_sigmaIetaIeta_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_sigmaEtaEta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_sigmaEtaEta_${SAMPLE}.png">

<h2> Photon2 HCAL & ECAL Isolation </h2>

<A HREF=${WEBDIR}/h_Photon2_hadOverEm_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_hadOverEm_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_hcalIso04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_hcalIso04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon2_hcalIso03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_hcalIso03_${SAMPLE}.png">-->
<A HREF=${WEBDIR}/h_Photon2_ecalIso04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_ecalIso04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon2_ecalIso03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_ecalIso03_${SAMPLE}.png">-->


<h2> Photon2 Track Isolation</h2>

<!--
<A HREF=${WEBDIR}/h_Photon2_trkIsoSumPtHollow03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoSumPtHollow03_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_trkIsoSumPtSolid03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoSumPtSolid03_${SAMPLE}.png">
-->
<A HREF=${WEBDIR}/h_Photon2_trkIsoSumPtHollow04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoSumPtHollow04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon2_trkIsoSumPtSolid04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoSumPtSolid04_${SAMPLE}.png">-->

<!--
<A HREF=${WEBDIR}/h_Photon2_trkIsoNtrksHollow03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoNtrksHollow03_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_trkIsoNtrksSolid03_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoNtrksSolid03_${SAMPLE}.png">
-->
<A HREF=${WEBDIR}/h_Photon2_trkIsoNtrksHollow04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoNtrksHollow04_${SAMPLE}.png">
<!--<A HREF=${WEBDIR}/h_Photon2_trkIsoNtrksSolid04_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_trkIsoNtrksSolid04_${SAMPLE}.png">-->

</A>

<h2> Photon2 Spike Rejection</h2>

<A HREF=${WEBDIR}/h_Photon2_swisscross_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_swisscross_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_e2e9_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_e2e9_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_e2x2e4x4_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_e2x2e4x4_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_severityLevel_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_severityLevel_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_recHitFlag_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_recHitFlag_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_maxRecHitTime_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_maxRecHitTime_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_maxRecHitTime_wide_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_maxRecHitTime_wide_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Photon2_e2e9_v_maxRecHitTime_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Photon2_e2e9_v_maxRecHitTime_${SAMPLE}.png">

<h2> DiPhoton plots</h2>

<A HREF=${WEBDIR}/h_Diphoton_Minv_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_Minv_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_Minv_log_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_Minv_log_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_Minv_high_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_Minv_high_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_Minv_low_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_Minv_low_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_qt_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_qt_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_deltaPhi_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_deltaPhi_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_deltaEta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_deltaEta_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_deltaR_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_deltaR_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_Diphoton_cosThetaStar_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_Diphoton_cosThetaStar_${SAMPLE}.png">

<h2> Cumulative Background plot </h2>

<A HREF=${WEBDIR}/h_CumulativeBackground_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_CumulativeBackground_${SAMPLE}.png">


<h2> Fake Rate plots </h2>
<A HREF=${WEBDIR}/h_FakeRate_minv_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_minv_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_minv_log_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_minv_log_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_minv_high_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_minv_high_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_qt_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_qt_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_deltaPhi_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_deltaPhi_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_deltaEta_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_deltaEta_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_deltaR_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_deltaR_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_cosThetaStar_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_cosThetaStar_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_pt1_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_pt1_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_pt2_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_pt2_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_pt1_log_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_pt1_log_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_pt2_log_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_pt2_log_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_eta1_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_eta1_${SAMPLE}.png">
<A HREF=${WEBDIR}/h_FakeRate_eta2_${SAMPLE}.png><img height="300" src="${WEBDIR}/h_FakeRate_eta2_${SAMPLE}.png">

</FONT>
</BODY>
</HTML>

EOF

echo ${WEBDIR};