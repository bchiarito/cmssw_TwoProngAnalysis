#!/bin/bash                                                                                                

SAMPLE=$1;
DATA="data";

eval `scramv1 runtime -sh`;

if [ ${DATA} == "data" ]; then
  ISDATA="kTRUE";
  TITLE="Diphoton: ${SAMPLE}"
else
  ISDATA="kFALSE";
  TITLE="Diphoton: ${SAMPLE} Background MC"
fi

echo $ISDATA;

root -b <<!
 .x make_diphoton_plots.C("${SAMPLE}",kTRUE,${ISDATA});
 .q
!

rfcp histograms_${SAMPLE}.root  /castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/${DATA}/${SAMPLE};

# now make webpge

WEBDIR=/afs/cern.ch/user/t/torimoto/www/public_html/physics/diphoton/${DATA}/$1;
mkdir $WEBDIR;
mv *$1*.png $WEBDIR;

cat > ${WEBDIR}/index.html <<EOF

<HTML>

<HEAD><TITLE>${TITLE}</TITLE></HEAD>

<h1><FONT color="Blue">${TITLE}</FONT></h1>
 
<BODY link="Red">
<FONT color="Red"> 

<h2> Trigger Info </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/trigger_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/trigger_${SAMPLE}.png">
</A>

<h2> Photon1 Kinematics </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/kinematics_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/kinematics_photon1_${SAMPLE}.png">
</A>

<h2> Photon1 HCAL & ECAL Isolation </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/caloIso_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/caloIso_photon1_${SAMPLE}.png">
</A>

<h2> Photon1 Track Isolation</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/trackIso_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/trackIso_photon1_${SAMPLE}.png">
</A>

<h2> Photon1 Spike Rejection</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/spike_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/spike_photon1_${SAMPLE}.png">
</A>

<h2> Photon2 Kinematics </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/kinematics_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/kinematics_photon2_${SAMPLE}.png">
</A>

<h2> Photon2 HCAL & ECAL Isolation </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/caloIso_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/caloIso_photon2_${SAMPLE}.png">
</A>

<h2> Photon2 Track Isolation</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/trackIso_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/trackIso_photon2_${SAMPLE}.png">
</A>

<h2> Photon2 Spike Rejection</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/spike_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/spike_photon2_${SAMPLE}.png">
</A>

<h2> DiPhoton plots</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/diphoton_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/${DATA}/${SAMPLE}/diphoton_${SAMPLE}.png">
</A>

</FONT>
</BODY>
</HTML>

EOF
