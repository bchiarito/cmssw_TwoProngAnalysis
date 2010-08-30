#!/bin/bash                                                                                                
SAMPLE=$1;

eval `scramv1 runtime -sh`;

root -b <<!
.x make_diphoton_plots.C("${SAMPLE}");
.q
!

rfcp histograms_${SAMPLE}.root  /castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/mc/${SAMPLE};

# now make webpge

WEBDIR=/afs/cern.ch/user/t/torimoto/www/public_html/physics/diphoton/mc/$1;
mkdir $WEBDIR;
mv *$1*.png $WEBDIR;

cat > ${WEBDIR}/index.html <<EOF

<HTML>

<HEAD><TITLE>Diphoton: ${SAMPLE} Background MC</TITLE></HEAD>

<h1><FONT color="Blue">Diphoton: ${SAMPLE} Background MC</FONT></h1>
 
<BODY link="Red">
<FONT color="Red"> 

<h2> Trigger Info </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/trigger_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/trigger_${SAMPLE}.png">
</A>

<h2> Photon1 Kinematics </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/kinematics_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/kinematics_photon1_${SAMPLE}.png">
</A>

<h2> Photon1 HCAL & ECAL Isolation </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/caloIso_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/caloIso_photon1_${SAMPLE}.png">
</A>

<h2> Photon1 Track Isolation</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/trackIso_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/trackIso_photon1_${SAMPLE}.png">
</A>

<h2> Photon1 Spike Rejection</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/spike_photon1_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/spike_photon1_${SAMPLE}.png">
</A>

<h2> Photon2 Kinematics </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/kinematics_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/kinematics_photon2_${SAMPLE}.png">
</A>

<h2> Photon2 HCAL & ECAL Isolation </h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/caloIso_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/caloIso_photon2_${SAMPLE}.png">
</A>

<h2> Photon2 Track Isolation</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/trackIso_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/trackIso_photon2_${SAMPLE}.png">
</A>

<h2> Photon2 Spike Rejection</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/spike_photon2_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/spike_photon2_${SAMPLE}.png">
</A>

<h2> DiPhoton plots</h2>
<A HREF=http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/diphoton_${SAMPLE}.png>
<img height="600" src="http://torimoto.web.cern.ch/torimoto/physics/diphoton/mc/${SAMPLE}/diphoton_${SAMPLE}.png">
</A>

</FONT>
</BODY>
</HTML>

EOF
