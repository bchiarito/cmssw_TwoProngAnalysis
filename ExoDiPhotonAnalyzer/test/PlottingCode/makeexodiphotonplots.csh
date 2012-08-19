#!/bin/csh
       
   set SAMPLE = $argv[1]
   set SAMPLETYPE = $argv[2]                       
   set MERGED = $argv[3]                        
   set LUMI = $argv[4]                   
   set JSON = $argv[5]    
                     
set TreeDIR=/afs/cern.ch/work/j/jcarson/private/DiPhotonTrees;
set HistoDIR=/afs/cern.ch/work/j/jcarson/private/DiPhotonTrees/histograms/${SAMPLE};

mkdir $HistoDIR


if (${SAMPLETYPE} == "data") then
    root -b<<!
    .L  CreateHistogramFiles.C+
    CreateHistogramFiles("${SAMPLE}","${SAMPLETYPE}","${JSON}")  
    makeplots("${SAMPLE}","${LUMI}","${JSON}","${SAMPLETYPE}")
   
   .q
   !
      endif

if (${SAMPLETYPE} == "mc") then    
   if (${MERGED} == "0") then  
     root -b <<!
	.L CreateHistogramFiles.C+ 
   CreateHistogramFiles("${SAMPLE}","${SAMPLETYPE}","${JSON}")
   makeplots("${SAMPLE}","${LUMI}","${JSON}","${SAMPLETYPE}")
   .q
    !
  endif 
    if (${MERGED} == "1") then 
        root -b <<!
        .L CreateHistogramFiles.C+          
	MakeCombinedMCHistos();
        .q
          !
        endif
        endif 
    
if (${SAMPLETYPE} == "datamc") then
        root -b <<!
       .L CreateHistogramFiles.C+ 
    OverlayMCandData("${SAMPLE}","${LUMI}","${JSON}")
.q
!
endif

if (${SAMPLETYPE} == "datastitchmc") then
        root -b <<!
        .L CreateHistogramFiles.C+
 StitchBackgroundandMC("${SAMPLE}","${LUMI}","${JSON}")
.q
!
endif

if (${SAMPLETYPE} == "fake") then
        root -b <<!
        .L CreateHistogramFiles.C+
 fakeratehistos("${SAMPLE}","${LUMI}","${JSON}")
.q
!
endif


if (${SAMPLETYPE} == "datafakemc") then
        root -b <<!
        .L CreateHistogramFiles.C+
 OverlayMCFake("${SAMPLE}","${LUMI}","${JSON}")
.q
!
endif


EOF
EOF
set CASTORDIR=/castor/cern.ch/user/j/jcarson/DiPhotonTrees/Histograms/${SAMPLE};
echo $CASTORDIR;
nsmkdir $CASTORDIR
rfcp histograms_${SAMPLE}*.root /castor/cern.ch/user/j/jcarson/DiPhotonTrees/Histograms/${SAMPLE};
cp histograms_${SAMPLE}*.root /afs/cern.ch/work/j/jcarson/private/DiPhotonTrees/histograms
set WEBDIR=/afs/cern.ch/user/j/jcarson/www/physics/diphoton/${SAMPLETYPE}/${SAMPLE};


mkdir $WEBDIR
mv *${SAMPLE}*.png $WEBDIR;

   
EOF
