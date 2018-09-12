#!/bin/tcsh

# Thicknesses we are shooting
setenv thicknesses "eta1p6 eta2p0 eta2p5" # 300, 200, 120 
# Noise cut
setenv noisecuts "3 5 10"
 
# Starting the loop through all thicknesses
foreach thick ($thicknesses)
echo "===================================================================================="
echo "Thickness $thick"

foreach ncut  ($noisecuts)
echo "------------------------"
echo "Noise cut ${ncut}"

setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${thick}/${ncut}/jobs"

setenv runumberslist ` ls ${workpath} | grep .job `

foreach run  ($runumberslist)

#echo ${run}
chmod 755 ${workpath}/${run}

echo "Sending ${run}"
bsub -q 8nh -o /tmp/junk ${workpath}/${run}
#echo "bsub -q 2nd -o /tmp/junk ${workpath}/${run}"

end

end 

end




