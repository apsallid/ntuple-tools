#!/bin/tcsh

cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/CMSSW_10_3_0_pre4/src/
eval `scramv1 runtime -csh`
cd -



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

hadd output_partGun_PDGid22_x10_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${thick}/${ncut}/output/*.root 

end

end 





