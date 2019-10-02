#!/bin/tcsh

cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/hitcalibrationV11/fcpermip/CMSSW_11_0_X_2019-09-23-2300/src
eval `scramv1 runtime -csh`
cd -

# Thicknesses we are shooting
#setenv thicknesses "eta1p6 eta2p0 eta2p5" # 300, 200, 120 
#setenv thicknesses "CE_E_Front_120um CE_E_Front_200um CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um"
#setenv thicknesses "CE_H_Coarse_Scint"
setenv thicknesses "CE_H_Fine_Scint"
#setenv thicknesses "CE_E_Front_120um CE_E_Front_200um"
#setenv thicknesses "eta1p6 eta2p0" # 300, 200
#setenv thicknesses "eta2p0"
#setenv thicknesses "eta2p5" # 120
# Noise cut
setenv noisecuts "0 3 5 10"
#setenv noisecuts "0" 

#Options are zerostage, firststage, secondstage
setenv calibstage "zerostage"

# Starting the loop through all thicknesses
foreach thick ($thicknesses)
echo "===================================================================================="
echo "Thickness $thick"

foreach ncut  ($noisecuts)
echo "------------------------"
echo "Noise cut ${ncut}"

hadd -f ${calibstage}/output_partGun_PDGid22_x100_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/${thick}/${ncut}/output/*.root 

end

end 





