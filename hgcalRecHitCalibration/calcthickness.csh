#!/bin/tcsh

cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/CMSSW_10_3_0_pre4/src
eval `scramv1 runtime -csh`
cd -

#====>>> IMPORTANT: To avoid accidental overwrite you should put by hand 
#the latest version of hgcalRecHitCalibration.py in the ntuples directory. 

#Run mode True and then False if you want to check the initial histo 
setenv MODE "True"

# Thicknesses we are shooting
setenv thicknesses "eta1p6 eta2p0 eta2p5" # 300, 200, 120 
# Noise cut
setenv noisecuts "3 5 10"

# Starting the loop through all thicknesses
foreach thick ($thicknesses)
echo "===================================================================================="
echo "Thickness $thick"

if ( ${thick} == "eta1p6"  ) then 
setenv thicknum "300"
endif

if ( ${thick} == "eta2p0"  ) then
setenv thicknum "200"
endif

if ( ${thick} == "eta2p5"  ) then
setenv thicknum "120"
endif


foreach ncut  ($noisecuts)
echo "------------------------"
echo "Noise cut ${ncut}"

echo "python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/output_partGun_PDGid22_x10_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output --calcthickfactors $MODE --expectedthick $thicknum"

python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/output_partGun_PDGid22_x10_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output --calcthickfactors $MODE --expectedthick $thicknum

end 

end

