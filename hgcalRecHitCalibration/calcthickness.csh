#!/bin/tcsh

#cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/CMSSW_10_3_0_pre4/src
#cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/CMSSW_10_3_0_pre5/src
#cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/hitcalibrationV11/fcpermip/CMSSW_11_0_X_2019-09-23-2300/src
cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/hitcalibrationV11/fcpermip/CMSSW_11_0_X_2019-10-24-1100/src/
eval `scramv1 runtime -csh`
cd -

#====>>> IMPORTANT: To avoid accidental overwrite you should put by hand 
#the latest version of hgcalRecHitCalibration.py in the ntuples directory. 

#====>>> Also important: Look out for the name below "x100" part
#Run mode True and then False if you want to check the initial histo 
#EDIT: No need for false mode at the moment
setenv MODE "True"

# Thicknesses we are shooting
#setenv thicknesses "eta1p6 eta2p0 eta2p5" # 300, 200, 120 
#setenv thicknesses "CE_E_Front_120um CE_E_Front_200um CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um"
#setenv thicknesses "CE_H_Coarse_Scint"
#setenv thicknesses "CE_H_Fine_Scint"
setenv thicknesses "CE_H_Coarse_Scint_4285 CE_H_Coarse_Scint_4295 CE_H_Coarse_Scint_4305 CE_H_Coarse_Scint_4315 CE_H_Coarse_Scint_4325 CE_H_Coarse_Scint_4335 CE_H_Coarse_Scint_4345 CE_H_Coarse_Scint_4354 CE_H_Coarse_Scint_4364" 
#setenv thicknesses "CE_H_Coarse_Scint CE_H_Coarse_Scint_Var1 CE_H_Coarse_Scint_Var2"
#setenv thicknesses "CE_H_Fine_Scint CE_H_Fine_Scint_Var1 CE_H_Fine_Scint_Var2"
#setenv thicknesses "CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um"
#setenv thicknesses "CE_E_Front_120um CE_E_Front_200um"
#setenv thicknesses "Scint_R180_Z430 Scint_R80_Z430"
#setenv thicknesses "Scint_R80_Z430"
#setenv thicknesses "eta1p6 eta2p0" # 300, 200
#setenv thicknesses "eta2p0"
#setenv thicknesses "eta2p5" # 120
# Noise cut
setenv noisecuts "0 3 5 10"
#setenv noisecuts "0"

#Options are zerostage, firststage, secondstage
setenv calibstage "zerostage"
#if you are on first put the thickness factor and the corresponding thickness 
setenv firstthicknessfactor "0.771"
setenv firststagechoice "300"
#if you are on second stage put the thickness factor and the corresponding thickness 
setenv secondthicknessfactor "0.770"
setenv secondstagechoice "200"

# Starting the loop through all thicknesses
foreach thick ($thicknesses)
echo "===================================================================================="
echo "Thickness $thick"

if ( ${thick} == "eta1p6" || ${thick} == "CE_E_Front_300um" || ${thick} == "CE_H_Fine_300um" || ${thick} == "CE_H_Coarse_300um" ) then
setenv thicknum "300"
endif

if ( ${thick} == "eta2p0" || ${thick} == "CE_E_Front_200um"  || ${thick} == "CE_H_Fine_200um" ) then
setenv thicknum "200"
endif

if ( ${thick} == "eta2p5" || ${thick} == "CE_E_Front_120um"  || ${thick} == "CE_H_Fine_120um" ) then
setenv thicknum "120"
endif

if ( ${thick} == "CE_H_Coarse_Scint" || ${thick} == "CE_H_Fine_Scint" || ${thick} == "CE_H_Fine_Scint_Var1" || ${thick} == "CE_H_Fine_Scint_Var2" || ${thick} == "CE_H_Coarse_Scint_Var1" || ${thick} == "CE_H_Coarse_Scint_Var2" || ${thick} == "CE_H_Coarse_Scint_4285" || ${thick} == "CE_H_Coarse_Scint_4295" || ${thick} == "CE_H_Coarse_Scint_4305" || ${thick} == "CE_H_Coarse_Scint_4315" || ${thick} == "CE_H_Coarse_Scint_4325" || ${thick} == "CE_H_Coarse_Scint_4335" || ${thick} == "CE_H_Coarse_Scint_4345" || ${thick} == "CE_H_Coarse_Scint_4354" || ${thick} == "CE_H_Coarse_Scint_4364") then
setenv thicknum "-1"
endif

foreach ncut  ($noisecuts)
echo "------------------------"
echo "Noise cut ${ncut}"

if ( ${calibstage} == "zerostage" ) then 
echo "python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/output_partGun_PDGid22_x100_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output/${thick}/${calibstage} --calcthickfactors $MODE --expectedthick $thicknum  --region ${thick}"

python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/output_partGun_PDGid22_x100_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output/${thick}/${calibstage} --calcthickfactors $MODE --expectedthick $thicknum --region ${thick}

endif

if ( ${calibstage} == "firststage" ) then   
echo "python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/output_partGun_PDGid22_x100_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output/${thick}/${calibstage} --calcthickfactors $MODE --expectedthick $thicknum  --region ${thick} --firststage --firstthicknessfactor ${firstthicknessfactor} --firststagechoice ${firststagechoice}"

python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/output_partGun_PDGid22_x100_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output/${thick}/${calibstage} --calcthickfactors $MODE --expectedthick $thicknum --region ${thick} --firststage --firstthicknessfactor ${firstthicknessfactor} --firststagechoice ${firststagechoice}

endif

if ( ${calibstage} == "secondstage" ) then
echo "python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/output_partGun_PDGid22_x100_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output/${thick}/${calibstage} --calcthickfactors $MODE --expectedthick $thicknum  --region ${thick} --firststage --firstthicknessfactor ${firstthicknessfactor} --firststagechoice ${firststagechoice} --secondstage --secondthicknessfactor ${secondthicknessfactor} --secondstagechoice ${secondstagechoice}"

python /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration.py --input /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/output_partGun_PDGid22_x100_E60.0To60.0_NTUP_${thick}_ecut${ncut}.root --maxEvents -1 --ecut $ncut --verbosityLevel 2 --dependSensor True --output output.root --outDir output/${thick}/${calibstage} --calcthickfactors $MODE --expectedthick $thicknum --region ${thick} --firststage --firstthicknessfactor ${firstthicknessfactor} --firststagechoice ${firststagechoice} --secondstage --secondthicknessfactor ${secondthicknessfactor} --secondstagechoice ${secondstagechoice}

endif

end 

end

