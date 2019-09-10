#!/bin/bash

echo "The script starts now."

clusterid=${1}
procid=${2}
THEECUT=${3}
THETHICKNESS=${4}
THEINPUTFILE=${5}
THETHICKNUM=${6}
THECALIBSTAGE=${7}
THEFIRSTTHICKNESSFACTOR=${8}
THEFIRSTSTAGECHOICE=${9}
THESECONDTHICKNESSFACTOR=${10}
THESECONDSTAGECHOICE=${11}
echo "System: "
uname -a

cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/hitcalibrationV11/fcpermip/CMSSW_11_0_X_2019-09-04-2300/src
eval `scramv1 runtime -sh`
cd -

cp /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/*.py .
cp /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/hgcalRecHitCalibration.py . 

export PWD=`pwd`

export eospath="/eos/cms/store/user/apsallid/HGCal/Validation/Photons"
export targetdirs="CloseByParticleGunProducer_apsallid_PDGId22_nPart1_E60_eta1p4to4p0"
#The delta 2.5 is dummy. We do not shoot two particles. 
export date="Delta_2p5_20190905"

export eosmainpath="${eospath}/${targetdirs}_${THETHICKNESS}_${date}/NTUP"
#export eosmainpath="${eospath}/${targetdirs}_${THETHICKNESS}_${date}/NTUP_THICKNESS1"

#echo ${targetdirs} | grep CloseBy

#if [[ $? == 0 ]]; then 

#python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/${THEINPUTFILE} --maxEvents -1 --ecut ${THEECUT} --verbosityLevel 2 --dependSensor True --output output_${THEINPUTFILE} --outDir output --calcthickfactors False --expectedthick ${THETHICKNUM} --shootoneside --region ${THETHICKNESS}

#else 

if [ ${THECALIBSTAGE} == "zerostage" ]; then
echo "python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/${THEINPUTFILE} --maxEvents -1 --ecut ${THEECUT} --verbosityLevel 2 --dependSensor True --output output_${THEINPUTFILE} --outDir output --calcthickfactors False --expectedthick ${THETHICKNUM}"
python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/${THEINPUTFILE} --maxEvents -1 --ecut ${THEECUT} --verbosityLevel 2 --dependSensor True --output output_${THEINPUTFILE} --outDir output --calcthickfactors False --expectedthick ${THETHICKNUM} --region ${THETHICKNESS}

fi

if [ ${THECALIBSTAGE} == "firststage" ]; then
echo"======================================"
echo "python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/${THEINPUTFILE} --maxEvents -1 --ecut ${THEECUT} --verbosityLevel 2 --dependSensor True --output output_${THEINPUTFILE} --outDir output --calcthickfactors False --expectedthick ${THETHICKNUM} --firststage --firstthicknessfactor ${firstthicknessfactor} --firststagechoice ${firststagechoice}"
python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/${THEINPUTFILE} --maxEvents -1 --ecut ${THEECUT} --verbosityLevel 2 --dependSensor True --output output_${THEINPUTFILE} --outDir output --calcthickfactors False --expectedthick ${THETHICKNUM} --region ${THETHICKNESS} --firststage --firstthicknessfactor ${THEFIRSTTHICKNESSFACTOR} --firststagechoice ${THEFIRSTSTAGECHOICE}

fi

if [ ${THECALIBSTAGE} == "secondstage" ]; then
echo "python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/${THEINPUTFILE} --maxEvents -1 --ecut ${THEECUT} --verbosityLevel 2 --dependSensor True --output output_${THEINPUTFILE} --outDir output --calcthickfactors False --expectedthick ${THETHICKNUM} --firststage --firstthicknessfactor ${firstthicknessfactor} --firststagechoice ${firststagechoice} --secondstage --secondthicknessfactor ${secondthicknessfactor} --secondstagechoice ${secondstagechoice}"
python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/${THEINPUTFILE} --maxEvents -1 --ecut ${THEECUT} --verbosityLevel 2 --dependSensor True --output output_${THEINPUTFILE} --outDir output --calcthickfactors False --expectedthick ${THETHICKNUM} --region ${THETHICKNESS} --firststage --firstthicknessfactor ${THEFIRSTTHICKNESSFACTOR} --firststagechoice ${THEFIRSTSTAGECHOICE} --secondstage --secondthicknessfactor ${THESECONDTHICKNESSFACTOR} --secondstagechoice ${THESECONDSTAGECHOICE}

fi

cp output/output_${THEINPUTFILE} /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${THECALIBSTAGE}/${THETHICKNESS}/${THEECUT}/output/.



 
