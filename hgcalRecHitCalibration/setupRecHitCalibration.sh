#!/bin/bash

echo "The script starts now."

echo "System: "
uname -a


cd /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/CMSSW_10_3_0_pre4/src
eval `scramv1 runtime -sh`
cd -

cp /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/*.py .
cp /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/hgcalRecHitCalibration.py . 

export PWD=`pwd`

export eospath="/eos/cms/store/user/apsallid/HGCal"
export targetdirs="FlatRandomEGunProducer_apsallid_PDGId22_nPart1_E60"
export date="20180909"

export eosmainpath="${eospath}/${targetdirs}_THETHICKNESS_${date}/NTUP"


python hgcalRecHitCalibration.py --input root://eoscms.cern.ch/${eosmainpath}/THEINPUTFILE --maxEvents -1 --ecut THEECUT --verbosityLevel 2 --dependSensor True --output output_THEINPUTFILE --outDir output --calcthickfactors False --expectedthick THETHICKNUM

cp output/output_THEINPUTFILE /afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/THETHICKNESS/THEECUT/output/.



 
