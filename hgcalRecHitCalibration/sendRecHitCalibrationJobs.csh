#!/bin/tcsh

# Thicknesses we are shooting
#setenv thicknesses "eta1p6 eta2p0 eta2p5" # 300, 200, 120 
#setenv thicknesses "CE_E_Front_120um CE_E_Front_200um CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um"
#setenv thicknesses "CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um"
#setenv thicknesses "CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um"
#setenv thicknesses "CE_H_Coarse_Scint"
#setenv thicknesses "CE_H_Fine_Scint"
#setenv thicknesses "CE_H_Fine_Scint_Var1 CE_H_Fine_Scint_Var2 CE_H_Coarse_Scint_Var1 CE_H_Coarse_Scint_Var2"
setenv thicknesses "CE_H_Fine_Scint CE_H_Fine_Scint_Var1 CE_H_Fine_Scint_Var2"
#setenv thicknesses "CE_H_Fine_200um CE_H_Fine_300um"
#setenv thicknesses "CE_E_Front_120um CE_E_Front_200um"
#setenv thicknesses "eta1p6 eta2p0"
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

setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/PFCalStudies/CMS-HGCAL/calib/ntuple-tools/hgcalRecHitCalibration/${calibstage}/${thick}/${ncut}/jobs"

setenv runumberslist ` ls ${workpath} | grep .sub `

foreach run  ($runumberslist)

#echo ${run}
chmod 755 ${workpath}/${run}

echo "Sending ${run}"
#bsub -q 8nh -o /tmp/junk ${workpath}/${run}
condor_submit ${workpath}/${run}
echo "condor_submit ${workpath}/${run}"
#bsub -q 8nh -o ${workpath}/../logs/${run}.txt ${workpath}/${run}
#echo "bsub -q 8nh -o /tmp/junk ${workpath}/${run}"

end

end 

end




