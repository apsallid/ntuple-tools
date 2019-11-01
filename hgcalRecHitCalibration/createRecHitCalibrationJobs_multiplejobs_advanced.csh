#!/bin/tcsh

#====>>> IMPORTANT: To avoid accidental overwrite you should put by hand 
#the latest version of hgcalRecHitCalibration.py in the ntuples directory. 

#====>>> Also important: Look out for the name below "x100" part

#It is advantageous to submit multiple jobs as a single cluster because:
# 1. Only one copy of the checkpoint file is needed to represent all jobs 
#    in a cluster until they begin execution. 
# 2. There is much less overhead involved for Condor to start the next job 
#    in a cluster than for Condor to start a new cluster. This can make a big 
#    difference when submitting lots of short jobs.

# So, these are short jobs and it is in any case not doable to run each one in 
# separate clusters. 

# Thicknesses we are shooting
#setenv thicknesses "eta1p6 eta2p0 eta2p5" # 300, 200, 120 
#setenv thicknesses "CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um"
#setenv thicknesses "CE_E_Front_300um CE_H_Coarse_300um CE_H_Coarse_Scint CE_H_Fine_120um"
#setenv thicknesses "CE_H_Coarse_Scint"
#setenv thicknesses "CE_H_Fine_Scint"
#setenv thicknesses "CE_H_Fine_Scint_Var1 CE_H_Fine_Scint_Var2 CE_H_Coarse_Scint_Var1 CE_H_Coarse_Scint_Var2"
#setenv thicknesses "CE_H_Fine_200um CE_H_Fine_300um"
#---------------
#setenv thicknesses "CE_E_Front_120um CE_E_Front_200um CE_E_Front_300um"
#setenv thicknesses "CE_H_Fine_120um CE_H_Fine_200um CE_H_Fine_300um"
#setenv thicknesses "CE_H_Coarse_Scint CE_H_Coarse_Scint_Var1 CE_H_Coarse_Scint_Var2"
setenv thicknesses "CE_H_Coarse_Scint_4285 CE_H_Coarse_Scint_4295 CE_H_Coarse_Scint_4305 CE_H_Coarse_Scint_4315 CE_H_Coarse_Scint_4325 CE_H_Coarse_Scint_4335 CE_H_Coarse_Scint_4345 CE_H_Coarse_Scint_4354 CE_H_Coarse_Scint_4364" 
#setenv thicknesses "CE_H_Fine_Scint CE_H_Fine_Scint_Var1 CE_H_Fine_Scint_Var2"
#---------------
#setenv thicknesses "eta1p6 eta2p0"  # 300, 200
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

#When running on multiple jobs you should give the number of input files in eos
#because this will the queue variable. Do it in wc style in the future.
#procid goes to N-1  

setenv alljobs "100"
#This is for the number of jobs per clusterid
set jobsperclusterchoice=10
#Since all jobs are 100 and we want 10 jobs per cluster we split 
#100 jobs in 10 batches. 
# SET BATCHES ON THE LOOP BY HAND

setenv eospath "/eos/cms/store/user/apsallid/HGCal/Validation/Photons"
setenv targetdirs "CloseByParticleGunProducer_apsallid_PDGId22_nPart1_E60_eta1p4to4p0"
#The delta 2.5 is dummy. We do not shoot two particles. 
#setenv date "Delta_2p5_20191007"
setenv date "Delta_2p5_20191101" #CE_H_Coarse_Scint_Var1 CE_H_Coarse_Scint_Var2

setenv PWD `pwd`

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

setenv filelistpath "${eospath}/${targetdirs}_${thick}_${date}/NTUP"
#setenv filelistpath "${eospath}/${targetdirs}_${thick}_${date}/NTUP_THICKNESS1"

if ( ${thick} == "CE_H_Coarse_Scint" || ${thick} == "CE_H_Fine_Scint" || ${thick} == "CE_H_Fine_Scint_Var1" || ${thick} == "CE_H_Fine_Scint_Var2" || ${thick} == "CE_H_Coarse_Scint_Var1" || ${thick} == "CE_H_Coarse_Scint_Var2" || ${thick} == "CE_H_Coarse_Scint_4285" || ${thick} == "CE_H_Coarse_Scint_4295" || ${thick} == "CE_H_Coarse_Scint_4305" || ${thick} == "CE_H_Coarse_Scint_4315" || ${thick} == "CE_H_Coarse_Scint_4325" || ${thick} == "CE_H_Coarse_Scint_4335" || ${thick} == "CE_H_Coarse_Scint_4345" || ${thick} == "CE_H_Coarse_Scint_4354" || ${thick} == "CE_H_Coarse_Scint_4364") then
setenv thicknum "-1"
endif



foreach ncut  ($noisecuts)
echo "------------------------"
echo "Noise cut ${ncut}"

#Create local structure for the output
rm -rf ${calibstage}/${thick}/${ncut}/output ${calibstage}/${thick}/${ncut}/jobs ${calibstage}/${thick}/${ncut}/logs
mkdir -p ${calibstage}/${thick}/${ncut}/output ${calibstage}/${thick}/${ncut}/jobs ${calibstage}/${thick}/${ncut}/logs
chmod 755 -R ${calibstage}/${thick}/${ncut}/output ${calibstage}/${thick}/${ncut}/jobs ${calibstage}/${thick}/${ncut}/logs

#foreach file (`ls ${filelistpath}`)
foreach batch (`seq 0 9`)

#cat setupRecHitCalibration.sh > voodoo

#sed -e "s/THEECUT/$ncut/g" voodoo > voodoo1
#sed -e "s/THETHICKNESS/$thick/g" voodoo1 > voodoo2
#sed -e "s/THEINPUTFILE/$file/g" voodoo2 > voodoo3
#sed -e "s/THETHICKNUM/${thicknum}/g" voodoo3 > voodoo4

#mv voodoo4 ${thick}/${ncut}/jobs/rehitcalib_${file}.job
echo '+JobFlavour = "tomorrow" ' > rechitcalib_$batch.sub
echo ' ' >> rechitcalib_$batch.sub
echo "executable  = ${PWD}/setupRecHitCalibration.sh" >> rechitcalib_$batch.sub
#echo "arguments   = "'$(ClusterID) $(ProcId)'" ${ncut} ${thick} ${file} ${thicknum} " >> rechitcalib_${file}.sub
echo "arguments   = "'$(ClusterID) $(ProcId)'" ${ncut} ${thick} "'$(infile)'" ${thicknum} ${calibstage} ${firstthicknessfactor} ${firststagechoice} ${secondthicknessfactor} ${secondstagechoice}" >> rechitcalib_$batch.sub
echo "output      = ${PWD}/${calibstage}/${thick}/${ncut}/logs/rehitcalib_"'$(infile)'".out " >> rechitcalib_$batch.sub
echo "error       = ${PWD}/${calibstage}/${thick}/${ncut}/logs/rehitcalib_"'$(infile)'".err " >> rechitcalib_$batch.sub
echo "log         = ${PWD}/${calibstage}/${thick}/${ncut}/logs/rehitcalib_"'$(infile)'"_htc.log " >> rechitcalib_$batch.sub
echo 'requirements = (OpSysAndVer =?= "CentOS7") ' >> rechitcalib_${batch}.sub
echo 'max_retries = 1' >> rechitcalib_$batch.sub

rm voodoo
touch voodoo
#foreach jobspercluster (`seq 1 ${jobsperclusterchoice}`)
foreach jobspercluster (`seq 1 10`)
set filesnum=`expr ${batch} \* 10  + ${jobspercluster} ` 
echo -n "closeby_PDGid22_x100_E60.0To60.0_NTUP_${filesnum}.root " >> voodoo
end 

setenv batchfilelist `cat voodoo`
echo "queue infile in (${batchfilelist}) " >> rechitcalib_$batch.sub

mv rechitcalib_$batch.sub ${calibstage}/${thick}/${ncut}/jobs/rehitcalib_$batch.sub
chmod 755 ${calibstage}/${thick}/${ncut}/jobs/rehitcalib_$batch.sub

#echo rechitcalib_${file}.sub

#rm voodoo 
#rm voodoo1
#rm voodoo2
#rm voodoo3

end

end

end
~

