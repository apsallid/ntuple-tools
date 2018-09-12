#!/bin/tcsh

# Thicknesses we are shooting
setenv thicknesses "eta1p6 eta2p0 eta2p5" # 300, 200, 120 
# Noise cut
setenv noisecuts "3 5 10"

rm voodoo voodoo1 voodoo2 voodoo3 
setenv eospath "/eos/cms/store/user/apsallid/HGCal"
setenv targetdirs "FlatRandomEGunProducer_apsallid_PDGId22_nPart1_E60"
setenv date "20180909"

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

setenv filelistpath "${eospath}/${targetdirs}_${thick}_${date}/NTUP"

foreach ncut  ($noisecuts)
echo "------------------------"
echo "Noise cut ${ncut}"

#Create local structure for the output
rm -rf ${thick}/${ncut}/output ${thick}/${ncut}/jobs ${thick}/${ncut}/logs
mkdir -p ${thick}/${ncut}/output ${thick}/${ncut}/jobs ${thick}/${ncut}/logs
chmod 755 -R ${thick}/${ncut}/output ${thick}/${ncut}/jobs ${thick}/${ncut}/logs

foreach file (`ls ${filelistpath}`)

cat setupRecHitCalibration.sh > voodoo

sed -e "s/THEECUT/$ncut/g" voodoo > voodoo1
sed -e "s/THETHICKNESS/$thick/g" voodoo1 > voodoo2
sed -e "s/THEINPUTFILE/$file/g" voodoo2 > voodoo3
sed -e "s/THETHICKNUM/${thicknum}/g" voodoo3 > voodoo4

mv voodoo4 ${thick}/${ncut}/jobs/rehitcalib_${file}.job
chmod 755 ${thick}/${ncut}/jobs/rehitcalib_${file}.job

echo rehitcalib_${file}.job

rm voodoo 
rm voodoo1
rm voodoo2
rm voodoo3

end

end

end
~

