from RecHitCalibration import *
import numpy as np
import pandas as pd

RecHitCalib = RecHitCalibration()

thickIndex = [0,1,2] # 120, 200, 300 um

sigma_noiseMeV = np.zeros(shape=(53,3)) # 52 layers by 3 thicknesses in the order 120, 200, 300 um
sigma_noiseMIP = np.zeros(shape=(53,3)) # 52 layers by 3 thicknesses in the order 120, 200, 300 um
MeVperMIP = np.zeros(shape=(53,3)) # 52 layers by 3 thicknesses in the order 120, 200, 300 um
MIPperGeV = np.zeros(shape=(53,3)) # 52 layers by 3 thicknesses in the order 120, 200, 300 um

#First is zero so put it by hand
sigma_noiseMeV[1] = [0,0,0]
sigma_noiseMIP[1] = [0,0,0]
MeVperMIP[1] = [0,0,0]
MIPperGeV[1] = [0,0,0]

#Loop over layers #starting from 2 since 1 is the zero in dedxweights
for layer in range(1,53): 
    sigma_noiseMeV[layer] = [RecHitCalib.sigmaNoiseMeV(layer,0), RecHitCalib.sigmaNoiseMeV(layer,1), RecHitCalib.sigmaNoiseMeV(layer,2)]
    sigma_noiseMIP[layer] = [RecHitCalib.sigmaNoiseMIP(layer,0), RecHitCalib.sigmaNoiseMIP(layer,1), RecHitCalib.sigmaNoiseMIP(layer,2)]
    MeVperMIP[layer] = [RecHitCalib.MeVperMIP(layer,0), RecHitCalib.MeVperMIP(layer,1), RecHitCalib.MeVperMIP(layer,2)]
    MIPperGeV[layer] = [RecHitCalib.MIPperGeV(layer,0), RecHitCalib.MIPperGeV(layer,1), RecHitCalib.MIPperGeV(layer,2)]
    

datasetsigma_noiseMeV = pd.DataFrame(data=sigma_noiseMeV, columns = ["120 um","200 um","300 um"] )
datasetsigma_noiseMIP = pd.DataFrame(data=sigma_noiseMIP, columns = ["120 um","200 um","300 um"] )
datasetMeVperMIP = pd.DataFrame(data=MeVperMIP, columns = ["120 um","200 um","300 um"] )
datasetMIPperGeV = pd.DataFrame(data=MIPperGeV, columns = ["120 um","200 um","300 um"] )
frames = [datasetsigma_noiseMeV, datasetsigma_noiseMIP, datasetMeVperMIP, datasetMIPperGeV]

dataset = pd.concat(frames, axis=1)

print(dataset)
