# Perform the RecHit Calibration to obtain the thickness correction factors. 

# Example command of a basic setup: 
# python hgcalRecHitCalibration.py --input root://eoscms.cern.ch//eos/cms/store/user/apsallid/HGCal/FlatRandomEGunProducer_apsallid_PDGId22_nPart1_E60_eta1p6_20180904/NTUP/partGun_PDGid22_x10_E60.0To60.0_NTUP_25.root --maxEvents -1 --ecut 3 --verbosityLevel 2 --dependSensor True --output output.root --outDir output --calcthickfactors False

#---------------------------------------------------------------------------------------------------
#Necessary imports
from __future__ import print_function
import ROOT
import os,sys
import numpy as np
from HGCalImagingAlgo import *
from NtupleDataFormat import HGCalNtuple

#---------------------------------------------------------------------------------------------------
#Options from command line
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--input",help="The input hgcal ntuple file coming from reco-ntuples")
parser.add_option("--maxEvents",help="The maximum number of events you want to process. -1 for all events. ", type=int)
parser.add_option("--ecut",help="The rechit energy cut relative to noise",type=float)
parser.add_option("--verbosityLevel",help=" 0 - only basic info (default); 1 - additional info; 2 - detailed info printed, histograms produced",type=int)
parser.add_option("--dependSensor",default=True,help="Introduces a sensor dependened noise threshold. Dummy for now for thicknesscorrection. Look RecHitCalibration.py")
parser.add_option("--output",default="output.root",help="Name of output root file to keep the histos (default: %default)")
parser.add_option("--outDir",default="output",help="Output directory with all plots and files (default: %default)")
parser.add_option("--calcthickfactors",default=False,help="If True the fitting will be performed and the thickness correction factors will be calculated")
(options,args)=parser.parse_args()

#---------------------------------------------------------------------------------------------------
def getRecHitDetIds(rechits):
    recHitsList = []
    for rHit in rechits:
        recHitsList.append(rHit.detid())
    # print("RecHits -"*10)
    # print(recHitsList)
    recHits = np.array(recHitsList)
    return recHits

#---------------------------------------------------------------------------------------------------
def getHitList(simClus, recHitDetIds):
    sClusHitsList = []
    for DetId in simClus.hits():
        sClusHitsList.append(DetId)
    sClusHits = np.array(sClusHitsList)
    # thanks to http://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array
    recHitIndices = np.nonzero(np.in1d(recHitDetIds, sClusHits))
    return recHitIndices

#---------------------------------------------------------------------------------------------------
# get list of rechist associated to sim-cluster hits
def getRecHitsSimAssoc(rechits_raw, simcluster):
    # get sim-cluster associations
    nSimClus = 0
    simClusHitAssoc = []
    recHitDetIds = getRecHitDetIds(rechits_raw)
    for simClusIndex, simClus in enumerate(simcluster):
        simClusHitAssoc.append(getHitList(simClus, recHitDetIds))
        nSimClus += 1
    # get list of rechist associated to simhits
    rHitsSimAssoc = [[] for k in range(0,nSimClus)]
    for simClusIndex, simCl in enumerate(simcluster):
        if (options.verbosityLevel>=1): print( "Sim-cluster index: ",simClusIndex, ", pT: ",simCl.pt(), ", E: ",simCl.energy(), ", phi: ",simCl.phi(), ", eta: ",simCl.eta())
        # loop over sim clusters and then rechits
        rHitsSimAssocTemp = []
        for hitIndexArray in simClusHitAssoc[simClusIndex]:
            for hitIndex in hitIndexArray:
                thisHit = rechits_raw[hitIndex]
                if(not recHitAboveThreshold(thisHit, options.ecut, options.dependSensor)[1]): continue
                # independent of sim cluster, after cleaning
                rHitsSimAssocTemp.append(thisHit)
        rHitsSimAssoc[simClusIndex]= rHitsSimAssocTemp
    return rHitsSimAssoc

#---------------------------------------------------------------------------------------------------
# 1D histograming of given list of values
def histValue1D(fValues, histDict, tag = "hist1D_", title = "hist 1D", axunit = "a.u.", binsRangeList = [10, -1, 1], ayunit = "a.u."):
    # sanity check
    if (histDict == None): return

    # define event-level hists
    histDict[tag]  = ROOT.TH1F(tag, title+";"+axunit+";"+ayunit, binsRangeList[0], binsRangeList[1], binsRangeList[2])
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset()*1.5)
    # loop over all values
    if (options.verbosityLevel>=3): print( "tag: ", tag, ", fValues: ", fValues)
    for value in fValues:
        histDict[tag].Fill(value)
    return histDict

#---------------------------------------------------------------------------------------------------
# print/save all histograms
def histPrintSaveAll(histDict, outDir):
    imgType = "pdf"
    outfile = ROOT.TFile("{}/{}".format(outDir, options.output), "recreate")
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)
    if (options.verbosityLevel>=3): print( "histDict.items(): ", histDict.items())
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPadTopMargin(0.05)
        ROOT.gStyle.SetPadBottomMargin(0.12)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadRightMargin(0.02)
        if type(item) == ROOT.TH1F:
            item.Draw("hist0")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            item.Draw("colz")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    return

def main():
    inFile = options.input
    ntuple = HGCalNtuple(inFile)

    # init output stuff
    outDir = options.outDir
    if not os.path.exists(outDir): os.makedirs(outDir)
    histDict = {}


    if options.calcthickfactors == "False":

        # In each event, from the rechits associated to sim-cluster hits 
        # will calculate sum of energy
        rHitsSimAssocSumEoverEgen = []
        # The same but just applying the cleaning to all raw rechits. 
        rHitsCleanedSumEoverEgen = []

        # start event loop
        for event in ntuple:
            if event.entry() >= options.maxEvents and options.maxEvents != -1 : break
            if (options.verbosityLevel>=1): print( "\nCurrent event: ", event.entry())

            # get collections of raw rechits, sim clusters, caloparticles.
            recHitsRaw = event.recHits()
            simClusters = event.simClusters()
            caloParts = event.caloParticles()

            # get flat list of rechist associated to sim-cluster hits. Cleaning with 
            # ecut threshold inside function. 
            rHitsSimAssoc = getRecHitsSimAssoc(recHitsRaw, simClusters)
            # get flat list of raw rechits which satisfy treshold condition
            rHitsCleaned = [rechit for rechit in recHitsRaw if recHitAboveThreshold(rechit, options.ecut, options.dependSensor)[1]]

            #Energies
            # shooting also anti-particle
            rHitsSimAssocEp = [] # eta > 0
            rHitsSimAssocEm = [] # eta <= 0
            for simClusIndex in range(0,len(rHitsSimAssoc)):
                # loop over assoc. rec hits
                for thisHit in rHitsSimAssoc[simClusIndex]:
                    # for each sim cluster
                    if thisHit.eta() > 0.: rHitsSimAssocEp.append(thisHit.energy())
                    elif thisHit.eta() <= 0.: rHitsSimAssocEm.append(thisHit.energy())
                    else: print("You shouldn't be here.")

            #The sum of all rechits energies of the event that are associated to simclusters
            eventrHitsSimAssocSumEp = np.asarray(rHitsSimAssocEp).sum()
            eventrHitsSimAssocSumEm = np.asarray(rHitsSimAssocEm).sum()

            Egenp = 0.
            Egenm = 0.
            for cp in caloParts: 
                if cp.eta() > 0. : Egenp = cp.energy()
                if cp.eta() <= 0.: Egenm = cp.energy()
  
            rHitsSimAssocSumEoverEgen.append( eventrHitsSimAssocSumEp / Egenp )
            rHitsSimAssocSumEoverEgen.append( eventrHitsSimAssocSumEm / Egenm )
            print ("Egen plus  ", Egenp, " GeV " , " Sum(E_i) ", eventrHitsSimAssocSumEp, " sum(E_i)/Egen ", (eventrHitsSimAssocSumEp / Egenp))
            print ("Egen minus ", Egenm, " GeV " , " Sum(E_i) ", eventrHitsSimAssocSumEm, " sum(E_i)/Egen ", (eventrHitsSimAssocSumEm / Egenm))

        # histograms
        histDict = histValue1D(rHitsSimAssocSumEoverEgen, histDict, tag = "SumEoverEgen", title = "Reconstructed hits energy over generated energy",   axunit = "#sum E_{i}/E_{gen}",    binsRangeList = [200, 0, 2], ayunit = "N(events)")
        histPrintSaveAll(histDict, outDir)

    #Now in the calculation of the thickness correction factors. 
    if options.calcthickfactors == "True":
        #print("Not yet in here. ")
        finalmergedfile = ROOT.TFile("output_partGun_PDGid22_x10_E60.0To60.0_NTUP.root", "recreate")
        #The histo we want to fit 
        sumEoverEgen = finalmergedfile.Get("SumEoverEgen")
        #The canvas to plot the histo before fitting
        mycE = ROOT.TCanvas(outDir, outDir, 500, 500)
        sumEoverEgen.Draw("");

        meanE = sumEoverEgen.GetMean();
        rmsE  = sumEoverEgen.GetRMS();

        sumEoverEgen.Fit("gaus","LR0","", meanE-2*rmsE, meanE+2*rmsE);
        fitResult = sumEoverEgen.GetFunction("gaus");
 
if __name__ == '__main__':
    main()
