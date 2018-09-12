# Perform the RecHit Calibration to obtain the thickness correction factors. 

# Example command of a basic setup: 
# python hgcalRecHitCalibration.py --input root://eoscms.cern.ch//eos/cms/store/user/apsallid/HGCal/FlatRandomEGunProducer_apsallid_PDGId22_nPart1_E60_eta1p6_20180904/NTUP/partGun_PDGid22_x10_E60.0To60.0_NTUP_25.root --maxEvents -1 --ecut 3 --verbosityLevel 2 --dependSensor True --output output.root --outDir output --calcthickfactors False --expectedthick 300

#---------------------------------------------------------------------------------------------------
#Necessary imports
from __future__ import print_function
import ROOT
import math
import os,sys
import numpy as np
from HGCalImagingAlgo import *
from NtupleDataFormat import HGCalNtuple
from hgcalHistHelpers import *

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
parser.add_option("--expectedthick",help="The silicon thickness of the cells you wish and expect to have", type=int)
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
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TProfile2D:
            item.Draw("")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    return

#---------------------------------------------------------------------------------------------------
def profValues2D(fValues, histDict, tag="hist2D_", title="hist 2D", axunit="a.u.", binsBoundariesX=[10, -1, 1], ayunit="a.u.", binsBoundariesY=[10, -1, 1], binsBoundariesZ=[-1., 1], weighted2D=False, weight=1, verbosityLevel=0):
    """2D profile of given list of values"""
    # sanity check for hists
    if histDict is None:
        return
    # sanity check for boundaries
    if (len(binsBoundariesX) != len(binsBoundariesY)):
        return
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2):
        return
    # define event-level hists
    elif len(binsBoundariesX) == 3:  # bondaries in format [nbins, low, high]
        histDict[tag] = ROOT.TProfile2D(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2], binsBoundariesY[0], binsBoundariesY[1], binsBoundariesY[2], binsBoundariesZ[0], binsBoundariesZ[1])
    # set some properties
    histDict[tag].GetXaxis().SetTitleOffset(histDict[tag].GetXaxis().GetTitleOffset() * 1.0)
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset() * 3.0)
    # loop over all values
    if (verbosityLevel >= 3):
        print ("tag: ", tag, ", fValues: ", fValues)
    if (not weighted2D):
        for (valueX, valueY, valueZ) in fValues:
            print (valueX, valueY, valueZ)
            histDict[tag].Fill(valueX, valueY, valueZ)
    else:
        for (valueX, valueY, valueZ) in fValues:
            print (valueX, valueY, valueZ)
            histDict[tag].Fill(valueX, valueY, valueZ, weight)
    return histDict

#---------------------------------------------------------------------------------------------------
def acustompalette():
    NRGBs = 7
    NCont = 100
    ncolors = array('i', [])
    ROOT.gStyle.SetNumberContours(NCont);
    stops   = [ 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 ]
    red     = [ 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 ]
    green   = [ 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 ]
    blue    = [ 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 ]
    stopsArray = array('d', stops)
    redArray = array('d', red)
    greenArray = array('d', green)
    blueArray = array('d', blue)
    first_color_number = ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont);
    ROOT.gStyle.SetNumberContours(NCont)

    palsize = NCont
    palette = []
    for i in range(palsize):
        palette.append(first_color_number+i)
        palarray = array('i',palette)

    ROOT.gStyle.SetPalette(palsize,palarray)




def main():

    # init output stuff
    outDir = options.outDir
    if not os.path.exists(outDir): os.makedirs(outDir)
    histDict = {}

    EtaBoundaries = []
    if options.expectedthick == 300: 
        EtaBoundaries = [80, 1.2, 2.0]
    if options.expectedthick == 200: 
        EtaBoundaries = [160, 1.2, 2.8]
    if options.expectedthick == 120: 
        EtaBoundaries = [200, 1.5, 3.5]

    print (EtaBoundaries)


    if options.calcthickfactors == "False":

        inFile = options.input
        ntuple = HGCalNtuple(inFile)

        # In each event, from the rechits associated to sim-cluster hits 
        # will calculate sum of energy
        rHitsSimAssocSumEoverEgen = []
        # The same but just applying the cleaning to all raw rechits. 
        rHitsCleanedSumEoverEgen = []
        #For the 2D maps
        RechitRvsEtavsThickness = []
        RechitRvsLayervsThickness = []
        #For the rechits thickenesses
        rHitthickallmatched = []
        rHitthicktherest = []
        rHitthick120 = []
        rHitthick200 = []
        rHitthick300 = []

        # start event loop
        for event in ntuple:
            if event.entry() >= options.maxEvents and options.maxEvents != -1 : break
            if (options.verbosityLevel>=1): print( "\nCurrent event: ", event.entry())

            # We will only look for unconverted photons 
            genParticles = event.genParticles()
            skipEvent = False
            for particle in genParticles:
                if particle.reachedEE() != 2:
                    print("particle converted before reaching EE -- skipping the event!!")
                    skipEvent = True
                    break
            if skipEvent: continue

            # get collections of raw rechits, sim clusters, caloparticles.
            recHitsRaw  = event.recHits()
            simClusters = event.simClusters()
            caloParts   = event.caloParticles()
 
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
                    #------------------------------------------------------------------
                    # Reproducing Arabella's plots
                    # Rechit R vs Eta vs thickness
                    # If there are also hits in scintillator, they have enormous values. We do not 
                    # throw them away but instead put them with a value of 400. 
                    if thisHit.thickness() > 10000. : 
                        RechitRvsEtavsThickness.append( ( abs(thisHit.eta()), math.sqrt( thisHit.x()**2. + thisHit.y()**2. ), 400. ) )
                        RechitRvsLayervsThickness.append( ( thisHit.layer(), math.sqrt( thisHit.x()**2. + thisHit.y()**2. ), 400. ) )
                    else:
                        RechitRvsEtavsThickness.append( ( abs(thisHit.eta()), math.sqrt( thisHit.x()**2. + thisHit.y()**2. ), thisHit.thickness() ) )
                        RechitRvsLayervsThickness.append( ( thisHit.layer(), math.sqrt( thisHit.x()**2. + thisHit.y()**2. ), thisHit.thickness() ) )
                    #print ( math.sqrt( thisHit.x()**2. + thisHit.y()**2. ), thisHit.eta(), thisHit.thickness() )
                    #To find what percent of hits was good.
                    rHitthickallmatched.append( 1. )
                    if(thisHit.thickness() > 119. and thisHit.thickness() < 121.):    rHitthick120.append( 1. )
                    elif(thisHit.thickness() > 199. and thisHit.thickness() < 201.):  rHitthick200.append( 1. )
                    elif(thisHit.thickness() > 299. and thisHit.thickness() < 301.):  rHitthick300.append( 1. )
                    else: rHitthicktherest.append( 1. )
                    #------------------------------------------------------------------
                    # for each sim cluster
                    if thisHit.eta() > 0.: rHitsSimAssocEp.append(thisHit.energy())
                    elif thisHit.eta() <= 0.: rHitsSimAssocEm.append(thisHit.energy())
                    else: print("You shouldn't be here.")
                    #------------------------------------------------------------------

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
        histDict = histValue1D(rHitsSimAssocSumEoverEgen, histDict, tag = "SumEoverEgen", title = "Reconstructed hits energy over generated energy",   axunit = "#sum E_{i}/E_{gen}",    binsBoundariesX = [200, 0, 2], ayunit = "N(events)", verbosityLevel=options.verbosityLevel)
        histDict = profValues2D(RechitRvsEtavsThickness, histDict, tag = "RvsEtavsThickness", title = "R vs Eta vs Thickness", axunit = "|#eta|", binsBoundariesX = EtaBoundaries, ayunit = "R (cm)", binsBoundariesY=[100, 0., 300.], binsBoundariesZ=[0.,400.], weighted2D=False, verbosityLevel=options.verbosityLevel)
        histDict = profValues2D(RechitRvsLayervsThickness, histDict, tag = "RvsLayervsThickness", title = "R vs Layer vs Thickness", axunit = "Layer", binsBoundariesX = [53, 0., 53.], ayunit = "R (cm)", binsBoundariesY=[100, 0., 300.], binsBoundariesZ=[0.,400.], weighted2D=False, verbosityLevel=options.verbosityLevel)
        #Checking hit thicknesses
        histDict = histValue1D(rHitthick120, histDict, tag = "rHitthick120", title = "Matched reconstructed hits 120 #{mu}m",   axunit = "#hits in 120 #{mu}m",    binsBoundariesX = [2, 0, 2], ayunit = "N(hits)", verbosityLevel=options.verbosityLevel)
        histDict = histValue1D(rHitthick200, histDict, tag = "rHitthick200", title = "Matched reconstructed hits 200 #{mu}m",   axunit = "#hits in 200 #{mu}m",    binsBoundariesX = [2, 0, 2], ayunit = "N(hits)", verbosityLevel=options.verbosityLevel)
        histDict = histValue1D(rHitthick300, histDict, tag = "rHitthick300", title = "Matched reconstructed hits 300 #{mu}m",   axunit = "#hits in 300 #{mu}m",    binsBoundariesX = [2, 0, 2], ayunit = "N(hits)", verbosityLevel=options.verbosityLevel)
        histDict = histValue1D(rHitthickallmatched, histDict, tag = "rHitthickallmatched", title = "All reconstructed hits matched from simhits from simcluters",   axunit = " ",    binsBoundariesX = [2, 0, 2], ayunit = "N(hits)", verbosityLevel=options.verbosityLevel)
        histDict = histValue1D(rHitthicktherest, histDict, tag = "rHitthicktherest", title = "All the rest matched reconstructed hits without a logical thickness",   axunit = " ",  binsBoundariesX = [2, 0, 2], ayunit = "N(hits)", verbosityLevel=options.verbosityLevel)


        histPrintSaveAll(histDict, outDir)

    #Now in the calculation of the thickness correction factors. 
    if options.calcthickfactors == "True":
        ROOT.gStyle.SetOptStat(0);
        
        finalmergedfile = ROOT.TFile(options.input, "open")
        #The histo we want to fit 
        sumEoverEgen = finalmergedfile.Get("SumEoverEgen")
        sumEoverEgen.Rebin(4)

        #The profile to plot hit maps
        RvsEtavsThickness = finalmergedfile.Get("RvsEtavsThickness")
        RvsLayervsThickness = finalmergedfile.Get("RvsLayervsThickness")

        #The histos to check what thicknesses we got. 
        rHitthickallmatched = finalmergedfile.Get("rHitthickallmatched")
        rHitthick120 = finalmergedfile.Get("rHitthick120")
        rHitthick200 = finalmergedfile.Get("rHitthick200")
        rHitthick300 = finalmergedfile.Get("rHitthick300")
        rHitthicktherest = finalmergedfile.Get("rHitthicktherest")

        if options.expectedthick == 300: thickhisto = rHitthick300
        elif options.expectedthick == 200: thickhisto = rHitthick200
        elif options.expectedthick == 120: thickhisto = rHitthick120


        finaloutput = ROOT.TFile("finaloutput_thick%d_ecut%d.root"%(options.expectedthick,options.ecut), "recreate")
        #The canvas to plot the histo before fitting
        mycE = ROOT.TCanvas("P22E60_thick%d_ecut%d"%(options.expectedthick,options.ecut), "P22E60_thick%d_ecut%d"%(options.expectedthick,options.ecut), 500, 500)
        sumEoverEgen.Draw("");

        meanE = sumEoverEgen.GetMean();
        rmsE  = sumEoverEgen.GetRMS();

        sumEoverEgen.Fit("gaus","LR0","", meanE-0.5*rmsE, meanE+1.5*rmsE);
        fitResult = sumEoverEgen.GetFunction("gaus");

        fitResult.SetLineColor(2);
        meanFitE = fitResult.GetParameter(1);
        rmsFitE = fitResult.GetParameter(2);
        meanFitEerr = fitResult.GetParError(1);
        rmsFitEerr = fitResult.GetParError(2);

        fitResult.Draw("same")
        
        latx =  sumEoverEgen.GetXaxis().GetXmin()+(sumEoverEgen.GetXaxis().GetXmax()-sumEoverEgen.GetXaxis().GetXmin())/20.
        laty =  sumEoverEgen.GetMaximum()

        lat = ROOT.TLatex()
        lat.DrawLatex(latx,laty*0.9, "<#sum E_{i}/E_{gen}> = %3.3f +/- %3.3f"%(fitResult.GetParameter(1),fitResult.GetParError(1))    );
        lat.DrawLatex(latx,laty*0.8, "RMSfit = %3.3f +/- %3.3f"%(fitResult.GetParameter(2),fitResult.GetParError(2))   );
        lat.DrawLatex(latx,laty*0.7, "RMS/meanfit = %3.3f"%(fitResult.GetParameter(2)/fitResult.GetParameter(1))   );
        lat.DrawLatex(latx,laty*0.6, "#chi^{2}/N = %3.3f/%d = %3.3f"%(fitResult.GetChisquare(),fitResult.GetNDF(),fitResult.GetChisquare()/fitResult.GetNDF()) )
        lat.DrawLatex(latx,laty*0.5, "S/N = %d "% options.ecut    );
        lat.DrawLatex(latx,laty*0.4, "# hits %d #mum = %3.3f %%"%(options.expectedthick,(thickhisto.GetEntries() * 100.) / rHitthickallmatched.GetEntries() ) );

        mycE.Update()

        mycE.SaveAs("P22E60_thick%d_ecut%d.png"%(options.expectedthick,options.ecut))
        mycE.Write()
        
        #----------------------------------------------------------
        #RvsEtavsThickness
        mycE1 = ROOT.TCanvas("P22E60_RvsEtavsThickness_thick%d_ecut%d"%(options.expectedthick,options.ecut), "P22E60_RvsEtavsThickness_thick%d_ecut%d"%(options.expectedthick,options.ecut), 500, 500)

        acustompalette()
        ex1 = ROOT.TExec("ex1","acustompalette();");
        ex1.Draw();

        RvsEtavsThickness.Draw("COLZ")
        #RvsEtavsThickness.GetXaxis().SetTitleOffset(0.95) 
        mycE1.Update()

        palette = RvsEtavsThickness.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.__class__ = ROOT.TPaletteAxis
            palette.SetX1NDC(0.85)
            palette.SetX2NDC(0.9)
            #palette.SetY1NDC(0.1)
            #palette.SetY2NDC(0.6)
            palette.GetAxis().SetTickSize(.01)
            palette.GetAxis().SetTitle("Si thick")
            palette.GetAxis().SetTitleOffset(0.8);
            #palette.GetAxis().LabelsOption("v")
            ROOT.gPad.Update()

        mycE1.SaveAs("P22E60_RvsEtavsThickness_thick%d_ecut%d.png"%(options.expectedthick,options.ecut))
        mycE1.Write()
        #----------------------------------------------------------
        #RvsLayervsThickness
        mycE2 = ROOT.TCanvas("P22E60_RvsLayervsThickness_thick%d_ecut%d"%(options.expectedthick,options.ecut), "P22E60_RvsLayervsThickness_thick%d_ecut%d"%(options.expectedthick,options.ecut), 500, 500)

        acustompalette()
        ex1 = ROOT.TExec("ex1","acustompalette();");
        ex1.Draw();

        RvsLayervsThickness.Draw("COLZ")
        mycE2.Update()

        palette = RvsLayervsThickness.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.__class__ = ROOT.TPaletteAxis
            palette.SetX1NDC(0.85)
            palette.SetX2NDC(0.9)
            #palette.SetY1NDC(0.1)
            #palette.SetY2NDC(0.6)
            palette.GetAxis().SetTickSize(.01)
            palette.GetAxis().SetTitle("Si thick")
            palette.GetAxis().SetTitleOffset(0.8);
            #palette.GetAxis().LabelsOption("v")
            ROOT.gPad.Update()

        mycE2.SaveAs("P22E60_RvsLayervsThickness_thick%d_ecut%d.png"%(options.expectedthick,options.ecut))
        mycE2.Write()

        #finaloutput.Close()
        
 
if __name__ == '__main__':
    main()
