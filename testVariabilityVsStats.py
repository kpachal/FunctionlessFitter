import sys
import os
import ROOT
import scipy
from decimal import *
import argparse

from HistWrapper import Dataset
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual
from PseudoExperimenter import PseudoExperimenter
from MathFunctions import makeHistFromVector

class RunFitter :

  def __init__(self) :
    # Make a fitter
    self.myFitter = FunctionlessFitter()
    self.myFitter.flatStartVal = 1.0

    self.setICHEPValues()

    self.parseCommandLineArgs()

    if self.args.dataset :
      if self.args.dataset == "EOYE" :
        self.setEOYEValues()
      elif self.args.dataset == "TLA" :
        self.setTLAValues()
      elif self.args.dataset == "TLAFull" :
        self.setTLAFull()
      elif self.args.dataset == "ICHEP" :
        self.setICHEPValues()
      else :
        print "Unrecognized dataset!"
        return -1

  def setEOYEValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/dataLikeHistograms_IBLOff.2015.root")
    self.nominalFitFile = ROOT.TFile("samples/Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","READ")
    self.hist = self.infile.Get("Nominal/mjj_Data_2015_3p57fb")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(1101)
    self.binHigh = -1

    self.outputFileName = "results/test/outputfile_testVariabilityVsStats_EOYE.root"

    #self.myFitter.derivativeConstraints = {0:-1, 1:1}
    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}

    self.myFitter.flatStartVal = 1.0

  def setTLAValues(self) :

    # Get a histogram we want to fit. This is TLA.
    self.infile = ROOT.TFile("samples/mjj_fullDataset.root")
    self.nominalFitFile = ROOT.TFile("samples/SearchResultData_UA2_fullDataset_from443.root","READ")
    self.hist = self.infile.Get("mjj")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(444)
    self.binHigh = self.hist.FindBin(1220) #1252

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1, 4:-1, 5:1, 6:-1, 7:1}

    self.myFitter.flatStartVal = 1.0
    
    self.outputFileName = "results/test/outputfile_testVariabilityVsStats_TLA_constTo6.root"

  def setTLAFull(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/mjj_fullDataset.root")
    self.nominalFitFile = ""
    self.hist = self.infile.Get("mjj")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(444)
    self.binHigh = -1

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}

    self.myFitter.flatStartVal = 1.0
    
    self.outputFileName = "results/test/outputfile_testVariabilityVsStats_TLAFull.root"

  def setICHEPValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/dataLikeHistograms.2016_DS2p1_Resonance_Fixed.root")
    self.nominalFitFile = ROOT.TFile("samples/Step1_SearchPhase_mjj_Data_2016_15p7fb.root","READ")
    self.hist = self.infile.Get("Nominal/mjj_Data_2016_15p7fb")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(1100)
    self.binHigh = -1

    self.myFitter.flatStartVal = 1.0

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}

    self.outputFileName = "results/test/outputfile_testVariabilityVsStats_ICHEP.root"
  
  def executeFit(self) :

    wInput = Dataset(self.hist)

    # Define various start values.
    startVals = "data"
    
    resultList = []
    
    # Compute start values for user based ones
    data = []
    startValHist = self.hist.Clone("startValHist")
    startValHist.SetName("startValHist")
    getStartVals = Dataset(startValHist)
    if self.binHigh < 0 :
      fullVals, xvals, widths, edges, w1, w2 = wInput.getSelectedBinInfo(self.binLow,getStartVals.lastBinWithData)
    else :
      fullVals, xvals, widths, edges, w1, w2 = wInput.getSelectedBinInfo(self.binLow,self.binHigh)
    index = -1
    for item in fullVals :
      index = index+1
      data.append(item/widths[index])
    
    #self.myFitter.startValFormat = "linear"
    self.myFitter.startValFormat = "user"
    self.myFitter.userStartVals = data

    outputFile = ROOT.TFile(self.outputFileName,"RECREATE")
    outputFile.cd()

    self.myFitter.nPEs = 100
    result = self.myFitter.fit(wInput,self.binLow,self.binHigh)#,errType = "Bootstrap")
    result.SetDirectory(0)
    result.SetName("result")
    result.SetTitle("result")

    # Write everything to a file
    self.hist.SetName("data")
    self.hist.SetTitle("data")
    self.hist.Write()
    result.Write()
    
    # Get the residual of the hist
    residual = getResidual(self.hist,result,self.myFitter.rangeLow,self.myFitter.rangeHigh)
    residual.SetDirectory(0)
    residual.SetName("residual")
    residual.SetTitle("residual")
    residual.Write()

    # Evaluate distribution of residuals.
    residualBins = []
    for bin in range(self.myFitter.rangeLow,self.myFitter.rangeHigh+1) :
      residualBins.append(residual.GetBinContent(bin))
    resDistHist = makeHistFromVector(residualBins,1.0)
    resDistHist.SetDirectory(0)
    resDistHist.Fit("gaus","0")
    fittedGauss = resDistHist.GetFunction("gaus")
    gausmean = fittedGauss.GetParameter(1)
    gauswidth = fittedGauss.GetParameter(2)
    histmean = resDistHist.GetMean()
    histRMS = resDistHist.GetRMS()
    print "Examined distribution of residuals."
    print "Gaussian mean and width were",gausmean,gauswidth
    resDistHist.SetName("pullDistribution")
    resDistHist.SetTitle("pullDistribution")
    resDistHist.Write()
    fittedGauss.Write("gausFitToResidual")

    # Compare to nominal fit result from the old code, if available
    if self.nominalFitFile :
      nominalResidual = self.nominalFitFile.Get("residualHist")
      nominalResidual.SetDirectory(0)
      nominalResidual.SetName("nominalResidual")
      nominalResidual.SetTitle("nominalResidual")
      nominalResidual.Write()

      # Evaluate distribution of residuals.
      residualBins = []
      for bin in range(self.myFitter.rangeLow,self.myFitter.rangeHigh+1) :
        residualBins.append(nominalResidual.GetBinContent(bin))
      nomresDistHist = makeHistFromVector(residualBins,1.0)
      nomresDistHist.SetDirectory(0)
      nomresDistHist.Fit("gaus","0")
      nomfittedGauss = nomresDistHist.GetFunction("gaus")
      nomgausmean = nomfittedGauss.GetParameter(1)
      nomgauswidth = nomfittedGauss.GetParameter(2)
      nomhistmean = nomresDistHist.GetMean()
      nomhistRMS = nomresDistHist.GetRMS()
      print "By comparison, in nominal distribution of residuals:."
      print "Gaussian mean and width were",nomgausmean,nomgauswidth
      nomresDistHist.SetName("nominalPullDistribution")
      nomresDistHist.SetTitle("nominalPullDistribution")
      nomresDistHist.Write()
      nomfittedGauss.Write("gausFitToNominalResidual")

    outputFile.Close()
    print "Made file",self.outputFileName

  def parseCommandLineArgs(self) :
  
    parser = argparse.ArgumentParser(description='Specify config and overwrite values if desired.')
    parser.add_argument('--dataset')
    self.args = parser.parse_args()

if __name__ == "__main__":

  fitter = RunFitter()
  fitter.executeFit()

