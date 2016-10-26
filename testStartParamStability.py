import sys
import os
import ROOT
import scipy
from decimal import *
import argparse

from HistWrapper import WrappedHist
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual
from BumpHunter import BumpHunter
from Chi2Test import Chi2Test
from LogLikelihoodTest import LogLikelihoodTest
from PseudoExperimenter import PseudoExperimenter

class RunFitter :

  def __init__(self) :
    # Make a fitter
    self.myFitter = FunctionlessFitter()
    self.myFitter.flatStartVal = 1.0

    #self.setEOYEValues()
    #self.setTLAValues()
    self.setTLAFull()
    #self.setICHEPValues()

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

    self.outputFileName = "results/test/outputfile_testStartParamStability_EOYE.root"

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
    #self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}

    self.myFitter.flatStartVal = 1.0
    
    self.outputFileName = "results/test/outputfile_testStartParamStability_TLA_upTo7.root"

  def setTLAFull(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/mjj_fullDataset.root")
    self.nominalFitFile = ""
    self.hist = self.infile.Get("mjj")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(444)
    self.binHigh = -1

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1, 4:-1, 5:1}
    #self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1}

    self.myFitter.flatStartVal = 1.0

    self.outputFileName = "results/test/outputfile_testStartParamStability_TLAFull_upTo5.root"

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

    self.outputFileName = "results/test/outputfile_testStartParamStability_ICHEP.root"
  
  def executeFit(self) :

    wInput = WrappedHist(self.hist)

    # Working: "data","linear","prelimFit"
    startVals = ["data"]#,"dataP5","exp","flat","linear","prelimFit"]
    
    # Make a bump hunter
    bumpHunter = BumpHunter()
    bumpHunter.minBinsInBump = 2
    bumpHunter.useSidebands = False
    
    resultList = []
    
    # Compute start values for user based ones
    startValDict = {}
    data = []
    datap5 = []
    startValHist = self.hist.Clone("startValHist")
    startValHist.SetName("startValHist")
    getStartVals = WrappedHist(startValHist)
    if self.binHigh < 0 :
      fullVals, xvals, widths, edges, w1, w2 = wInput.getSelectedBinInfo(self.binLow,getStartVals.lastBinWithData)
    else :
      fullVals, xvals, widths, edges, w1, w2 = wInput.getSelectedBinInfo(self.binLow,self.binHigh)
    index = -1
    for item in fullVals :
      index = index+1
      data.append(item/widths[index])
      datap5.append(Decimal(1.05)*(item/widths[index]))
    startValDict["data"] = data
    startValDict["dataP5"] = datap5
    
    outputFile = ROOT.TFile(self.outputFileName,"RECREATE")
    outputFile.cd()

    # Loop over start values.
    for valType in startVals :
    
      if valType == "exp" or valType == "flat" or valType == "linear" :
        self.myFitter.startValFormat = valType
      elif "data" in valType :
        self.myFitter.startValFormat = "user"
        self.myFitter.userStartVals = startValDict[valType]
      elif "prelimFit" in valType :
        # take result of previous fit as input
        self.myFitter.startValFormat = "user"
        self.myFitter.userStartVals = self.myFitter.result
      else :
        print "Unknown start value specified!"
        return
    
      result = self.myFitter.fit(wInput,self.binLow,self.binHigh)
      result.SetDirectory(0)
      result.SetName("result_"+valType)
      result.SetTitle("result_"+valType)
      resultList.append(result)

      # Write everything to a file
      self.hist.Write("data")
      result.Write()
    
      # Get the residual of the hist
      residual = getResidual(self.hist,result,self.binLow,self.binHigh)
      residual.SetDirectory(0)
      residual.SetName("residual_"+valType)
      residual.SetTitle("residual_"+valType)
      residual.Write()

    # Compare to nominal fit result from the old code, if available
    if self.nominalFitFile :
      nominalFit = self.nominalFitFile.Get("basicBkgFrom4ParamFit")
      nominalFit.SetDirectory(0)
      nominalFit.SetName("nominal")
      nominalFit.SetTitle("nominal")
      nominalFit.Write()

    outputFile.Close()
    print "Made file",self.outputFileName

  def parseCommandLineArgs(self) :
  
    parser = argparse.ArgumentParser(description='Specify config and overwrite values if desired.')
    parser.add_argument('--dataset')
    self.args = parser.parse_args()

if __name__ == "__main__":

  fitter = RunFitter()
  fitter.executeFit()

