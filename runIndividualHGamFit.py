import sys
import os
import ROOT
import scipy
import argparse
from sympy import *
import numpy

import MathFunctions
from HistWrapper import Dataset
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

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}
    self.myFitter.skipSmooth = True

    self.derivativeFuncs = {}
    self.lowestVal = 1E8
    self.highestVal = -1E8

  def setValues(self,string) :
    if string == "HighStat" :
      self.setHighStatValues()
    elif string == "LowStat" :
      self.setLowStatValues()
    else :
      print "Unrecognized request for start values!"
      return

  def setValues(self,histName) :

    self.infile = ROOT.TFile("samples/diphoton/data1516_13TeV.xsec_hists.root")
    self.hist = self.infile.Get(histName)
    self.hist.SetDirectory(0)
    self.fitHist = self.infile.Get(histName+"_fit")
    self.fitHist.SetDirectory(0)
    self.infile.Close()
    
    self.binVals = []
    for bin in range(1,self.hist.GetNbinsX()+1) :
      self.binVals.append(self.hist.GetBinCenter(bin))

    self.outputFileName = "results/diphoton/outputfile_{0}.root".format(histName)

    self.firstVal = -1
    self.lastVal = -1

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal+1)
    self.binHigh = self.hist.FindBin(self.lastVal-1)

  def executeFit(self) :

    # Fit with basic function to get weighted bin vals
#    fitfunc = ROOT.TF1("nom","[0]/(x^[1])",self.hist.GetBinLowEdge(self.hist.FindBin(self.firstVal)),self.hist.GetBinLowEdge(self.hist.FindBin(self.lastVal)))
#    self.hist.Fit("nom")
#    fitfunc.Print("all")

    wInput = Dataset(self.hist,binSpecifier=self.binVals)
#    wInput = Dataset(self.hist,function = fitfunc)
    
    if self.binHigh < 1 :
      self.binHigh = wInput.lastBinWithData
    if self.binLow < 1 :
      self.binLow = wInput.firstBinWithData
    
    histForStartVals = self.fitHist
    getStartVals = Dataset(histForStartVals,binSpecifier=self.binVals)
#    getStartVals = Dataset(histForStartVals,function = fitfunc)
    print "Retrieving values for bins [",self.binLow,self.binHigh,"]"
    print "Thus selected bin xvals should have length",self.binHigh - self.binLow + 1
    fullVals, xvals, widths, edges, w1, w2 = getStartVals.getSelectedBinInfo(self.binLow,self.binHigh)
    index = -1
    startVals = []
    for item in fullVals :
      index = index+1
      startVals.append(item/widths[index])
    self.myFitter.startValFormat = "user"
    self.myFitter.userStartVals = startVals

    self.signalWindowLow = 121.0
    self.signalWindowHigh = 129.0
    
    firstSignalBin = self.hist.FindBin(self.signalWindowLow)
    lastSignalBin = self.hist.FindBin(self.signalWindowHigh)-1
    
    self.myFitter.excludeWindow = True
    self.myFitter.firstBinInWindow = firstSignalBin-1
    self.myFitter.lastBinInWindow = lastSignalBin
    
    # Actually fit it
    result = self.myFitter.fit(wInput,self.binLow,self.binHigh)
    result.SetDirectory(0)
    wResult = Dataset(result,binSpecifier=self.binVals)
#    wResult = Dataset(result,function = fitfunc)

    # Write everything to a file
    outputFile = ROOT.TFile(self.outputFileName,"RECREATE")
    outputFile.cd()
    self.hist.Write("data")
    result.Write("bkgFunctionlessFit")
    self.fitHist.Write("bkgFromChrisFit")
    
    # Record time taken
    fitTime = self.myFitter.lastFitTime
    print "From main body: fit time was",fitTime
    vec = ROOT.TVectorD(1)
    vec[0] = fitTime
    vec.Write("fitCPUTime")
    
    # Get derivatives of functionless fit, up to 3
    wResult.graphUpToNthDerivatives(4)
    firstDerivative = wResult.der1
    secondDerivative = wResult.der2
    thirdDerivative = wResult.der3
    fourthDerivative = wResult.der4
    firstDerivative.Write("ffFirstDerivative")
    secondDerivative.Write("ffSecondDerivative")
    thirdDerivative.Write("ffThirdDerivative")
    fourthDerivative.Write("ffFourthDerivative")

    outputFile.Close()


if __name__ == "__main__":

  fitter = RunFitter()

  parser = argparse.ArgumentParser(description='Specify config and overwrite values if desired.')
  parser.add_argument("--hist",help="Histogram name to use")
  parser.add_argument("-b",help="Run in batch mode (no graphics)",action='store_true')
  args = parser.parse_args()
  histName = args.hist

  # Use to target a specific histogram
  fitter.setValues(histName)
  fitter.executeFit()





