import sys
import os
import ROOT
import scipy
from decimal import *
import argparse

from HistWrapper import WrappedHist
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

    self.outputFileName = "results/test/outputfile_testUncertainties_EOYE.root"

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

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}

    self.myFitter.flatStartVal = 1.0
    
    self.outputFileName = "results/test/outputfile_testUncertainties_TLA.root"

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
    
    self.outputFileName = "results/test/outputfile_testUncertainties_TLAFull.root"

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

    self.outputFileName = "results/test/outputfile_testUncertainties_ICHEP.root"
  
  def executeFit(self) :

    wInput = WrappedHist(self.hist)

    # Define various start values.
    startVals = "data"
    
    resultList = []
    
    # Compute start values for user based ones
    data = []
    startValHist = self.hist.Clone("startValHist")
    startValHist.SetName("startValHist")
    getStartVals = WrappedHist(startValHist)
    if self.binHigh < 0 :
      fullVals, xvals, widths, w1, w2 = wInput.getSelectedBinInfo(self.binLow,getStartVals.lastBinWithData)
    else :
      fullVals, xvals, widths, w1, w2 = wInput.getSelectedBinInfo(self.binLow,self.binHigh)
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
    result = self.myFitter.fit(wInput,self.binLow,self.binHigh,errType = "Bootstrap")
    result.SetDirectory(0)
    result.SetName("result")
    result.SetTitle("result")

    # Write everything to a file
    result.Write()
    
    # Get the residual of the hist
    residual = getResidual(self.hist,result,self.binLow,self.binHigh)
    residual.SetDirectory(0)
    residual.SetName("residual")
    residual.SetTitle("residual")
    residual.Write()

    # Fetch bin content vectors.
    binErrVecs = self.myFitter.lastUncertaintyBinContents
    
    # What were success rates for uncertainty fits?
    statusVec = self.myFitter.lastUncertaintyFitStatuses
    thisHist = makeHistFromVector(statusVec)
    thisHist.SetName("uncertaintyFitStatuses")
    thisHist.SetTitle("uncertaintyFitStatuses")
    thisHist.Write()

    # Make into histograms comparable to nominal ones.
    binindex = -1
    for binContents in binErrVecs :
      binindex = binindex + 1
      thisHist = makeHistFromVector(binContents)
      thisHist.SetDirectory(0)
      thisHist.SetName("NewBinErr_{0}".format(binindex))
      thisHist.SetTitle("NewBinErr_{0}".format(binindex))
      thisHist.Write()

    # Compare to nominal fit result from the old code, if available
    if self.nominalFitFile :
      nominalFit = self.nominalFitFile.Get("basicBkgFrom4ParamFit")
      nominalFit.SetDirectory(0)
      nominalFit.SetName("nominalWStatErr")
      nominalFit.SetTitle("nominalWFuncChoiceErr")
      nominalFit.Write()
      fitWFuncErr = self.nominalFitFile.Get("nomOnDataWithSymmetricRMSScaleFuncChoiceErr")
      fitWFuncErr.SetDirectory(0)
      fitWFuncErr.SetName("nominalWFuncChoiceErr")
      fitWFuncErr.SetTitle("nominalWFuncChoiceErr")
      fitWFuncErr.Write()
      binindex = -1
      for binContents in binErrVecs :
        binindex = binindex + 1
        toGet = "BinErr_{0}".format(binindex)
        thisHist = self.nominalFitFile.Get(toGet)
        thisHist.SetDirectory(0)
        thisHist.SetName("OldBinErr_{0}".format(binindex))
        thisHist.SetTitle("OldBinErr_{0}".format(binindex))
        thisHist.Write()

    outputFile.Close()
    print "Made file",self.outputFileName

  def parseCommandLineArgs(self) :
  
    parser = argparse.ArgumentParser(description='Specify config and overwrite values if desired.')
    parser.add_argument('--dataset')
    self.args = parser.parse_args()

if __name__ == "__main__":

  fitter = RunFitter()
  fitter.executeFit()

