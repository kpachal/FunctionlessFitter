import sys
import os
import ROOT
import scipy
import argparse
import ConfigParser

from HistWrapper import WrappedHist
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual

class RunSearchPhase :

  def __init__(self) :
  
    # Make a fitter
    self.myFitter = FunctionlessFitter()
    
    # Options which can be set by command line
    self.noDataErr = False
    self.useScaled = False
    self.thresholdPVal_FindSignal = 0.01
    self.thresholdPVal_RemoveSignal = 0.01
    self.minBHMass = -1
    self.maxBHMass = -1
    self.inputFileName = ""
    self.outputFileName = ""
    self.inputHistDir = ""
    self.dataMjjHistoName = ""
  
    self.parseCommandLineArgs()

  def execute(self) :

    # Start counting time
    totaltime = ROOT.TStopwatch()
    totaltime.Start()
      
    # Read configuration file to get
    # basic IO information and fit settings
    self.readConfiguration()
    
    # Access input file and get histogram to fit
    basichist = self.readFile()
    # Make wrapped version which we will use for many tests
    self.theHistogram = WrappedHist(basichist)
    
    # Get range for fit
    if self.minX > self.theHistogram.histogram.GetBinLowEdge(self.theHistogram.lastBinWithData+1)\
          or self.minX < 0 :
      self.firstBinFit = self.theHistogram.firstBinWithData
    else :
      self.firstBinFit = self.theHistogram.histogram.FindBin(self.minX)
    if self.maxX < self.theHistogram.histogram.GetBinLowEdge(self.theHistogram.firstBinWithData)\
          or self.maxX < 0 :
      self.lastBinFit = self.theHistogram.lastBinWithData
    else :
      self.lastBinFit = self.theHistogram.histogram.FindBin(self.maxX)

    # Set options in the fitter according to configuration
    self.myFitter.startValFormat = self.startVals
    self.myFitter.flatStartVal = 100.0
    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}

    # Fit the histogram
    prelim_result = self.myFitter.fit(self.theHistogram,self.firstBinFit,self.lastBinFit)
    prelim_result.SetDirectory(0)
    wResult = WrappedHist(prelim_result)

    # Bump hunt it

    # Get the residual of the hist
    residual = getResidual(basichist,prelim_result,self.firstBinFit,self.lastBinFit)
    residual.SetDirectory(0)

    ## Plot histogram's derivatives
    wResult.graphUpToNthDerivatives(4)
    firstDerivative = wResult.der1
    secondDerivative = wResult.der2
    thirdDerivative = wResult.der3
    fourthDerivative = wResult.der4

    # Write everything to a file
    print "Making file",self.outputFileName
    outputFile = ROOT.TFile(self.outputFileName,"RECREATE")
    outputFile.cd()
    basichist.Write("basicData")
    prelim_result.Write("basicBkg")
    residual.Write("residual")
    firstDerivative.Write("firstDerivative")
    secondDerivative.Write("secondDerivative")
    thirdDerivative.Write("thirdDerivative")
    fourthDerivative.Write("fourthDerivative")
    outputFile.Close()

  def parseCommandLineArgs(self) :
  
    parser = argparse.ArgumentParser(description='Specify config and overwrite values if desired.')
    parser.add_argument("--config",help = "Specify configuration file to use")
    parser.add_argument("--file",help="Input file (will overwrite config value)")
    parser.add_argument("--outputfile",help="Output file (will overwrite config value)")
    args = parser.parse_args()
    self.config = args.config
    if args.file: self.inputFileName = args.file
    if args.outputfile: self.outputFileName = args.outputfile

  def readConfiguration(self) :

    # Default values
    defaults = {
                "minXForFit": -1.0,
                "maxXForFit": -1.0,
                "nPseudoExpFit": 100,
                "startVals": "exp",
                "nPseudoExp": 1000,
                "permitWindow": "True"}

    configReader = ConfigParser.RawConfigParser(defaults)
    configReader.read(self.config)
    
    # Settings which can be overwritten by command line :
    section = "IO"
    if not self.inputFileName :
      self.inputFileName = configReader.get(section, "inputFileName")
    if not self.outputFileName :
      self.outputFileName = configReader.get(section, "outputFileName")
    
    # Settings which must be in the config file
    self.dataHistoName = configReader.get(section, "dataHist")
    
    # Settings which are optional and have default values
    section = "Fitting"
    self.minX = configReader.getfloat(section,"minXForFit")
    self.maxX = configReader.getfloat(section,"maxXForFit")
    self.nPseudoExpFit = configReader.getint(section,"nPseudoExpFit")
    self.startVals = configReader.get(section, "startVals")
    if '"' in self.startVals :
      self.startVals = self.startVals[1:-1]
    section = "General"
    self.nPseudoExpBH = configReader.getint(section, "nPseudoExp")
    self.permitWindow = configReader.getboolean(section, "permitWindow")

  def readFile(self) :
  
    # Get a histogram we want to fit
    print type(self.inputFileName)
    print "Retrieving histogram",self.dataHistoName,"from file",self.inputFileName
    infile = ROOT.TFile.Open(self.inputFileName)
    hist = infile.Get(self.dataHistoName)
    hist.SetDirectory(0)
    infile.Close()
    return hist

# Run
if __name__ == "__main__":
  
  fitter = RunSearchPhase()
  fitter.execute()

