import sys
import os
import ROOT
import scipy
import argparse
import ConfigParser
from operator import itemgetter

from HistWrapper import Dataset
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual,getRelativeDifference,getSignificanceOfDifference
from BumpHunter import BumpHunter
from Chi2Test import Chi2Test
from LogLikelihoodTest import LogLikelihoodTest
from PseudoExperimenter import PseudoExperimenter
from MathFunctions import makeHistFromVector

class HGamSideBand :

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

#    ROOT.gROOT.LoadMacro("AtlasStyle.C") 
#    ROOT.SetAtlasStyle()

    # Start counting time
    totaltime = ROOT.TStopwatch()
    totaltime.Start()
      
    # Read configuration file to get
    # basic IO information and fit settings
    self.readConfiguration()
    
    print "Making file",self.outputFileName
    outputFile = ROOT.TFile(self.outputFileName,"RECREATE")

    c = ROOT.TCanvas("c", "c", 800, 600)

    for var in self.varList :
      nBins = len(self.varRanges[var])-1
      print "about to start nBin loop, with nBins",nBins
      
      for b in range(nBins) :
        binName = var + "_bin" + str(b)
        print "in nBin loop"
        
        # Access input file and get histogram to fit
        basicHist = self.readFile(binName)
        if not basicHist :
          continue
        histForStartVals = ROOT.TH1D()

        if self.startVals == "fromHistFit" :
          histForStartVals = self.readFile(binName+"_fit")
        elif self.startVals == "fromHist" :
          histForStartVals = self.readFile(self.histNameForStartVals)

        binVals = []
        for bin in range(1,histForStartVals.GetNbinsX()+1) :
          binVals.append(histForStartVals.GetBinCenter(bin))
      
        # Make wrapped version which we will use for many tests
        self.theHistogram = Dataset(basicHist,binSpecifier = binVals)
        
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
          if self.theHistogram.histogram.GetBinLowEdge(self.lastBinFit) == self.maxX :
            self.lastBinFit = self.lastBinFit - 1

        # Set options in the fitter according to configuration
        print "Setting start vals from config:",self.startVals
        if self.startVals == "fromHistFit" or self.startVals == "fromHist":
          startVals = []
          histForStartVals.SetName("histForStartVals")
          getStartVals = Dataset(histForStartVals,binSpecifier = binVals)
          # basicHist.Print("all")
          # histForStartVals.Print("all")
          # print self.firstBinFit,self.lastBinFit
          fullVals, xvals, widths, edges, w1, w2 = getStartVals.getSelectedBinInfo(self.firstBinFit,self.lastBinFit)
          index = -1
          for item in fullVals :
            index = index+1
            #startVals.append(item)
            startVals.append(item/widths[index])
          if self.startVals == "fromHistFit" :
            self.myFitter.startValFormat = "user"
            self.myFitter.skipSmooth = True
          elif self.startVals == "fromHist" :
            self.myFitter.startValFormat = "user"
          self.myFitter.userStartVals = startVals
        else :
          self.myFitter.startValFormat = self.startVals
          self.myFitter.flatStartVal = 1.0

        # self.myFitter.derivativeConstraints = { 0:-1, 1:1, 2:-1, 3:1, 4:-1}#, 5:1, 6:-1, 7:1}
        self.myFitter.derivativeConstraints = { 0:-1, 1:1, 2:-1, 3:1}#, 4:-1}#, 5:1, 6:-1, 7:1}

        # # Fit the histogram
        # prelim_result = self.myFitter.fit(self.theHistogram,self.firstBinFit,self.lastBinFit)
        # prelim_result.SetName("prelim_result")
        # wTempResult = Dataset(prelim_result,binSpecifier = binVals)

        firstSignalBin = self.theHistogram.histogram.FindBin(self.signalWindowLow)
        lastSignalBin = self.theHistogram.histogram.FindBin(self.signalWindowHigh)-1

        # Exclude the signal region (indexed from 0, not 1 like histogram bins)
        self.myFitter.excludeWindow = True
        self.myFitter.firstBinInWindow = firstSignalBin-1
        self.myFitter.lastBinInWindow = lastSignalBin

        for index in range(firstSignalBin, lastSignalBin+1) :
          # print "Setting data bin",index,"to 0"
          self.theHistogram.histogram.SetBinContent(index, 0)
          self.theHistogram.histogram.SetBinError(index, 0)

        # print "Excluding bins",firstSignalBin,"to",lastSignalBin,"out of",self.theHistogram.histogram.GetNbinsX(),"total bins"

        # self.myFitter.startValFormat = "user"
        # self.myFitter.userStartVals = self.myFitter.result

        # Fit...
        print "About to fit",binName
        final_result = self.myFitter.fit(self.theHistogram,self.firstBinFit,self.lastBinFit)#,errType="Bootstrap")
        print "Done fitting",binName
        final_result.SetName("final_result")
        wResult = Dataset(final_result,binSpecifier = binVals)

        # # Get the residual of the hist
        # residual = getResidual(basichist,final_result,self.firstBinFit,self.lastBinFit)
        # residual.SetDirectory(0)

        # # Get other significance plots
        # relativeDiffHist = getRelativeDifference(basichist,final_result,self.firstBinFit,self.lastBinFit)
        # #sigOfDiffHist = getSignificanceOfDifference(basichist,final_result,self.firstBinFit,self.lastBinFit)

        # # Evaluate distribution of residuals.
        # residualBins = []
        # for bin in range(self.firstBinFit, self.lastBinFit+1) :
        #   residualBins.append(residual.GetBinContent(bin))
        # resDistHist = makeHistFromVector(residualBins,1.0)
        # resDistHist.Fit("gaus","0")
        # fittedGauss = resDistHist.GetFunction("gaus")
        # gausmean = fittedGauss.GetParameter(1)
        # gauswidth = fittedGauss.GetParameter(2)
        # histmean = resDistHist.GetMean()
        # histRMS = resDistHist.GetRMS()
        # print "Examined distribution of residuals."
        # print "Gaussian mean and width were",gausmean,gauswidth

        # ## Plot histogram's derivatives
        # wResult.graphUpToNthDerivatives(4)
        # firstDerivative = wResult.der1
        # secondDerivative = wResult.der2
        # thirdDerivative = wResult.der3
        # fourthDerivative = wResult.der4

        # Write everything to a file
        outputFile.cd()
        
        basicHist.Write(binName)
        if self.startVals == "fromHistFit" :
          histForStartVals.Write(binName+"_fit")
        final_result.Write(binName+"_functionlessFit")

        basicHist.GetXaxis().SetTitle("m_{#gamma#gamma} [GeV]")
        basicHist.GetYaxis().SetTitle("Events")
        basicHist.Draw("pe")

        if self.startVals == "fromHistFit" :
          histForStartVals.SetLineColor(ROOT.kBlue)
          histForStartVals.SetLineStyle(2)
          histForStartVals.Draw("hist same")

        final_result.SetLineColor(ROOT.kRed)
        final_result.Draw("hist same")

        la = ROOT.TLatex()
        la.SetNDC(ROOT.kTRUE)
        la.SetTextSize(0.05)
        la.DrawLatex(0.7, 0.85, binName)

        c.SaveAs("results/dump/"+binName+".png")

        # residual.Write("residual")
        # resDistHist.Write("residualDistribution")
        # fittedGauss.Write("gausFitToResidualDist")
        # relativeDiffHist.Write("relativeDiffHist")
        # sigOfDiffHist.Write("sigOfDiffHist")

        # statPValErr = ROOT.TVectorD(4)
        # BHDict["statHist"].Write("bumpHunterStatHist")
        # statPValErr[0] = BHDict["pValue"]
        # statPValErr[1] = BHDict["stat"]
        # statPValErr[2] = BHDict["furtherInformation"][0]
        # statPValErr[3] = BHDict["furtherInformation"][1]
        # statPValErr.Write("bumpHunterPValStat")
        # BHDict["furtherInformation"][2].Write("bumpHunterTomography")
        # statPValErr = ROOT.TVectorD(2)
        # Chi2Dict["statHist"].Write("chi2StatHist")
        # statPValErr[0] = Chi2Dict["pValue"]
        # statPValErr[1] = Chi2Dict["stat"]
        # statPValErr.Write("chi2PValStat")
        # LogLDict["statHist"].Write("logLikelihoodStatHist")
        # statPValErr[0] = LogLDict["pValue"]
        # statPValErr[1] = LogLDict["stat"]
        # statPValErr.Write("logLPValStat")

        # firstDerivative.Write("firstDerivative")
        # secondDerivative.Write("secondDerivative")
        # thirdDerivative.Write("thirdDerivative")
        # fourthDerivative.Write("fourthDerivative")
        print "Done with",binName
        return
      print "Done with",var

    outputFile.Close()

    print "Process complete."
    totaltime.Stop()
    totaltime.Print()

  def parseCommandLineArgs(self) :
  
    parser = argparse.ArgumentParser(description='Specify config and overwrite values if desired.')
    parser.add_argument("--config",help = "Specify configuration file to use")
    parser.add_argument("--file",help="Input file (will overwrite config value)")
    parser.add_argument("--outputfile",help="Output file (will overwrite config value)")
    parser.add_argument("-b",help="Run in batch mode (no graphics)",action='store_true')
    args = parser.parse_args()
    self.config = args.config
    if args.file: self.inputFileName = args.file
    if args.outputfile: self.outputFileName = args.outputfile

  def readConfiguration(self) :

    # Read configuration options from input TEnv file
    env = ROOT.TEnv(self.config)

    if not self.inputFileName :
      self.inputFileName = env.GetValue("DataHistFile", "DefaultInputFile")
    if not self.outputFileName :
      self.outputFileName = env.GetValue("DataFitFile", "DefaultInputFile")
    self.outputPlotFile = env.GetValue("OutputFitPlots", "DefaultOutputFitPlots")

    # Which differential bins to consider
    variables = env.GetValue("Variables", "DefaultVariables")
    self.varList = variables.split()
    
    binFile = env.GetValue("BinningFile", "DefaultBinningFile")
    binEnv = ROOT.TEnv(binFile)

    self.varRanges = {}
    for var in self.varList :
      binning = binEnv.GetValue(var, "")
      binningStr = binning.split()
      self.varRanges[var] = map(float, binningStr)

    # Fitting options
    self.minX = env.GetValue("FitRangeLow", 105.0)
    self.maxX = env.GetValue("FitRangeHigh", 160.0)

    self.signalWindowLow = env.GetValue("SignalWindowLow", 121.0)
    self.signalWindowHigh = env.GetValue("SignalWindowHigh", 129.0)

    self.startVals = env.GetValue("InitialValues", "DefaultInitialValues")
    if (self.startVals == "fromHist") :
      self.histNameForStartVals = env.GetValue("InitialValuesHist", "DefaultInitialValuesHist")

    self.nPseudoExpFit = env.GetValue("nPseudoExpFit", 1)

  def readFile(self, name) :
  
    # Get a histogram we want to fit
    print "Retrieving histogram",name,"from file",self.inputFileName
    infile = ROOT.TFile.Open(self.inputFileName, "READ")
    hist = infile.Get(name)
    if not hist :
      return None
    hist.SetDirectory(0)
    infile.Close()
    return hist

# Run
if __name__ == "__main__":
  
  fitter = HGamSideBand()
  fitter.execute()

