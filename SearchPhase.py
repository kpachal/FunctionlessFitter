import sys
import os
import ROOT
import scipy
import argparse
import ConfigParser
from operator import itemgetter

from HistWrapper import WrappedHist
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual,getRelativeDifference,getSignificanceOfDifference
from BumpHunter import BumpHunter
from Chi2Test import Chi2Test
from LogLikelihoodTest import LogLikelihoodTest
from PseudoExperimenter import PseudoExperimenter
from MathFunctions import makeHistFromVector

class RunSearchPhase :

  def __init__(self) :
  
    # Make a fitter
    self.myFitter = FunctionlessFitter()
    #self.myFitter.minAlg = "COBYLA"
    
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
    self.theHistogram = WrappedHist(basichist)#,scaleBy=1)
    #self.minX = self.minX
    #self.maxX = self.maxX
    
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
    if not self.startVals == "fromHist" :
      self.myFitter.startValFormat = self.startVals
      self.myFitter.flatStartVal = 1.0
    else :
      startVals = []
      self.histForStartVals.SetName("histForStartVals")
      getStartVals = WrappedHist(self.histForStartVals)
      fullVals, xvals, widths, w1, w2 = getStartVals.getSelectedBinInfo(self.firstBinFit,self.lastBinFit)
      index = -1
      for item in fullVals :
        index = index+1
        #startVals.append(item)
        startVals.append(item/widths[index])
      self.myFitter.startValFormat = "user"
      self.myFitter.userStartVals = startVals
    self.myFitter.derivativeConstraints = { 0:-1, 1:1, 2:-1, 3:1, 4:-1}#, 5:1, 6:-1, 7:1}

    # Fit the histogram
    prelim_result = self.myFitter.fit(self.theHistogram,self.firstBinFit,self.lastBinFit)
    prelim_result.SetName("prelim_result")
    wTempResult = WrappedHist(prelim_result)

    # Make a bump hunter
    bumpHunter = BumpHunter()
    bumpHunter.minBinsInBump = 2
    bumpHunter.useSidebands = False
    bumpFound = False
    
    # Use it to get a p-value for the spectrum:
    PEMaker = PseudoExperimenter()
    prelimPEDict = PEMaker.getPseudoexperiments(self.theHistogram,wTempResult,bumpHunter,self.firstBinFit,self.lastBinFit,nPEs=min(self.nPseudoExpBH,100))
    BHPVal = prelimPEDict[0]["pValue"]
    firstBinInWindow = prelimPEDict[0]["furtherInformation"][0]
    lastBinInWindow = prelimPEDict[0]["furtherInformation"][1]
    
    print "Found bumphunter p-value", BHPVal,"and excludewindow is",self.permitWindow
    
    # Exclude window if necessary, and continue expanding it until satisfied
    while BHPVal < 0.05 and self.permitWindow :

      bumpFound = True
      print "Found bump; excluding window!"
      print "At this point firstBinInWindow and lastBinInWindow are",firstBinInWindow,lastBinInWindow,"while tests run between bins",self.firstBinFit,"and",self.lastBinFit

      # Set start values to use last solution so it's more stable
      self.myFitter.startValFormat = "user"
      self.myFitter.userStartVals = self.myFitter.result

      # We don't want to use windows that touch the edge of the
      # fit range. If that happens take the second biggest discrepancy
      # and exclude a window there instead.
      if firstBinInWindow == self.firstBinFit :
      
        # Need to re-BumpHunt data to get rid of result coming from pseudoexperiments
        bumpHunter.doTest(self.theHistogram, wTempResult, self.firstBinFit,self.lastBinFit)
        allStats = bumpHunter.bumpInfoList
        sortedStats = sorted(allStats, key=itemgetter("prob"))
        print "Originally found most discrepant prob.", sortedStats[0]["prob"]
        print "Using the second most extreme instead:", sortedStats[1]["prob"]
        firstBinInWindow = sortedStats[1]["binlow"]
        lastBinInWindow = sortedStats[1]["binhigh"]

      # If window is larger than half the spectrum, we don't want it getting larger
      if (lastBinInWindow - firstBinInWindow+1) > (self.lastBinFit - self.firstBinFit+1)/2.0 : break
      
      # Now re-fit excluding this window.
      self.myFitter.excludeWindow = True
      self.myFitter.firstBinInWindow = firstBinInWindow
      self.myFitter.lastBinInWindow = lastBinInWindow
      
      print "About to do fit"
      prelim_result = self.myFitter.fit(self.theHistogram,self.firstBinFit,self.lastBinFit)
      prelim_result.SetName("intermediate_result")
      print "After fit"
      wTempResult = WrappedHist(prelim_result)

      # Check the result, after excluding matching window from BH
      bumpHunter.excludeWindow = True
      bumpHunter.firstBinToExclude = firstBinInWindow
      bumpHunter.lastBinToExclude = lastBinInWindow
      prelimPEDict = PEMaker.getPseudoexperiments(self.theHistogram,wTempResult,bumpHunter,self.firstBinFit,self.lastBinFit,nPEs=min(self.nPseudoExpBH,100))

      # Update p-value.
      # If it is good enough this will be the (nearly) final window location.
      # Otherwise, adjust window edges and prepare for next iteration.
      BHPVal = prelimPEDict[0]["pValue"]
      print "Window and p-val are now [",firstBinInWindow,lastBinInWindow,"] and",BHPVal,"\n"
      if BHPVal < 0.05 :
      
        # Returning this to old version. just replacing with the new window is silly.
        newFirstBin = prelimPEDict[0]["furtherInformation"][0]
        newLastBin = prelimPEDict[0]["furtherInformation"][1]
        if newLastBin == firstBinInWindow - 1 :
          firstBinInWindow = firstBinInWindow - 1
        elif newFirstBin == lastBinInWindow + 1 :
          lastBinInWindow = lastBinInWindow + 1
        else :
          firstBinInWindow = firstBinInWindow - 1
          lastBinInWindow = lastBinInWindow + 1
            
      # If this is the final window location we want additionally to add
      # one more bin to the left hand side of the exclusion window
      else :
        firstBinInWindow = firstBinInWindow - 1

    # Now have our final window location, if any. Re-fit with full statistics for the background
    # estimate uncertainty. Update the fitter to make sure it uses the
    # final window or do not make unnecessary exclusions.
    if bumpFound :
      self.myFitter.excludeWindow = True
      self.myFitter.firstBinInWindow = firstBinInWindow
      self.myFitter.lastBinInWindow = lastBinInWindow
    else :
      self.myFitter.excludeWindow = False

    # Fit...
    final_result = self.myFitter.fit(self.theHistogram,self.firstBinFit,self.lastBinFit)#,errType="Bootstrap")
    final_result.SetName("final_result")
    wResult = WrappedHist(final_result)

    # Other tools, however, should not exclude this window because we want an estimate of signal.
    bumpHunter.excludeWindow = False
    mychi2Test = Chi2Test()
    mychi2Test.excludeWindow = False
    myLogLTest = LogLikelihoodTest()
    myLogLTest.excludeWindow = False

    # Bump hunt fitted data
    finalPEDict = PEMaker.getPseudoexperiments(self.theHistogram,wResult,[bumpHunter,mychi2Test,myLogLTest],self.firstBinFit,self.lastBinFit,nPEs=self.nPseudoExpBH)
    # Collect statistics
    BHDict = finalPEDict[0]
    Chi2Dict = finalPEDict[1]
    LogLDict = finalPEDict[2]

    # Get the residual of the hist
    residual = getResidual(basichist,final_result,self.firstBinFit,self.lastBinFit)
    residual.SetDirectory(0)

    # Get other significance plots
    relativeDiffHist = getRelativeDifference(basichist,final_result,self.firstBinFit,self.lastBinFit)
    #sigOfDiffHist = getSignificanceOfDifference(basichist,final_result,self.firstBinFit,self.lastBinFit)

    # Evaluate distribution of residuals.
    residualBins = []
    for bin in range(self.firstBinFit, self.lastBinFit+1) :
      residualBins.append(residual.GetBinContent(bin))
    resDistHist = makeHistFromVector(residualBins,1.0)
    resDistHist.Fit("gaus","0")
    fittedGauss = resDistHist.GetFunction("gaus")
    gausmean = fittedGauss.GetParameter(1)
    gauswidth = fittedGauss.GetParameter(2)
    histmean = resDistHist.GetMean()
    histRMS = resDistHist.GetRMS()
    print "Examined distribution of residuals."
    print "Gaussian mean and width were",gausmean,gauswidth

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
    final_result.Write("basicBkg")
    residual.Write("residual")
    resDistHist.Write("residualDistribution")
    fittedGauss.Write("gausFitToResidualDist")
    relativeDiffHist.Write("relativeDiffHist")
    #sigOfDiffHist.Write("sigOfDiffHist")

    statPValErr = ROOT.TVectorD(4)
    BHDict["statHist"].Write("bumpHunterStatHist")
    statPValErr[0] = BHDict["pValue"]
    statPValErr[1] = BHDict["stat"]
    statPValErr[2] = BHDict["furtherInformation"][0]
    statPValErr[3] = BHDict["furtherInformation"][1]
    statPValErr.Write("bumpHunterPValStat")
    BHDict["furtherInformation"][2].Write("bumpHunterTomography")
    statPValErr = ROOT.TVectorD(2)
    Chi2Dict["statHist"].Write("chi2StatHist")
    statPValErr[0] = Chi2Dict["pValue"]
    statPValErr[1] = Chi2Dict["stat"]
    statPValErr.Write("chi2PValStat")
    LogLDict["statHist"].Write("logLikelihoodStatHist")
    statPValErr[0] = LogLDict["pValue"]
    statPValErr[1] = LogLDict["stat"]
    statPValErr.Write("logLPValStat")

    firstDerivative.Write("firstDerivative")
    secondDerivative.Write("secondDerivative")
    thirdDerivative.Write("thirdDerivative")
    fourthDerivative.Write("fourthDerivative")
    outputFile.Close()

    print "Process complete."

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
    self.histNameForStartVals = configReader.get(section,"histForStartVals")
    section = "General"
    self.nPseudoExpBH = configReader.getint(section, "nPseudoExp")
    self.permitWindow = configReader.getboolean(section, "permitWindow")

  def readFile(self) :
  
    # Get a histogram we want to fit
    print "Retrieving histogram",self.dataHistoName,"from file",self.inputFileName
    infile = ROOT.TFile.Open(self.inputFileName)
    hist = infile.Get(self.dataHistoName)
    hist.SetDirectory(0)
    if self.histNameForStartVals :
      self.histForStartVals = infile.Get(self.histNameForStartVals)
      self.histForStartVals.SetDirectory(0)
    infile.Close()
    return hist

# Run
if __name__ == "__main__":
  
  fitter = RunSearchPhase()
  fitter.execute()

