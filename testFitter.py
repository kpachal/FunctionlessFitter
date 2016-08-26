import sys
import os
import ROOT
import scipy

from HistWrapper import WrappedHist
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual
from BumpHunter import BumpHunter
from Chi2Test import Chi2Test
from LogLikelihoodTest import LogLikelihoodTest
from PseudoExperimenter import PseudoExperimenter

def getNthDerivative(func,degree,x,h) :

  val = 0
  for i in range(degree+1) :
    term = pow(-1.0,i) * scipy.special.binom(degree,i) * func.Eval(x + (degree/2.0 - i)*h)
    val = val + term
  return val/pow(h,degree)

def getNthDerivativeGraphFromFunc(template, func, degree, h) :

  graph = ROOT.TGraph()
  graph.SetName("fromTF1_Der{0}".format(degree))

  index = -1
  for bin in range(1,template.GetNbinsX()+1) :
    index = index+1
    val = template.GetBinCenter(bin)
    if template.GetBinLowEdge(bin) < 1101 :
      graph.SetPoint(index,val,0.0)
      continue
    D = getNthDerivative(func,degree,val,h)
    graph.SetPoint(index,val,D)
  return graph

class RunFitter :

  def __init__(self) :
    # Make a fitter
    self.myFitter = FunctionlessFitter()

    self.setEOYEValues()
    #self.setTLAValues()
    #self.setICHEPValues()

  def setEOYEValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/dataLikeHistograms_IBLOff.2015.root")
    self.nominalFitFile = ROOT.TFile("samples/Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","READ")
    self.hist = self.infile.Get("Nominal/mjj_Data_2015_3p57fb")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(1101)
    self.binHigh = -1

    #self.outputFileName = "results/test/outputfile_0thOrderConstraint.root"
    #self.outputFileName = "results/test/outputfile_1stOrderConstraint.root"
    #self.outputFileName = "results/test/outputfile_2ndOrderConstraint.root"
    self.outputFileName = "results/test/outputfile_3rdOrderConstraint.root"
    #self.outputFileName = "results/test/outputfile_4thOrderConstraint.root"

    #self.myFitter.derivativeConstraints = {0:-1}
    #self.myFitter.derivativeConstraints = {0:-1, 1:1}
    #self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1}
    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}
    #self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1, 4:-1}

    self.myFitter.startValFormat = "flat"
    self.myFitter.flatStartVal = 1.0

  def setTLAValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/mjj_fullDataset.root")
    self.nominalFitFile = ROOT.TFile("samples/SearchResultData_UA2_fullDataset_yStar0p3_from394_permitWindow.root","READ")
    self.hist = self.infile.Get("mjj")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(395)
    self.binHigh = self.hist.FindBin(1100) #1252

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}

    self.myFitter.startValFormat = "flat"
    self.myFitter.flatStartVal = 5E4
    
    self.outputFileName = "results/test/outputfile_TLA.root"

  def setICHEPValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/dataLikeHistograms.2016_DS2p1_Resonance_Fixed.root")
    #self.nominalFitFile = ROOT.TFile("samples/SearchResultData_UA2_fullDataset_yStar0p3_from394_permitWindow.root","READ")
    self.hist = self.infile.Get("Nominal/mjj_Data_2016_15p7fb")
    self.hist.SetDirectory(0)
    self.infile.Close()

    self.binLow = self.hist.FindBin(1100)
    self.binHigh = -1

    self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}
    #self.myFitter.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1, 4:-1}

    self.outputFileName = "results/test/outputfile_ICHEP.root"
  
  def executeFit(self) :

    wInput = WrappedHist(self.hist)
    result = self.myFitter.fit(wInput,self.binLow,self.binHigh)
    result.SetDirectory(0)
    wResult = WrappedHist(result)

    # Get the residual of the hist
    residual = getResidual(self.hist,result,self.binLow,self.binHigh)
    residual.SetDirectory(0)

    ## Plot histogram's derivatives
    wResult.graphUpToNthDerivatives(4)
    firstDerivative = wResult.der1
    secondDerivative = wResult.der2
    thirdDerivative = wResult.der3
    fourthDerivative = wResult.der4

    # Make a bump hunter
    bumpHunter = BumpHunter()
    bumpHunter.minBinsInBump = 2
    bumpHunter.useSidebands = False
    
    # Bump hunt fitted data
    bhStat = bumpHunter.doTest(wInput, wResult, self.binLow,self.binHigh)

    # Compare to nominal fit result from the old code
    nominalFit = self.nominalFitFile.Get("basicBkgFrom4ParamFit")
    nominalFit.SetDirectory(0)
    wNominal = WrappedHist(nominalFit)
    wNominal.graphUpToNthDerivatives(4)
    firstDerNom = wNominal.der1
    secondDerNom = wNominal.der2
    thirdDerNom = wNominal.der3
    fourthDerNom = wNominal.der4

    # Bump hunt with the old fit
    bhStatNom = bumpHunter.doTest(wInput, wNominal, self.binLow,self.binHigh)

    # Do PEs and get an actual p-value.
    # Make some other tests to perform, too.
    PEMaker = PseudoExperimenter()
    mychi2Test = Chi2Test()
    myLogLTest = LogLikelihoodTest()
    PEDict = PEMaker.getPseudoexperiments(wInput,wResult,[bumpHunter,mychi2Test,myLogLTest],self.binLow,self.binHigh,nPEs=50)
    BHPVal = PEDict[0]["pValue"]
    statHist = PEDict[0]["statHist"]
    
    # ... and old fit
    NomPEDict = PEMaker.getPseudoexperiments(wInput,wNominal,[bumpHunter,mychi2Test,myLogLTest],self.binLow,self.binHigh,nPEs=1000)

    # How does it compare to my C++ BumpHunter?
    bhStatNomOfficial = self.nominalFitFile.Get("bumpHunterStatOfFitToData")[0]
    bhPValNomOfficial = self.nominalFitFile.Get("bumpHunterStatOfFitToData")[1]
    print "BH stat from nominal fit, original:",bhStatNomOfficial
    print "BH stat from nominal fit, my BH:",bhStatNom
    print "BH pval from nominal fit, original:",bhPValNomOfficial
    print "BH pval from nominal fit, my PE:",NomPEDict[0]["pValue"]
    print "BH stat from functionless fit:",bhStat
    print "BH pval from functionless fit:",PEDict[0]["pValue"]

    # What about using the function?
    nominalFitTF1 = self.nominalFitFile.Get("theFitFunction")
    
    # Check if my residual calculator is working correctly:
    # do I reproduce the one in the nominal result correctly?
    nominalResidual = self.nominalFitFile.Get("residualHist")
    nominalResidual.SetDirectory(0)
    self.nominalFitFile.Close()
    reproducedResidual = getResidual(self.hist,nominalFit,self.binLow,self.binHigh)

    # Make derivative histograms
    firstDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 1, 55.0)
    secondDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 2, 55.0)
    thirdDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 3, 55.0)
    fourthDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 4, 55.0)

    # Write everything to a file
    print "Making file",self.outputFileName
    outputFile = ROOT.TFile(self.outputFileName,"RECREATE")
    outputFile.cd()
    self.hist.Write("basicData")
    result.Write("basicBkg")
    residual.Write("residual")
    firstDerivative.Write("firstDerivative")
    secondDerivative.Write("secondDerivative")
    thirdDerivative.Write("thirdDerivative")
    fourthDerivative.Write("fourthDerivative")
    firstDerNom.Write("firstDer_nominalFit")
    secondDerNom.Write("secondDer_nominalFit")
    thirdDerNom.Write("thirdDer_nominalFit")
    nominalResidual.Write("residual_nominalFit")
    reproducedResidual.Write("residual_nominalFit_myfunc")
    firstDerFromTF1.Write("firstDer_fromTF1")
    secondDerFromTF1.Write("secondDer_fromTF1")
    thirdDerFromTF1.Write("thirdDer_fromTF1")
    fourthDerFromTF1.Write("fourthDer_fromTF1")
    outputFile.Close()


if __name__ == "__main__":

  fitter = RunFitter()
  fitter.executeFit()

