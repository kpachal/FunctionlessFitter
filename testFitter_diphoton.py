import sys
import os
import ROOT
import scipy
from sympy import *
from art.morisot import Morisot
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
    self.myPainter = Morisot()

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

  def setHighStatValues(self) :

    self.infile = ROOT.TFile("samples/diphoton/data1516_13TeV.xsec_hists.root")
    self.hist = self.infile.Get("pT_yy_bin1")
    self.hist.SetDirectory(0)
    self.fitHist = self.infile.Get("pT_yy_bin1_fit")
    self.fitHist.SetDirectory(0)
    self.infile.Close()
    
    self.binVals = []
    for bin in range(1,self.hist.GetNbinsX()+1) :
      self.binVals.append(self.hist.GetBinCenter(bin))

    self.outputFileName = "results/diphoton/outputfile_test_highstats.root"

    self.firstVal = -1
    self.lastVal = -1

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal+1)
    self.binHigh = self.hist.FindBin(self.lastVal-1)

  def setLowStatValues(self) :

    self.infile = ROOT.TFile("samples/diphoton/data1516_13TeV.xsec_hists.root")
    self.hist = self.infile.Get("pT_yy_bin2")
    print "Getting hist","pT_yy_bin2"
    self.hist.SetDirectory(0)
    self.fitHist = self.infile.Get("pT_yy_bin2_fit")
    self.fitHist.SetDirectory(0)
    self.infile.Close()

    self.binVals = []
    for bin in range(1,self.hist.GetNbinsX()+1) :
      self.binVals.append(self.hist.GetBinCenter(bin))
    
    self.outputFileName = "results/diphoton/outputfile_test_lowstats.root"

    # -1, -1 for full test
    self.firstVal = -1 #130
    self.lastVal = -1 # 150

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal+1)
    self.binHigh = self.hist.FindBin(self.lastVal-1)

  def executeFit(self,name) :

    # Fit with basic function to get weighted bin vals
    fitfunc = ROOT.TF1("nom","[0]/(x^[1])",self.hist.GetBinLowEdge(self.hist.FindBin(self.firstVal)),self.hist.GetBinLowEdge(self.hist.FindBin(self.lastVal)))
    self.hist.Fit("nom")
    fitfunc.Print("all")

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
    self.hist.Write("basicData")
    result.Write("basicBkg")
    self.fitHist.Write("otherBkg")

    #for handle, func in self.derivativeFuncs[name].iteritems() :
      #func.Write()
    outputFile.Close()

  def plotFunction(self, f, x, vallow, valhigh, title) :

    firstDer = f.diff(x)
    secondDer = firstDer.diff(x)
    thirdDer = secondDer.diff(x)
    fourthDer = thirdDer.diff(x)
    fifthDer = fourthDer.diff(x)
    sixthDer = fifthDer.diff(x)
    seventhDer = sixthDer.diff(x)
    

    if not title in self.derivativeFuncs.keys() :
      self.derivativeFuncs[title] = {}
    
    # Make plots: plotFunction(self, funcString, vallow, valhigh, title)
    self.callPlotter("{0}".format(firstDer),vallow,valhigh,"firstDerivativeFromTF1_"+title,1)
    self.callPlotter("{0}".format(secondDer),vallow,valhigh,"secondDerivativeFromTF1_"+title,2)
    self.callPlotter("{0}".format(thirdDer),vallow,valhigh,"thirdDerivativeFromTF1_"+title,3)
    self.callPlotter("{0}".format(fourthDer),vallow,valhigh,"fourthDerivativeFromTF1_"+title,4)
    self.callPlotter("{0}".format(fifthDer),vallow,valhigh,"fifthDerivativeFromTF1_"+title,5)
    self.callPlotter("{0}".format(sixthDer),vallow,valhigh,"sixthDerivativeFromTF1_"+title,6)
    self.callPlotter("{0}".format(seventhDer),vallow,valhigh,"seventhDerivativeFromTF1_"+title,7)

  def callPlotter(self, funcString, vallow, valhigh, title, derivativeOrder=0) :
  
    funcString = funcString.replace("**","^")
    myFunc = ROOT.TF1(title+"_f",funcString, vallow, valhigh)
    myFunc.SetName(title)

    tokens = title.split("_")
    self.derivativeFuncs[tokens[-1]][tokens[0]] = myFunc

    if derivativeOrder < 1 :
      yname = "Value"
    else :
      yname = "f"+"'"*derivativeOrder
    self.myPainter.drawBasicFunction(myFunc, vallow,valhigh,"m_{jj}",yname,"plotting/plots/testGlobalFitBehaviours/"+title,makeCanvas=True,doLogY=False,doLogX=True,lineColour = ROOT.kCyan+2,doRectangular = False)


if __name__ == "__main__":

  fitter = RunFitter()
  for result in ["LowStat"] : #"HighStat",
    fitter.setValues(result)
    fitter.executeFit(result)

#  for order in fitter.derivativeFuncs["EOYE"].keys() :
#    functions = []
#    for result in fitter.derivativeFuncs.keys() :
#      functions.append(fitter.derivativeFuncs[result][order])
#    fitter.myPainter.drawBasicFunction(functions,fitter.lowestVal, fitter.highestVal,"m_{jj}","Value","plotting/plots/testGlobalFitBehaviours/"+order+"_all",legendlines = fitter.derivativeFuncs.keys(),makeCanvas=True,doLogY=False,doLogX=True,lineColour = ROOT.kCyan+2,doRectangular = False)




