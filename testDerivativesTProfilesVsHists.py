import sys
import os
import ROOT
import scipy
from sympy import *
from art.morisot import Morisot
import numpy
from array import array

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

    self.derivativeFuncs = {}
    self.lowestVal = 1E8
    self.highestVal = -1E8
    self.extraString = ""

  def setValues(self,string) :
    if string == "analyticalTest" :
      self.setSimpleTestValues()
    elif string == "analyticalTest_variableBins" :
      self.setVariableBinTestValues()
    elif string == "analyticalTest_randomBins" :
      self.setRandomBinTestValues()
    elif string == "Moriond2017" :
      self.setMoriond2017DijetValues()
    elif string == "Moriond2017NoRoot" :
      self.setMoriond2017DijetValuesWithoutRoot()
    elif string == "Moriond2017EqualBins" :
      self.setMoriond2017DijetValuesEqualBins()
    else :
      print "Unrecognized request for start values!"
      return

  def setSimpleTestValues(self) :
  
    print "Setting simple exponential values!"
    name = "analyticalTest"
    
    self.firstVal = 1.0
    self.lastVal = 20.0
    self.binLow = 1
    self.binHigh = 38
    
    # Make my function with easy analytically solvable derivatives.
    # I pick 1000*e^(-x), where the 1e3 exists in order to ensure that my numbers
    # are large enough to get values near 1 easily in higher order derivatives.
    # For bin low edge a and bin high edge b, analytical solution for bin "center" is:
    # center = - ln(1/2) - ln(e^(-b) + e^(-a))
    x = Symbol('x')
    y = 3000.0 * exp(-x)
    self.plotFunction(y, x, self.firstVal, self.lastVal, name)
  
    self.fitTProfile = ROOT.TProfile("exponentialProf","exponentialProf",38,1.0,20.0)
    self.binLow = 1
    self.binHigh = self.fitTProfile.GetNbinsX()
    # Fill the tprofile
    for bin in range (self.binLow,self.binHigh+1) :
      lowEdge = self.fitTProfile.GetBinLowEdge(bin)
      highEdge = self.fitTProfile.GetBinLowEdge(bin+1)
      xval = - numpy.log(0.5) - numpy.log(numpy.exp(- highEdge) + numpy.exp(- lowEdge))
      integral = - 3000.0 * (numpy.exp(- highEdge) - numpy.exp(- lowEdge))
      self.fitTProfile.Fill(xval,xval,integral)
  
    self.hist = self.fitTProfile.ProjectionX("fakeData","B")
    self.nominalFit = self.fitTProfile.ProjectionX("exponentialHist","B")

    self.extraString = ""

  def setVariableBinTestValues(self) :
  
    print "Setting simple exponential values with variable bin widths!"
    name = "analyticalTest_variableBins"
    
    self.firstVal = 1.0
    self.lastVal = 20.0
    bins = [self.firstVal]
    val = self.firstVal
    i = 0.1
    while val < self.lastVal :
      val = val + i
      i = 0.1*(1+len(bins))
      bins.append(val)
    self.binLow = 1
    self.binHigh = len(bins)
    print "Using bins",bins
    
    # Make my function with easy analytically solvable derivatives.
    # I pick 1000*e^(-x), where the 1e3 exists in order to ensure that my numbers
    # are large enough to get values near 1 easily in higher order derivatives.
    # For bin low edge a and bin high edge b, analytical solution for bin "center" is:
    # center = - ln(1/2) - ln(e^(-b) + e^(-a))
    x = Symbol('x')
    y = 3000.0 * exp(-x)
    self.plotFunction(y, x, self.firstVal, self.lastVal, name)
  
    self.fitTProfile = ROOT.TProfile("exponentialProf","exponentialProf",len(bins)-1,array('d',bins))
    self.binLow = 1
    self.binHigh = self.fitTProfile.GetNbinsX()
    # Fill the tprofile
    for bin in range (self.binLow,self.binHigh+1) :
      lowEdge = self.fitTProfile.GetBinLowEdge(bin)
      highEdge = self.fitTProfile.GetBinLowEdge(bin+1)
      xval = - numpy.log(0.5) - numpy.log(numpy.exp(- highEdge) + numpy.exp(- lowEdge))
      integral = - 3000.0 * (numpy.exp(- highEdge) - numpy.exp(- lowEdge))
      self.fitTProfile.Fill(xval,xval,integral)
  
    self.hist = self.fitTProfile.ProjectionX("fakeData","B")
    self.nominalFit = self.fitTProfile.ProjectionX("exponentialHist","B")

    self.extraString = ""


  def setRandomBinTestValues(self) :
  
    print "Setting simple exponential values with variable bin widths!"
    name = "analyticalTest_randomBins"
    
    self.firstVal = 1.0
    self.lastVal = 20.0
    bins = [self.firstVal]
    val = self.firstVal
    i = 0.1
    bin = 0
    while val < self.lastVal :
      bin = bin+1
      val = val + i
      if bin%2==0 :
        i = 0.1*(1+len(bins))
      else :
        i = 0.2
      bins.append(val)
    self.binLow = 1
    self.binHigh = len(bins)
    print "Using bins",bins
    
    # Make my function with easy analytically solvable derivatives.
    # I pick 1000*e^(-x), where the 1e3 exists in order to ensure that my numbers
    # are large enough to get values near 1 easily in higher order derivatives.
    # For bin low edge a and bin high edge b, analytical solution for bin "center" is:
    # center = - ln(1/2) - ln(e^(-b) + e^(-a))
    x = Symbol('x')
    y = 3000.0 * exp(-x)
    self.plotFunction(y, x, self.firstVal, self.lastVal, name)
  
    self.fitTProfile = ROOT.TProfile("exponentialProf","exponentialProf",len(bins)-1,array('d',bins))
    self.binLow = 1
    self.binHigh = self.fitTProfile.GetNbinsX()
    # Fill the tprofile
    for bin in range (self.binLow,self.binHigh+1) :
      lowEdge = self.fitTProfile.GetBinLowEdge(bin)
      highEdge = self.fitTProfile.GetBinLowEdge(bin+1)
      xval = - numpy.log(0.5) - numpy.log(numpy.exp(- highEdge) + numpy.exp(- lowEdge))
      integral = - 3000.0 * (numpy.exp(- highEdge) - numpy.exp(- lowEdge))
      self.fitTProfile.Fill(xval,xval,integral)
  
    self.hist = self.fitTProfile.ProjectionX("fakeData","B")
    self.nominalFit = self.fitTProfile.ProjectionX("exponentialHist","B")

    self.extraString = ""

  def setMoriond2017DijetValues(self) :
  
    print "Setting Moriond 2017 values!"
    name = "Moriond2017"

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/allInputs_4ParGlobalFit_DijetMoriond2017.root","READ")
    self.hist = self.infile.Get("basicData")
    self.hist.SetDirectory(0)
    
    self.nominalFit = self.infile.Get("basicBkgFrom4ParamFit")
    self.nominalFit.SetDirectory(0)

    # Note: we want to do this function analytically
    # so that we can do infinitely precise derivatives
    parvec = self.infile.Get("fittedParameters")
    par0 = parvec[0]
    par1 = parvec[1]
    par2 = parvec[2]
    par3 = parvec[3]

    fitRange = self.infile.Get("FitRange")
    self.firstVal = fitRange[0]
    self.lastVal = fitRange[1]

    self.infile.Close()

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal)
    self.binHigh = self.hist.FindBin(self.lastVal)
  
    x = Symbol('x')
    y = par0*((1-x/13000.0)**par1)*1.0/((x/13000.0)**(par2 + par3*log(x/13000.0)))
    self.plotFunction(y, x, self.firstVal, self.lastVal, name)

    # Now that analytical functions computed, can make a tprofile out of it
    nDivisions = 100000
    func = self.derivativeFuncs[name]["TF1"]
    self.fitTProfile = ROOT.TProfile("fit_profile","fit_profile",self.hist.GetNbinsX(),array('d',self.hist.GetXaxis().GetXbins()))
    for bin in range(self.binLow,self.binHigh+1) :
      print "Bin",bin
      xlow = self.hist.GetBinLowEdge(bin)
      xhigh = self.hist.GetBinLowEdge(bin+1)
      integral = func.Integral(xlow,xhigh)
      print "\tFinding definite integral:"
      print "\t",integral
      for x in range(int(0.35*nDivisions),int(0.5*nDivisions)) :
        xval = float(x*(xhigh-xlow))/float(nDivisions) + xlow
        thisint = func.Integral(xlow,xval)
        if not thisint > 0 :
          xval = self.hist.GetBinCenter(bin)
          break
        if thisint > 0.5*integral :
          break
      print "\t",xval
      self.fitTProfile.Fill(xval,xval,integral)

    self.extraString = "unequalBins"

  def setMoriond2017DijetValuesEqualBins(self) :
  
    print "Setting Moriond 2017 values in an equally-binned histogram!"
    name = "Moriond2017EqualBins"

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/allInputs_4ParGlobalFit_DijetMoriond2017.root","READ")

    # Note: we want to do this function analytically
    # so that we can do infinitely precise derivatives
    parvec = self.infile.Get("fittedParameters")
    par0 = parvec[0]
    par1 = parvec[1]
    par2 = parvec[2]
    par3 = parvec[3]

    fitRange = self.infile.Get("FitRange")
    self.firstVal = fitRange[0]
    self.lastVal = fitRange[1]

    self.infile.Close()

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    
    # Compute bins for this
    self.lowestVal = numpy.ceil(self.lowestVal)
    self.highestVal = numpy.floor(self.highestVal)
    step = 130.0
    nBins = 0
    val = self.lowestVal
    while val < self.highestVal :
      nBins = nBins+1
      val = val+step
    self.highestVal = val
    
    self.binLow = 1
    self.binHigh = nBins

    x = Symbol('x')
    y = par0*((1-x/13000.0)**par1)*1.0/((x/13000.0)**(par2 + par3*log(x/13000.0)))
    self.plotFunction(y, x, self.firstVal, self.lastVal, name)

    # Now that analytical functions computed, can make a tprofile out of it
    nDivisions = 100000
    func = self.derivativeFuncs[name]["TF1"]
    self.fitTProfile = ROOT.TProfile("fit_profile","fit_profile",int(nBins),self.lowestVal,self.highestVal)
    print "Made tprofile with:",
    for bin in range(1,self.fitTProfile.GetNbinsX()+2) : print self.fitTProfile.GetBinLowEdge(bin)," ",
    print ""

    for bin in range(self.binLow,self.binHigh+1) :
      print "Bin",bin
      xlow = self.fitTProfile.GetBinLowEdge(bin)
      xhigh = self.fitTProfile.GetBinLowEdge(bin+1)
      integral = func.Integral(xlow,xhigh)
      print "\tFinding definite integral:"
      print "\t",integral
      for x in range(int(0.35*nDivisions),int(0.5*nDivisions)) :
        xval = float(x*(xhigh-xlow))/float(nDivisions) + xlow
        thisint = func.Integral(xlow,xval)
        if not thisint > 0 :
          xval = self.fitTProfile.GetBinCenter(bin)
          break
        if thisint > 0.5*integral :
          break
      print "\t",xval
      self.fitTProfile.Fill(xval,xval,integral)
    print "Finished filling profile"

    self.hist = self.fitTProfile.ProjectionX("dataLikeHist","B")
    self.nominalFit = self.fitTProfile.ProjectionX("fitHist","B")
    print "Hist is"
    print self.hist

    self.extraString = "equalBins"

  def setMoriond2017DijetValuesWithoutRoot(self) :
  
    print "Setting Moriond 2017 values, without using ROOT.TF1.Integral()!"
    name = "Moriond2017NoRoot"

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/allInputs_4ParGlobalFit_DijetMoriond2017.root","READ")
    self.hist = self.infile.Get("basicData")
    self.hist.SetDirectory(0)
    
    self.nominalFit = self.infile.Get("basicBkgFrom4ParamFit")
    self.nominalFit.SetDirectory(0)

    # Note: we want to do this function analytically
    # so that we can do infinitely precise derivatives
    parvec = self.infile.Get("fittedParameters")
    par0 = parvec[0]
    par1 = parvec[1]
    par2 = parvec[2]
    par3 = parvec[3]

    fitRange = self.infile.Get("FitRange")
    self.firstVal = fitRange[0]
    self.lastVal = fitRange[1]

    self.infile.Close()

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal)
    self.binHigh = self.hist.FindBin(self.lastVal)

    eqString = "par0*((1-x/13000.0)**par1)*1.0/((x/13000.0)**(par2 + par3*log(x/13000.0)))"
    x = Symbol('x')
    myString = "y = "+ eqString
    exec(myString)
    self.plotFunction(y, x, self.firstVal, self.lastVal, name)

    # Now that analytical functions computed, can make a tprofile out of it
    func = self.derivativeFuncs[name]["TF1"]
    self.fitTProfile = ROOT.TProfile("fit_profile","fit_profile",self.hist.GetNbinsX(),array('d',self.hist.GetXaxis().GetXbins()))
    for bin in range(self.binLow,self.binHigh+1) :
      print "Bin",bin
      xlow = self.hist.GetBinLowEdge(bin)
      xhigh = self.hist.GetBinLowEdge(bin+1)
      print "\tFinding definite integral:"
      coreEq = eqString.replace("par0","{0}".format(par0))
      coreEq = coreEq.replace("par1","{0}".format(par1))
      coreEq = coreEq.replace("par2","{0}".format(par2))
      coreEq = coreEq.replace("par3","{0}".format(par3))
      code = "integral, err = scipy.integrate.quad(lambda x: {0}, xlow, xhigh)".format(coreEq)
      exec(code)
      print "\t",integral
      code = """def func(xval) :
        return scipy.integrate.quad(lambda x: {0}, {1}, xval)[0] - 0.5*{2}""".format(coreEq,xlow,integral)
      exec(code)
      xval = scipy.optimize.fsolve(func,1.0)[0]
      print "\t",xval
      self.fitTProfile.Fill(xval,xval,integral)

    self.extraString = "{0}".format("analytical")


  def executeTest(self,name) :

    binCenters = [self.nominalFit.GetBinCenter(i) for i in range(1,self.nominalFit.GetNbinsX()+1)]
    wInput = Dataset(self.nominalFit,binSpecifier=binCenters)

    if self.binHigh < 0 :
      self.binHigh = wInput.lastBinWithData

    # Compare to histogram approximated from fit
    # This section: approximate using histogram and bin centers
    wHistogram = Dataset(self.nominalFit,binSpecifier=binCenters)
    print "About to graph, using self.binLow =",self.binLow,"and self.binHigh =",self.binHigh
    wHistogram.graphUpToNthDerivatives(4,self.binLow,self.binHigh)
    firstDerHist = wHistogram.der1
    secondDerHist = wHistogram.der2
    thirdDerHist = wHistogram.der3
    fourthDerHist = wHistogram.der4

    # This section: approximate using true tprofile
    newTProfile = Dataset(self.fitTProfile)
    newTProfile.graphUpToNthDerivatives(4,self.binLow,self.binHigh)
    firstDerAppx = newTProfile.der1
    secondDerAppx = newTProfile.der2
    thirdDerAppx = newTProfile.der3
    fourthDerAppx = newTProfile.der4

    # Compare 3 lines: analytical derivative, approximation from hist, approximation from tprofile
    toPlotLists = [[self.derivativeFuncs[name]["firstDerivativeFromTF1"],firstDerAppx,firstDerHist],
                   [self.derivativeFuncs[name]["secondDerivativeFromTF1"],secondDerAppx,secondDerHist],
                   [self.derivativeFuncs[name]["thirdDerivativeFromTF1"],thirdDerAppx,thirdDerHist],
                   [self.derivativeFuncs[name]["fourthDerivativeFromTF1"],fourthDerAppx,fourthDerHist]]
    toPlotNames = 4*[["Analytical","TProfile", "Histogram"]]

    toPlotLines = ["First derivative","Second derivative","Third derivative","Fourth derivative"]
    saveNames = ["1stDer","2ndDer","3rdDer","4thDer"]
    saveNames = ["{0}_".format(name)+i for i in saveNames]
    if self.extraString :
      saveNames = [i+"_{0}".format(self.extraString) for i in saveNames]

    # Make plots
    for plotList, nameList, extraLines, saveName in zip(toPlotLists,toPlotNames,toPlotLines,saveNames) :
      lowy = plotList[0].GetMinimum()
      highy = plotList[0].GetMaximum()
      self.myPainter.drawBasicFunction(plotList, self.firstVal, self.lastVal,"m_{jj}","Value","plotting/plots/testGlobalFitBehaviours/funcOnHist_compareMethods_"+saveName,legendlines = nameList,ylow=(0.8*lowy if lowy > 0 else 1.2*lowy),yhigh=(1.2*highy if highy > 0 else 0.8*highy), makeCanvas=True,doLogY=False,doLogX=True,lineColour = ROOT.kCyan+2,doRectangular = False,doLegendHigh=True,extraLegendLines = [extraLines])

    print "Done"

  def plotFunction(self, f, x, vallow, valhigh, key) :

    zerothDer = f
    firstDer = f.diff(x)
    secondDer = firstDer.diff(x)
    thirdDer = secondDer.diff(x)
    fourthDer = thirdDer.diff(x)
    fifthDer = fourthDer.diff(x)
    sixthDer = fifthDer.diff(x)
    seventhDer = sixthDer.diff(x)
    

    if not key in self.derivativeFuncs.keys() :
      self.derivativeFuncs[key] = {}
    
    # Make plots: plotFunction(self, funcString, vallow, valhigh, title)
    self.callPlotter("{0}".format(zerothDer),vallow,valhigh,key,"TF1",0)
    self.callPlotter("{0}".format(firstDer),vallow,valhigh,key,"firstDerivativeFromTF1",1)
    self.callPlotter("{0}".format(secondDer),vallow,valhigh,key,"secondDerivativeFromTF1",2)
    self.callPlotter("{0}".format(thirdDer),vallow,valhigh,key,"thirdDerivativeFromTF1",3)
    self.callPlotter("{0}".format(fourthDer),vallow,valhigh,key,"fourthDerivativeFromTF1",4)
    self.callPlotter("{0}".format(fifthDer),vallow,valhigh,key,"fifthDerivativeFromTF1",5)
    self.callPlotter("{0}".format(sixthDer),vallow,valhigh,key,"sixthDerivativeFromTF1",6)
    self.callPlotter("{0}".format(seventhDer),vallow,valhigh,key,"seventhDerivativeFromTF1",7)

  def callPlotter(self, funcString, vallow, valhigh, key, title, derivativeOrder=0) :
  
    funcString = funcString.replace("**","^")
    myFunc = ROOT.TF1(key+"_"+title+"_f",funcString, vallow, valhigh)
    myFunc.SetName(key+"_"+title)

    self.derivativeFuncs[key][title] = myFunc

    if derivativeOrder < 1 :
      yname = "Value"
    else :
      yname = "f"+"'"*derivativeOrder
    self.myPainter.drawBasicFunction(myFunc, vallow,valhigh,"m_{jj}",yname,"plotting/plots/testGlobalFitBehaviours/"+title+"_"+key,makeCanvas=True,doLogY=False,doLogX=True,lineColour = ROOT.kCyan+2,doRectangular = False,doLegendHigh=True)


if __name__ == "__main__":

  fitter = RunFitter()
#  for result in ["analyticalTest","analyticalTest_variableBins"] :
  for result in ["Moriond2017EqualBins","analyticalTest_randomBins"] : # "Moriond2017", "Moriond2017NoRoot"
    fitter.setValues(result)
    fitter.executeTest(result)

  for order in fitter.derivativeFuncs[fitter.derivativeFuncs.keys()[0]].keys() :
    functions = []
    for result in fitter.derivativeFuncs.keys() :
      functions.append(fitter.derivativeFuncs[result][order])
    fitter.myPainter.drawBasicFunction(functions,fitter.lowestVal, fitter.highestVal,"m_{jj}","Value","plotting/plots/testGlobalFitBehaviours/"+order+"_all",legendlines = fitter.derivativeFuncs.keys(),makeCanvas=True,doLogY=False,doLogX=True,lineColour = ROOT.kCyan+2,doRectangular = False,doLegendHigh=True)




