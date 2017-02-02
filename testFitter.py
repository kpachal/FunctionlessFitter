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

    self.derivativeFuncs = {}
    self.lowestVal = 1E8
    self.highestVal = -1E8

  def setValues(self,string) :
    if string == "TLA" :
      self.setTLAValues()
    elif string == "ICHEP" :
      self.setICHEPValues()
    elif string == "EOYE" :
      self.setEOYEValues()
    elif string == "TEST" :
      self.setTestValues()
    else :
      print "Unrecognized request for start values!"
      return

  def setEOYEValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/dataLikeHistograms_IBLOff.2015.root")
    #self.nominalFitFile = ROOT.TFile("samples/Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","READ")
    self.nominalFitFile = ROOT.TFile("samples/newNominal_EOYE.root","READ")
    self.hist = self.infile.Get("Nominal/mjj_Data_2015_3p57fb")
    self.hist.SetDirectory(0)
    parvec = self.nominalFitFile.Get("fittedParameters")
    par0 = parvec[0]
    par1 = parvec[1]
    par2 = parvec[2]
    self.infile.Close()

    #self.outputFileName = "results/test/outputfile_0thOrderConstraint.root"
    #self.outputFileName = "results/test/outputfile_1stOrderConstraint.root"
    #self.outputFileName = "results/test/outputfile_2ndOrderConstraint.root"
    self.outputFileName = "results/test/outputfile_3rdOrderConstraint.root"

    self.firstVal = 1100
    #self.lastVal = 7052
    self.lastVal = 1500

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal+1)
    self.binHigh = self.hist.FindBin(self.lastVal-1)
    
    x = Symbol('x')
    y = par0*((1-x/13000.0)**par1)*((x/13000.0)**par2)
    self.plotFunction(y, x, self.firstVal, self.lastVal, "EOYE")
    # Four par is:
    # y = par0*((1-x/13000.0)**par1)*1.0/((x/13000.0)**(par2 + par3*log(x/13000.0)))

  def setTLAValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/mjj_fullDataset.root")
    self.nominalFitFile = ROOT.TFile("samples/SearchResultData_UA2_fullDataset_from443.root","READ")
    self.hist = self.infile.Get("mjj")
    self.hist.SetDirectory(0)
    parvec = self.nominalFitFile.Get("fittedParameters")
    par0 = parvec[0]
    par1 = parvec[1]
    par2 = parvec[2]
    par3 = parvec[3]
    self.infile.Close()
    self.outputFileName = "results/test/outputfile_TLA.root"

    self.firstVal = 444
    self.lastVal = 1260

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal)
    self.binHigh = self.hist.FindBin(self.lastVal) #1252

    # UA2
    x = Symbol('x')
    y = par0*((x/13000.0)**(- par1))*exp(-par2*x/13000.0 - par3*((x/13000.0)**2))
    self.plotFunction(y, x, self.firstVal, self.lastVal, "TLA")

  def setICHEPValues(self) :

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/dataLikeHistograms.2016_DS2p1_Resonance_Fixed.root")
    self.nominalFitFile = ROOT.TFile("samples/Step1_SearchPhase_mjj_Data_2016_15p7fb.root","READ")
    self.hist = self.infile.Get("Nominal/mjj_Data_2016_15p7fb")
    self.hist.SetDirectory(0)
    parvec = self.nominalFitFile.Get("fittedParameters")
    par0 = parvec[0]
    par1 = parvec[1]
    par2 = parvec[2]

    self.infile.Close()

    self.firstVal = 1100
    self.lastVal = 7000

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal)
    self.binHigh = -1

    self.outputFileName = "results/test/outputfile_ICHEP.root"
  
    x = Symbol('x')
    y = par0*((1-x/13000.0)**par1)*((x/13000.0)**par2)
    self.plotFunction(y, x, self.firstVal, self.lastVal, "ICHEP")
    # Four par is:
    # y = par0*((1-x/13000.0)**par1)*1.0/((x/13000.0)**(par2 + par3*log(x/13000.0)))

  def setTestValues(self) :
  
    print "Setting test values!"

    # Get a histogram we want to fit. I'll use dijets EOYE.
    self.infile = ROOT.TFile("samples/toyTinyHists.root")
    self.nominalFitFile = ROOT.TFile("samples/toyTinyHists.root","READ")
    self.hist = self.infile.Get("newHist_varyingBins")
    self.hist.SetDirectory(0)
    parvec = self.nominalFitFile.Get("fittedParameters")
    par0 = parvec[0]
    par1 = parvec[1]
    par2 = parvec[2]

    self.infile.Close()

    self.firstVal = 4370
    self.lastVal = 4782.00

    if self.firstVal < self.lowestVal :
      self.lowestVal = self.firstVal
    if self.lastVal > self.highestVal :
      self.highestVal = self.lastVal
    self.binLow = self.hist.FindBin(self.firstVal)
    self.binHigh = -1

    self.outputFileName = "results/test/outputfile_TEST.root"
  
    x = Symbol('x')
    y = par0*((1-x/13000.0)**par1)*((x/13000.0)**par2)
    self.plotFunction(y, x, self.firstVal, self.lastVal, "TEST")
    # Four par is:
    # y = par0*((1-x/13000.0)**par1)*1.0/((x/13000.0)**(par2 + par3*log(x/13000.0)))


  def executeFit(self,name) :

    wInput = Dataset(self.hist)
    
    if self.binHigh < 0 :
      self.binHigh = wInput.lastBinWithData
    
    histForStartVals = self.hist.Clone()
    histForStartVals.SetName("histForStartVals")
    getStartVals = Dataset(histForStartVals)
    fullVals, xvals, widths, edges, w1, w2 = getStartVals.getSelectedBinInfo(self.binLow,self.binHigh)
    index = -1
    startVals = []
    for item in fullVals :
      index = index+1
      startVals.append(item/widths[index])
    self.myFitter.startValFormat = "user"
    self.myFitter.userStartVals = startVals
    
#    result = self.myFitter.fit(wInput,self.binLow,self.binHigh)
#    result.SetDirectory(0)
#    wResult = Dataset(result)

    # Get the residual of the hist
#    residual = getResidual(self.hist,result,self.binLow,self.binHigh)
#    residual.SetDirectory(0)

    ## Plot histogram's derivatives
#    wResult.graphUpToNthDerivatives(4)
#    firstDerivative = wResult.der1
#    secondDerivative = wResult.der2
#    thirdDerivative = wResult.der3
#    fourthDerivative = wResult.der4

#    # Make a bump hunter
#    bumpHunter = BumpHunter()
#    bumpHunter.minBinsInBump = 2
#    bumpHunter.useSidebands = False
#    
#    # Bump hunt fitted data
#    bhStat = bumpHunter.doTest(wInput, wResult, self.binLow,self.binHigh)

    # Compare to nominal fit result from the old code
    #nominalFit = self.nominalFitFile.Get("basicBkgFrom4ParamFit")
    nominalFit = self.nominalFitFile.Get("nominal_simple")
    nominalFit.SetDirectory(0)
    wNominal = Dataset(nominalFit)
    print "About to graph, using self.binLow =",self.binLow,"and self.binHigh =",self.binHigh
    wNominal.graphUpToNthDerivatives(4,self.binLow,self.binHigh)
    firstDerNom = wNominal.der1
    secondDerNom = wNominal.der2
    thirdDerNom = wNominal.der3
    fourthDerNom = wNominal.der4
    
#    fifthDerNom = wNominal.der5
#    sixthDerNom = wNominal.der6
#    seventhDerNom = wNominal.der7

#    # Bump hunt with the old fit
#    bhStatNom = bumpHunter.doTest(wInput, wNominal, self.binLow,self.binHigh)
#
#    # Do PEs and get an actual p-value.
#    # Make some other tests to perform, too.
#    PEMaker = PseudoExperimenter()
#    mychi2Test = Chi2Test()
#    myLogLTest = LogLikelihoodTest()
#    PEDict = PEMaker.getPseudoexperiments(wInput,wResult,[bumpHunter,mychi2Test,myLogLTest],self.binLow,self.binHigh,nPEs=50)
#    BHPVal = PEDict[0]["pValue"]
#    statHist = PEDict[0]["statHist"]
#    
#    # ... and old fit
#    NomPEDict = PEMaker.getPseudoexperiments(wInput,wNominal,[bumpHunter,mychi2Test,myLogLTest],self.binLow,self.binHigh,nPEs=1000)

#    # How does it compare to my C++ BumpHunter?
#    bhStatNomOfficial = self.nominalFitFile.Get("bumpHunterStatOfFitToData")[0]
#    bhPValNomOfficial = self.nominalFitFile.Get("bumpHunterStatOfFitToData")[1]
#    print "BH stat from nominal fit, original:",bhStatNomOfficial
#    print "BH stat from nominal fit, my BH:",bhStatNom
#    print "BH pval from nominal fit, original:",bhPValNomOfficial
#    print "BH pval from nominal fit, my PE:",NomPEDict[0]["pValue"]
#    print "BH stat from functionless fit:",bhStat
#    print "BH pval from functionless fit:",PEDict[0]["pValue"]
    
    # Check if my residual calculator is working correctly:
    # do I reproduce the one in the nominal result correctly?
#    nominalResidual = self.nominalFitFile.Get("residualHist")
#    nominalResidual.SetDirectory(0)
#    self.nominalFitFile.Close()
#    reproducedResidual = getResidual(self.hist,nominalFit,self.binLow,self.binHigh)

    # Make plots overlaying fit derivatives from func on those from histograms
    for thishist, thisname in [[firstDerNom,"firstDerivative"],[secondDerNom,"secondDerivative"],[thirdDerNom,"thirdDerivative"],[fourthDerNom,"fourthDerivative"]] : #,[fifthDerNom,"fifthDerivative"],[sixthDerNom,"sixthDerivative"],[seventhDerNom,"seventhDerivative"]] :
      for index in range(thishist.GetN()) :
        xval = thishist.GetX()[index]
        print "index",index," at x",xval,": comp.",thishist.GetY()[index],",true",self.derivativeFuncs[name][thisname+"FromTF1"].Eval(xval)
      lowy = self.derivativeFuncs[name][thisname+"FromTF1"].GetMinimum()
      highy = self.derivativeFuncs[name][thisname+"FromTF1"].GetMaximum()
      self.myPainter.drawBasicFunction([self.derivativeFuncs[name][thisname+"FromTF1"],thishist], self.firstVal, self.lastVal,"m_{jj}","Value","plotting/plots/testGlobalFitBehaviours/funcOnHist_"+thisname+"_"+name,legendlines = ["analytical","approximate"],ylow=(0.8*lowy if lowy > 0 else 1.2*lowy),yhigh=(1.2*highy if highy > 0 else 0.8*highy), makeCanvas=True,doLogY=False,doLogX=True,lineColour = ROOT.kCyan+2,doRectangular = False)

    # Write everything to a file
    outputFile = ROOT.TFile(self.outputFileName,"RECREATE")
    outputFile.cd()
    self.hist.Write("basicData")
#    result.Write("basicBkg")
#    residual.Write("residual")
#    firstDerivative.Write("firstDerivative")
#    secondDerivative.Write("secondDerivative")
#    thirdDerivative.Write("thirdDerivative")
#    fourthDerivative.Write("fourthDerivative")
    firstDerNom.Write("firstDer_nominalFit")
    secondDerNom.Write("secondDer_nominalFit")
    thirdDerNom.Write("thirdDer_nominalFit")
    fourthDerNom.Write("fourthDer_nominalFit")
#    fifthDerNom.Write("fifthDer_nominalFit")
#    sixthDerNom.Write("sixthDer_nominalFit")
#    seventhDerNom.Write("seventhDer_nominalFit")
#    nominalResidual.Write("residual_nominalFit")
#    reproducedResidual.Write("residual_nominalFit_myfunc")
    for handle, func in self.derivativeFuncs[name].iteritems() :
      func.Write()
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
  for result in ["EOYE"] :
  #for result in ["TEST"] : #,"TLA","ICHEP"] :
    fitter.setValues(result)
    fitter.executeFit(result)

  #for order in fitter.derivativeFuncs["TEST"].keys() :
  for order in fitter.derivativeFuncs["EOYE"].keys() :
    functions = []
    for result in fitter.derivativeFuncs.keys() :
      functions.append(fitter.derivativeFuncs[result][order])
    fitter.myPainter.drawBasicFunction(functions,fitter.lowestVal, fitter.highestVal,"m_{jj}","Value","plotting/plots/testGlobalFitBehaviours/"+order+"_all",legendlines = fitter.derivativeFuncs.keys(),makeCanvas=True,doLogY=False,doLogX=True,lineColour = ROOT.kCyan+2,doRectangular = False)




