import sys
import os
import ROOT
import scipy

from HistWrapper import WrappedHist
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual

def getNthDerivative(func,degree,x,h) :

  val = 0
  for i in range(degree+1) :
    term = pow(-1.0,i) * scipy.special.binom(degree,i) * func.Eval(x + (degree/2.0 - i)*h)
    if degree > 1:
      val = val + term*100.0
    else :
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

# Make a fitter
myFitter = FunctionlessFitter()

# Get a histogram we want to fit. I'll use dijets EOYE.
infile = ROOT.TFile("samples/dataLikeHistograms_IBLOff.2015.root")
hist = infile.Get("Nominal/mjj_Data_2015_3p57fb")
hist.SetDirectory(0)
infile.Close()

binLow = hist.FindBin(1101)
binHigh = -1

#binLow = 77
#binHigh = 95
#binHigh = 95
#binLow = 110
#binHigh = 120
#binHigh = 52
#binHigh = 130

result = myFitter.fit(hist,binLow,binHigh)
wResult = WrappedHist(result)

# Get the residual of the hist
residual = getResidual(hist,result,binLow,binHigh)

# Plot histogram's first and second derivatives
firstDerivative = wResult.graphFirstDerivatives()
secondDerivative = wResult.graphSecondDerivatives()
thirdDerivative = wResult.graphThirdDerivatives()

# Now bump hunt it



# Compare to nominal fit result from the old code
nominalFitFile = ROOT.TFile("samples/Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","READ")
nominalFit = nominalFitFile.Get("basicBkgFrom4ParamFit")
nominalFit.SetDirectory(0)
wNominal = WrappedHist(nominalFit)
firstDerNom = wNominal.graphFirstDerivatives()
secondDerNom = wNominal.graphSecondDerivatives()
thirdDerNom = wNominal.graphThirdDerivatives()

# What about using the function?
nominalFitTF1 = nominalFitFile.Get("theFitFunction")
nominalFitFile.Close()

# Make derivative histograms
firstDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 1, 55.0)
secondDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 2, 55.0)
thirdDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 3, 55.0)
fourthDerFromTF1 = getNthDerivativeGraphFromFunc(residual, nominalFitTF1, 4, 55.0)

# Write everything to a file

outputFile = ROOT.TFile("outputfile_3Constraints.root","RECREATE")
#outputFile = ROOT.TFile("outputfile_1stAnd2ndConstraints.root","RECREATE")
#outputFile = ROOT.TFile("outputfile_only1stConstraint.root","RECREATE")
outputFile.cd()
hist.Write("basicData")
result.Write("basicBkg")
residual.Write("residual")
firstDerivative.Write("firstDerivative")
secondDerivative.Write("secondDerivative")
firstDerNom.Write("firstDer_nominalFit")
secondDerNom.Write("secondDer_nominalFit")
firstDerFromTF1.Write("firstDer_fromTF1")
secondDerFromTF1.Write("secondDer_fromTF1")
thirdDerFromTF1.Write("thirdDer_fromTF1")
fourthDerFromTF1.Write("fourthDer_fromTF1")
outputFile.Close()
