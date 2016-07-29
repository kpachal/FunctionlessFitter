import sys
import os
import ROOT

from HistWrapper import WrappedHist
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual

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

for bin in range(1,result.GetNbinsX()) :
  print "x, data, fit:",hist.GetBinCenter(bin),hist.GetBinContent(bin),result.GetBinContent(bin)

# Get the residual of the hist
residual = getResidual(hist,result,binLow,binHigh)

# Plot histogram's first and second derivatives
firstDerivative = wResult.graphFirstDerivatives()
secondDerivative = wResult.graphSecondDerivatives()

# Now bump hunt it



# Compare to nominal fit result from the old code
nominalFitFile = ROOT.TFile("samples/Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","READ")
nominalFit = nominalFitFile.Get("basicBkgFrom4ParamFit")
nominalFit.SetDirectory(0)
wNominal = WrappedHist(nominalFit)
firstDerNom = wNominal.graphFirstDerivatives()
secondDerNom = wNominal.graphSecondDerivatives()

# What about using the function?
nominalFitTF1 = nominalFitFile.Get("theFitFunction")
nominalFitFile.Close()

# Make derivative histogram
firstDerFromTF1 = ROOT.TGraph()
firstDerFromTF1.SetName("fromTF1_firstDer")
secondDerFromTF1 = ROOT.TGraph()
secondDerFromTF1.SetName("fromTF1_secondDer")

index = -1
nominalFitTF1.Print("all")
print "Examining TF1."
for bin in range(1,residual.GetNbinsX()+1) :
  index = index+1
  val = residual.GetBinCenter(bin)
  h = 55.0
  valp1 = val + h
  valm1 = val - h
  D1 = (nominalFitTF1.Eval(valp1) - nominalFitTF1.Eval(valm1))/(2.0*h)
  D2 = (nominalFitTF1.Eval(valp1) + nominalFitTF1.Eval(valm1) - 2.0*nominalFitTF1.Eval(val))*100/(h*h)
  #if val > 1900 and val < 2100 :
    #print "For x value",val
    #print "found first derivative (",nominalFitTF1.Eval(valp1),"-",nominalFitTF1.Eval(valm1),")/",2.0*(h),"=",D1
    #print "found second derivative",D2
  firstDerFromTF1.SetPoint(index,val,D1)
  secondDerFromTF1.SetPoint(index,val,D2)

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
outputFile.Close()
