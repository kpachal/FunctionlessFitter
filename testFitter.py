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

binLow = hist.FindBin(1101)
binHigh = -1

#binLow = 77
#binHigh = 95
#binHigh = 95
#binLow = 100
#binHigh = 120
#binHigh = 52
binHigh = 130

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


# Write everything to a file

outputFile = ROOT.TFile("outputfile_1stAnd2ndConstraints.root","RECREATE")
outputFile.cd()
hist.Write("basicData")
result.Write("basicBkg")
residual.Write("residual")
firstDerivative.Write("firstDerivative")
secondDerivative.Write("secondDerivative")
outputFile.Close()
