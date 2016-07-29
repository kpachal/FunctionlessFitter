import sys
import os
import ROOT

from HistWrapper import WrappedHist
from FunctionlessFitter import FunctionlessFitter
from SignificanceTests import getResidual

# Make a fitter
myFitter = FunctionlessFitter()

# Get a histogram we want to fit. I'll use dijets EOYE.
infile = ROOT.TFile("/Users/kpachal/RootFiles/TLA/Data/mjj_fullDataset.root")
hist = infile.Get("mjj")

binLow = hist.FindBin(444)
binHigh = -1

#binLow = 77
#binHigh = 95
#binHigh = 95
#binLow = 100
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


# Write everything to a file

outputFile = ROOT.TFile("outputfile_TLA.root","RECREATE")
#outputFile = ROOT.TFile("outputfile_only1stConstraint.root","RECREATE")
outputFile.cd()
hist.Write("basicData")
result.Write("basicBkg")
residual.Write("residual")
firstDerivative.Write("firstDerivative")
secondDerivative.Write("secondDerivative")
outputFile.Close()
