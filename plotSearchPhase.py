from art.morisot import Morisot
import ROOT

myPainter = Morisot()

#infile = ROOT.TFile("/Users/kpachal/Code/FunctionlessFitter/results/SearchPhase/SearchPhase_EPS2016_TLA_allrange.root","READ")
#ext = "TLA_full"

infile = ROOT.TFile("/Users/kpachal/Code/FunctionlessFitter/results/SearchPhase/SearchPhase_EPS2016_TLA.root","READ")
ext = "TLA"

#infile = ROOT.TFile("results/test/outputfile_ICHEP.root","READ")
#ext = "ICHEP"

data = infile.Get("basicData")
bkg = infile.Get("basicBkg")
residual = infile.Get("residual")

if "TLA" in ext :
  binLow = data.FindBin(395)
  if "full" in ext :
    binHigh = data.FindBin(7000)
  else :
    binHigh = data.FindBin(1252)

minX = data.GetBinLowEdge(binLow)
maxX = data.GetBinLowEdge(binHigh+1)

myPainter.drawDataAndFitOverSignificanceHist(data,bkg,residual,"m_{jj}","Events","Sign.","plots/figure1"+ext,3640,13,binLow,binHigh,doBumpLimits=False,bumpLow=0,bumpHigh=0,extraLegendLines=[],doLogX=True,doRectangular=False,setYRange=[])

