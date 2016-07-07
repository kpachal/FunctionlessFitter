from art.morisot import Morisot
import ROOT

myPainter = Morisot()

#infile = ROOT.TFile("outputfile_only1stConstraint.root","READ")
#ext = "_only1stConstraint"
#infile = ROOT.TFile("outputfile_1stAnd2ndConstraints.root","READ")
#ext = "_1stAnd2ndConstraints"
infile = ROOT.TFile("outputfile_3Constraints.root","READ")
ext = "_3Constraints"

data = infile.Get("basicData")
bkg = infile.Get("basicBkg")
residual = infile.Get("residual")

firstDer = infile.Get("firstDerivative")
secondDer = infile.Get("secondDerivative")

binLow = data.FindBin(1100)
binHigh = 130

minX = data.GetBinLowEdge(binLow)
maxX = data.GetBinLowEdge(binHigh+1)

myPainter.drawDataAndFitOverSignificanceHist(data,bkg,residual,"m_{jj}","Events","Sign.","plots/figure1"+ext,3640,13,-1,-1,binLow,binHigh,doBumpLimits=False,bumpLow=0,bumpHigh=0,extraLegendLines=[],doLogX=True,doRectangular=False,setYRange=[],writeOnpval = False, pval = -999,doWindowLimits=False,windowLow=0,windowHigh=0)

myPainter.drawBasicHistogram(bkg,100,129,"m_{jj}","Events","plots/figure1_zoom"+ext,True,True,True,yLow=0.001,yHigh=12.0)

myPainter.drawSeveralObservedLimits([firstDer],["First derivative"],"plots/firstDerivative"+ext,"mjj","f'",3540,13,minX,maxX,-40,0.0,extraLegendLines = [], doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])

myPainter.drawSeveralObservedLimits([secondDer],["Second derivative"],"plots/secondDerivative"+ext,"mjj","f''",3540,13,minX,maxX,-1.0,3.0,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])
