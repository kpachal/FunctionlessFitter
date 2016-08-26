from art.morisot import Morisot
import ROOT

myPainter = Morisot()

#infile = ROOT.TFile("results/test/outputfile_only1stConstraint.root","READ")
#ext = "_only1stConstraint"
#infile = ROOT.TFile("results/test/outputfile_1stAnd2ndConstraints.root","READ")
#ext = "_1stAnd2ndConstraints"
#infile = ROOT.TFile("results/test/outputfile_3Constraints.root","READ")
#ext = "_3Constraints"
#
#infile = ROOT.TFile("results/test/outputfile_0thOrderConstraint.root","READ")
#ext = "_0thOrderConstraint"
#infile = ROOT.TFile("results/test/outputfile_1stOrderConstraint.root","READ")
#ext = "_1stOrderConstraint"
#infile = ROOT.TFile("results/test/outputfile_2ndOrderConstraint.root","READ")
#ext = "_2ndOrderConstraint"
#infile = ROOT.TFile("results/test/outputfile_3rdOrderConstraint.root","READ")
#ext = "_3rdOrderConstraint"
#infile = ROOT.TFile("results/test/outputfile_4thOrderConstraint.root","READ")
#ext = "_4thOrderConstraint"

infile = ROOT.TFile("results/test/outputfile_TLA.root","READ")
ext = "TLA"

#infile = ROOT.TFile("results/test/outputfile_ICHEP.root","READ")
#ext = "ICHEP"

data = infile.Get("basicData")
bkg = infile.Get("basicBkg")
residual = infile.Get("residual")

firstDer = infile.Get("firstDerivative")
secondDer = infile.Get("secondDerivative")
thirdDer = infile.Get("thirdDerivative")

firstDerNom = infile.Get("firstDer_nominalFit")
secondDerNom = infile.Get("secondDer_nominalFit")
thirdDerNom = infile.Get("thirdDer_nominalFit")

firstDerTF1 = infile.Get("firstDer_fromTF1")
secondDerTF1 = infile.Get("secondDer_fromTF1")
thirdDerTF1 = infile.Get("thirdDer_fromTF1")
fourthDerTF1 = infile.Get("fourthDer_fromTF1")

if "TLA" in ext :
  binLow = data.FindBin(395)
  binHigh = data.FindBin(1252)

else :
  binLow = data.FindBin(1100)
  binHigh = 130

minX = data.GetBinLowEdge(binLow)
maxX = data.GetBinLowEdge(binHigh+1)

myPainter.drawDataAndFitOverSignificanceHist(data,bkg,residual,"m_{jj}","Events","Sign.","plots/figure1"+ext,3640,13,-1,-1,binLow,binHigh,doBumpLimits=False,bumpLow=0,bumpHigh=0,extraLegendLines=[],doLogX=True,doRectangular=False,setYRange=[],writeOnpval = False, pval = -999,doWindowLimits=False,windowLow=0,windowHigh=0)

myPainter.drawSeveralObservedLimits([firstDer],["First derivative"],"plots/firstDerivative"+ext,"mjj","f'",3540,13,minX,maxX,-25,0.0,extraLegendLines = [], doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])

myPainter.drawSeveralObservedLimits([secondDer],["Second derivative"],"plots/secondDerivative"+ext,"mjj","f''",3540,13,minX,maxX,-0.02,0.2,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])

myPainter.drawSeveralObservedLimits([thirdDer],["Third derivative"],"plots/thirdDerivative"+ext,"mjj","f'''",3540,13,minX,maxX,-2E-4,2E-4,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])


myPainter.drawSeveralObservedLimits([firstDerNom],["First derivative"],"plots/firstDerivative_nominalFit"+ext,"mjj","f'",3540,13,minX,maxX,-25,0.0,extraLegendLines = [], doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])

myPainter.drawSeveralObservedLimits([secondDerNom],["Second derivative"],"plots/secondDerivative_nominalFit"+ext,"mjj","f''",3540,13,minX,maxX,-0.02,0.2,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])

myPainter.drawSeveralObservedLimits([thirdDerNom],["Third derivative"],"plots/thirdDerivative_nominalFit"+ext,"mjj","f'''",3540,13,minX,maxX,-2E-4,2E-4,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])

myPainter.drawSeveralObservedLimits([firstDer,firstDerNom],["First derivative","Nominal fit result"],"plots/firstDerivative_compare"+ext,"mjj","f''",3540,13,minX,maxX,-25,0.0,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])
myPainter.drawSeveralObservedLimits([secondDer,secondDerNom],["Second derivative","Nominal fit result"],"plots/secondDerivative_compare"+ext,"mjj","f''",3540,13,minX,maxX,-0.02,0.2,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])
myPainter.drawSeveralObservedLimits([thirdDer,thirdDerNom],["Third derivative","Nominal fit result"],"plots/thirdDerivative_compare"+ext,"mjj","f'''",3540,13,minX,maxX,-2E-4,2E-4,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])

myPainter.drawSeveralObservedLimits([firstDerTF1],["First derivative"],"plots/firstDerivative_fromTF1"+ext,"mjj","f'",3540,13,minX,maxX,-25,0.0,extraLegendLines = [], doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])
myPainter.drawSeveralObservedLimits([secondDerTF1],["Second derivative"],"plots/secondDerivative_fromTF1"+ext,"mjj","f''",3540,13,minX,maxX,-0.02,0.2,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])
myPainter.drawSeveralObservedLimits([thirdDerTF1],["Third derivative"],"plots/thirdDerivative_fromTF1"+ext,"mjj","f'''",3540,13,minX,maxX,-2E-4,2E-4,extraLegendLines = [], doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])
myPainter.drawSeveralObservedLimits([fourthDerTF1],["Fourth derivative"],"plots/fourthDerivative_fromTF1"+ext,"mjj","f''''",3540,13,minX,maxX,-4E-6,4E-6,extraLegendLines=[],doLogY=False,doLogX=True,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[])




