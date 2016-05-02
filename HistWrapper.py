import ROOT
import numpy

class WrappedHist :

  def __init__(self,inputHist, seed = 0) :

    self.histogram = inputHist
    self.getHistOutermostBinsWithData()
    self.bincontents = []
    self.binxvals = []
    self.randomNumberGenerator = ROOT.TRandom3(seed)
    
    return

  def getHistOutermostBinsWithData(self,epsilon = -1) :

    firstBin = 1
    while self.histogram.GetBinContent(firstBin)==0 and firstBin < self.histogram.GetNbinsX() :
      firstBin = firstBin+1
    lastBin = self.histogram.GetNbinsX()
    if epsilon < 0 :
      while self.histogram.GetBinContent(lastBin)==0 and lastBin > 1 : lastBin = lastBin-1
    else :
      while self.histogram.GetBinContent(lastBin) < epsilon and lastBin > 1 : lastBin = lastBin-1
    if firstBin == self.histogram.GetNbinsX() and lastBin == 1 :
      print "No data in histogram! Resetting limits to first and last bin."
      firstBin = 1
      lastBin = self.histogram.GetNbinsX()
      
    self.firstBinWithData = firstBin
    self.lastBinWithData = lastBin
    return

  def poissonFluctuateBinByBin(self) :

    pseudoHist = ROOT.TH1D(self.histogram)
    pseudoHist.SetName(self.histogram.GetName()+"_pseudo")
    pseudoHist.Reset()
    for bin in range(self.histogram.GetNbinsX()+1) :
      effExp = self.histogram.GetBinContent(bin)
      weight = 1.0 #fWeightsHistogram.GetBinContent(bin);
      pseudo = self.randomNumberGenerator.PoissonD(effExp)
      # Protect against future floating point errors
      if (numpy.isclose(weight,1.0)) :
        pseudoHist.SetBinContent(bin,pseudo)
        pseudoHist.SetBinError(bin,sqrt(pseudo))
      else :
        pseudoHist.SetBinContent(bin,pseudo*weight)
        pseudoHist.SetBinError(bin,sqrt(pseudo)*weight)

    return pseudoHist

  def graphFirstDerivatives(self) :

    graph = ROOT.TGraph()
    graph.SetName("firstDerivatives_"+self.histogram.GetName())
    index = -1
    for bin in range(1,self.histogram.GetNbinsX()) :
      index = index+1
      run = self.histogram.GetBinCenter(bin+1)-self.histogram.GetBinCenter(bin)
      rise = self.histogram.GetBinContent(bin+1)-self.histogram.GetBinContent(bin)
      slope = rise/run
      xval = self.histogram.GetBinCenter(bin) + 0.5*run
      graph.SetPoint(index,xval,slope)
    self.firstDer = graph
    return self.firstDer

  def graphSecondDerivatives(self) :

    graph = ROOT.TGraph()
    graph.SetName("secondDerivatives_"+self.histogram.GetName())
    index = -1
    for bin in range(1,self.histogram.GetNbinsX()-2) :
      index = index+1
      run1 = self.histogram.GetBinCenter(bin+1)-self.histogram.GetBinCenter(bin)
      rise1 = self.histogram.GetBinContent(bin+1)-self.histogram.GetBinContent(bin)
      slope1 = rise1/run1
      run2 = self.histogram.GetBinCenter(bin+2)-self.histogram.GetBinCenter(bin+1)
      rise2 = self.histogram.GetBinContent(bin+2)-self.histogram.GetBinContent(bin+1)
      slope2 = rise2/run2
      dSlopedX = (slope2-slope1)/(run1+run2)
      xval = self.histogram.GetBinCenter(bin+1)
      graph.SetPoint(index,xval,dSlopedX)
    self.secondDer = graph
    return self.secondDer
