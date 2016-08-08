import ROOT
import numpy

class WrappedHist() :

  def __init__(self,inputHist, seed = 0, binStructure="central") :

    self.histogram = inputHist
    self.getHistOutermostBinsWithData()
    self.bincontents = []
    self.binxvals = []
    self.randomNumberGenerator = ROOT.TRandom3(seed)
    self.recordBinXVals(binStructure)
    
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

  def recordBinXVals(self,type) :
    self.binxvals = []
    for bin in range(self.histogram.GetNbinsX()+2) :
      if type=="central" :
        self.binxvals.append(self.histogram.GetBinCenter(bin))
      elif type=="exp" :
        print "exp"
      else :
        raise Exception( "Unrecognized bin center definition!" )

  def getSelectedBinInfo(self,rangeLow,rangeHigh) :
    selectedbincontents = []
    selectedbinxvals = []
    selectedbinwidths = []
    for bin in range(rangeLow, rangeHigh+1) :
      selectedbincontents.append(self.histogram.GetBinContent(bin))
      selectedbinxvals.append(self.binxvals[bin])
      selectedbinwidths.append(self.histogram.GetBinWidth(bin))
    return selectedbincontents,selectedbinxvals,selectedbinwidths


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

    print "Making 1st der graph from histogram."
    graph = ROOT.TGraph()
    graph.SetName("firstDerivatives_"+self.histogram.GetName())
    index = -1
    for bin in range(1,self.histogram.GetNbinsX()) :
      index = index+1
      run = self.histogram.GetBinCenter(bin+1)-self.histogram.GetBinCenter(bin)
      #run = self.histogram.GetBinLowEdge(bin+2) - self.histogram.GetBinLowEdge(bin)
      rise = self.histogram.GetBinContent(bin+1)/self.histogram.GetBinWidth(bin+1)-self.histogram.GetBinContent(bin)/self.histogram.GetBinWidth(bin)
      slope = rise/run
      xval = self.histogram.GetBinLowEdge(bin+1)
      graph.SetPoint(index,xval,slope)
      if xval > 1900 and xval < 2100 :
        print "For x value",xval
        print "found first derivative (",self.histogram.GetBinContent(bin+1)/self.histogram.GetBinWidth(bin+1),"-",self.histogram.GetBinContent(bin)/self.histogram.GetBinWidth(bin),")/",run,"=",slope

    self.firstDer = graph
    return self.firstDer

  def graphSecondDerivatives(self) :

    # Bin content = height of function integrated over bin
    # so height of function ~ bin content/bin width

    # Use finite difference formulas to do this.
    # D = (y2 - y1)/((x2 - x1)*(x2 - x0)) - (y1 - y0)/((x1 - x0)*(x2-x0))
    
    print "Making 2nd der graph from histogram."
    graph = ROOT.TGraph()
    graph.SetName("secondDerivatives_"+self.histogram.GetName())
    index = -1
    for bin in range(1,self.histogram.GetNbinsX()-2) :
      index = index+1

      run1 = self.histogram.GetBinCenter(bin+1)-self.histogram.GetBinCenter(bin)
      rise1 = self.histogram.GetBinContent(bin+1)/self.histogram.GetBinWidth(bin+1)-self.histogram.GetBinContent(bin)/self.histogram.GetBinWidth(bin)
      run2 = self.histogram.GetBinCenter(bin+2)-self.histogram.GetBinCenter(bin+1)
      rise2 = self.histogram.GetBinContent(bin+2)/self.histogram.GetBinWidth(bin+2)-self.histogram.GetBinContent(bin+1)/self.histogram.GetBinWidth(bin+1)
      slope1 = rise1*100.0/run1
      slope2 = rise2*100.0/run2
      dSlopedX = (slope2-slope1)/(run1+run2)
      xval = self.histogram.GetBinCenter(bin+1)
      graph.SetPoint(index,xval,dSlopedX)

    self.secondDer = graph
    return self.secondDer
