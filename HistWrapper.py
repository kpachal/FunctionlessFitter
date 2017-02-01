import ROOT
import numpy
from decimal import *
getcontext().prec = 28
import MathFunctions

class Dataset() :

  def __init__(self,inputHist, seed=0, binStructure="central") :

    self.histogram = inputHist

    self.getHistOutermostBinsWithData()
    self.bincontents = []
    self.binxvals = []
    self.binedges = []
    self.binwidths = []
    self.randomNumberGenerator = ROOT.TRandom3(seed)
    self.recordBinXValsWidths(binStructure)

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

  def recordBinXValsWidths(self,type,slope=-1) :
    self.binxvals = []
    for bin in range(1, self.histogram.GetNbinsX()+1) :
      if type=="central" :
        self.binxvals.append(Decimal(self.histogram.GetBinCenter(bin)))
      elif type=="exp" :
        print "exp"
        xlow = self.histogram.GetBinLowEdge(bin)
        xhigh = self.histogram.GetBinWidth(bin)+xlow
#        ylow =
#        exp =
      elif type=="linear" :
        xlow = self.histogram.GetBinLowEdge(bin)
        binwidth = self.histogram.GetBinWidth(bin)
        if slope < 0 :
          self.binxvals.append(Decimal(xlow+(1.0-1.0/numpy.sqrt(2.0))*binwidth))
      else :
        raise Exception( "Unrecognized bin center definition!" )
      self.binedges.append(Decimal(self.histogram.GetBinLowEdge(bin)))
      self.binwidths.append(Decimal(self.histogram.GetBinWidth(bin)))
      self.bincontents.append(Decimal(self.histogram.GetBinContent(bin)))

  def getSelectedBinInfo(self,rangeLow,rangeHigh,windowLow=-1,windowHigh=-1) :
    selectedbincontents = []
    selectedbinxvals = []
    selectedbinwidths = []
    selectedbinedges = []
    indexLow = -1
    indexHigh = -1
    self.scaleFactors = []
    for bin in range(rangeLow, rangeHigh+1) :
      parVal = self.histogram.GetBinContent(bin)/self.histogram.GetBinWidth(bin)
      scale = 1.0
      while parVal > 10.0 :
        parVal = parVal/10.0
        scale = scale*10
      self.scaleFactors.append(Decimal(scale))

      selectedbincontents.append(Decimal(self.histogram.GetBinContent(bin)))
      selectedbinxvals.append(self.binxvals[bin-1])
      selectedbinwidths.append(Decimal(self.histogram.GetBinWidth(bin)))
      selectedbinedges.append(Decimal(self.histogram.GetBinLowEdge(bin)))
      if bin == windowLow :
        indexLow = len(selectedbincontents)-1
      if bin == windowHigh :
        indexHigh = len(selectedbincontents)-1
    return selectedbincontents,selectedbinxvals,selectedbinwidths,selectedbinedges,indexLow,indexHigh

  def poissonFluctuateBinByBin(self) :

    pseudoHist = ROOT.TH1D(self.histogram)
    pseudoHist.SetName(self.histogram.GetName()+"_pseudo")
    pseudoHist.Reset()
    for bin in range(1,self.histogram.GetNbinsX()+1) :
      effExp = self.histogram.GetBinContent(bin)
      weight = 1.0 #fWeightsHistogram.GetBinContent(bin);
      pseudo = self.randomNumberGenerator.PoissonD(effExp)
      # Protect against future floating point errors
      if (numpy.isclose(weight,1.0)) :
        pseudoHist.SetBinContent(bin,pseudo)
        pseudoHist.SetBinError(bin,numpy.sqrt(pseudo))
      else :
        pseudoHist.SetBinContent(bin,pseudo*weight)
        pseudoHist.SetBinError(bin,numpy.sqrt(pseudo)*weight)

    return pseudoHist

  def graphUpToNthDerivatives(self,degree,binLow=-1,binHigh=-1) :
    
    if binHigh < 0 : binHigh = len(self.binxvals)
    else : binHigh = binHigh - 1
    if binLow < 0 : binLow = 0
    else : binLow = binLow - 1
    
    # selectedbinvals = self.binxvals
    # scaleparsby = all 1's for representative result
    scaleparsby = MathFunctions.getFlatVector(len(self.binxvals),1.0)
    
    dividedDifferenceDatabase, jacobianDatabase = MathFunctions.computeDividedDifferences(degree,self.binxvals,self.binedges,scaleparsby)
    
    graphs = {}
    for order in range(1,degree+1) :
      graph = ROOT.TGraph()
      name = "Derivatives_new_"+self.histogram.GetName()
      if order == 1 :
        name = "1st"+name
      elif order == 2 :
        name = "2nd"+name
      elif order == 3 :
        name = "3rd"+name
      else :
        name = "{0}th".format(order)+name
      graph.SetName(name)
      graphs[order] = graph

    # Parameters are equivalent of bin values divided by binxvals
    pars = numpy.divide(self.bincontents,self.binwidths)
    print "Parameters are",pars
    for order in range(1,degree+1) :
      print "\nBeginning order",order

      thisorderdict = dividedDifferenceDatabase[int(order)]

      #print self.binedges
      for bin in range(binLow,binHigh - order) :
        
        # Bins used range from bin to bin+degree inclusive.
        # Even orders use center of middle bin
        if order%2==0 :
          xval = self.binxvals[bin+int(float(order)/2.0)]
        # Odd orders use bin boundary between central two bins
        else :
          xval = self.binedges[bin+int(float(order)/2.0)+1]

        #print "in bin",bin,":"
        #print thisorderdict[bin]
        code = """myfunc = lambda pars : {0}""".format(thisorderdict[bin])
        exec code
        graphs[order].SetPoint(bin-binLow,xval,myfunc(pars))
  
        print "At xval",xval,"assigning value",thisorderdict[bin]

    for order in range(1,degree+1) :
      code = "self.der{0} = graphs[{0}]".format(order)
      exec code



