import ROOT
import numpy

class WrappedHist() :

  def __init__(self,inputHist, seed=0, binStructure="central",scaleXAxis=False,scaleBy=-1) :

    self.original = inputHist
    self.histogram = inputHist
    self.scaleamount = 1
#    if scaleXAxis :
#      if scaleBy > 0 :
#        self.scaleamount = scaleBy
#      else :
#        self.scaleamount = 1
#        useWidth = self.histogram.GetBinWidth(self.histogram.GetNbinsX())
#        while useWidth > 10.0 :
#          useWidth = useWidth/10.0
#          self.scaleamount = self.scaleamount*10
#      a = self.histogram.GetXaxis()
#      if (a.GetXbins().GetSize()) :
#        oldbins = a.GetXbins()
#        newbins = []
#        for edge in oldbins : newbins.append(edge/self.scaleamount)
#        self.histogram = ROOT.TH1D(self.original.GetName()+"_rebinned",self.original.GetName()+"_rebinned",(len(newbins) - 1), numpy.array(newbins))
#      else :
#        self.histogram = ROOT.TH1D(self.original.GetName()+"_rebinned",self.original.GetName()+"_rebinned",a.GetNbins(), a.GetXmin()/self.scaleamount, a.GetXmax()/self.scaleamount )
      for bin in range(self.histogram.GetNbinsX()+2) :
        self.histogram.SetBinContent(bin,self.original.GetBinContent(bin))
        self.histogram.SetBinError(bin,self.original.GetBinError(bin))

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

  def getSelectedBinInfo(self,rangeLow,rangeHigh,windowLow=-1,windowHigh=-1) :
    selectedbincontents = []
    selectedbinxvals = []
    selectedbinwidths = []
    indexLow = -1
    indexHigh = -1
    for bin in range(rangeLow, rangeHigh+1) :
      selectedbincontents.append(self.histogram.GetBinContent(bin))
      selectedbinxvals.append(self.binxvals[bin])
      selectedbinwidths.append(self.histogram.GetBinWidth(bin))
      if bin == windowLow :
        indexLow = len(selectedbincontents)-1
      if bin == windowHigh :
        indexHigh = len(selectedbincontents)-1
    return selectedbincontents,selectedbinxvals,selectedbinwidths,indexLow,indexHigh

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
        pseudoHist.SetBinError(bin,numpy.sqrt(pseudo))
      else :
        pseudoHist.SetBinContent(bin,pseudo*weight)
        pseudoHist.SetBinError(bin,numpy.sqrt(pseudo)*weight)

    return pseudoHist

  def graphUpToNthDerivatives(self,degree) :
    
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

    dividedDifferenceDatabase = {}
    
    # 0th degree
    baseDict = {}
    for bin in range(len(self.binxvals)) :
      baseDict[bin] = self.histogram.GetBinContent(bin)/self.histogram.GetBinWidth(bin)
    dividedDifferenceDatabase[0] = baseDict

    for order in range(1,degree+1) :

      lastorderdict = dividedDifferenceDatabase[int(order-1)]
      thisorderdict = {}

      for bin in range(self.histogram.GetNbinsX()+2 - order) :
        fxa = lastorderdict[bin]
        fxb = lastorderdict[bin+1]
        diff = (fxb - fxa)/(self.binxvals[bin+order]-self.binxvals[bin])
        
        # Bins used range from bin to bin+degree inclusive.
        # Even orders use center of middle bin
        if order%2==0 :
          xval = self.binxvals[bin+int(float(order)/2.0)]
        # Odd orders use bin boundary between central two bins
        else :
          xval = self.histogram.GetBinLowEdge(bin+int(numpy.ceil(float(order)/2.0)))

        thisorderdict[bin] = diff
        graphs[order].SetPoint(bin,xval,diff)

      dividedDifferenceDatabase[int(order)] = thisorderdict

    for order in range(1,degree+1) :
      code = "self.der{0} = graphs[{0}]".format(order)
      exec code



