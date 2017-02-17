import ROOT
import numpy
from decimal import *
getcontext().prec = 28
import MathFunctions
from array import array

class Dataset() :

  def __init__(self, inputData, binSpecifier=None, function=None, seed=0) :

    if "TProfile" not in type(inputData).__name__ and "TH1" not in type(inputData).__name__ :
      raise TypeError("""\tThis is not a TProfile or a TH1.\n\t\tYou gave me a {0}. \n\t\tNo other data formats are accepted.""".format(type(inputProfile).__name__))

    if "TProfile" not in type(inputData).__name__ :
      if not binSpecifier and not function :
        raise TypeError("""\tThis is not a TProfile.\n\t\tYou gave me a {0}.\n\t\tAdd a function or a vector of bin centers to use this.""".format(type(inputProfile).__name__))
      if binSpecifier :
        if len(binSpecifier)!=inputData.GetNbinsX() :
          raise ValueError("""Vector of x-values must be same length
                       as number of bins in histogram!""")

        self.histogram = inputData
        self.binxvals = [Decimal(item) for item in binSpecifier]
        self.getTProfileFromHistAndVector()

      else :
        self.histogram = inputData
        self.getTProfileFromHistAndFunction(function)

    else :
      self.profile = inputData
      self.getHistFromTProfile()
      self.getBinXValsFromTProfile()

    self.getHistOutermostBinsWithData()
    self.randomNumberGenerator = ROOT.TRandom3(seed)
    
    return

  def getHistFromTProfile(self) :
    self.histogram = inputProfile.ProjectionX("_hist","B")
    return

  def getBinXValsFromTProfile(self) :
    self.binxvals = []
    meanHist = inputProfile.ProjectionX("_means")
    for bin in range(1,meanHist.GetNbinsX()+1) :
      val = meanHist.GetBinContent(bin)
      self.binxvals.push_back(Decimal(val))
    return

  def getTProfileFromHistAndVector(self) :
    print "Getting bins from vector"
    name = self.histogram.GetName()+"_profile"
    self.profile = self.makeTProfileFromBins(self.histogram, name)
    for xval, quantity in zip(self.binxvals,[self.histogram.GetBinContent(i) for i in range(1,self.histogram.GetNbinsX()+1)]) :
      self.profile.Fill(xval,xval,quantity)
    return

  def getTProfileFromHistAndFunction(self,function) :
    print "Getting bins from function"
    self.binxvals = []
    name = self.histogram.GetName()+"_profile"
    self.profile = self.makeTProfileFromBins(self.histogram,name)
    for bin in range(1,self.histogram.GetNbinsX()+1) :
      xlow = self.histogram.GetBinLowEdge(bin)
      xhigh = self.histogram.GetBinLowEdge(bin+1)
      integral = function.Integral(xlow,xhigh)
      #print "Func and low and high edges of bin",bin,"is",function.Eval(xlow),function.Eval(xhigh),". Integral",integral
      for x in range(0,100) :
        xval = x*(xhigh-xlow)/100.0 + xlow
        thisint = function.Integral(xlow,xval)
        #print "testing xval",xval,": integral",xlow,xval,"=",thisint,". Compare to",0.5*integral
        if thisint > 0.5*integral :
          break
      #print "for bin",bin,"kept xval",xval
      self.binxvals.append(Decimal(xval))
      self.profile.Fill(xval,xval,self.histogram.GetBinContent(bin))
    return

  def makeTProfileFromBins(self,hist,newname) :
    testBinArray = hist.GetXaxis().GetXbins()
    if len(testBinArray) > 0 :
      profile = ROOT.TProfile(newname,newname,hist.GetNbinsX(),array('d',testBinArray))
    else :
      profile = ROOT.TProfile(newname,newname,hist.GetNbinsX(),hist.GetBinLowEdge(1),hist.GetBinLowEdge(hist.GetNbinsX()+1))
    return profile

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
    
    if binHigh < 0 : binHigh = self.histogram.GetNbinsX()
    if binLow < 0 : binLow = 1
    
    # selectedbinvals = self.binxvals
    # scaleparsby = all 1's for representative result
    scaleparsby = MathFunctions.getFlatVector(binHigh-binLow+1,1.0)
    
    usebincontents,usebinxvals,usebinwidths,usebinedges,useindexLow,useindexHigh = self.getSelectedBinInfo(binLow,binHigh)
    
    dividedDifferenceDatabase, jacobianDatabase = MathFunctions.computeDividedDifferences(degree,usebinxvals,usebinedges,scaleparsby)
    
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
    #print "usebincontents:",usebincontents
    #print "usebinwidths:",usebinwidths
    pars = numpy.divide(usebincontents,usebinwidths)
    #print "Parameters are",pars
    for order in range(1,degree+1) :
      #print "\nBeginning order",order

      thisorderdict = dividedDifferenceDatabase[int(order)]

      for bin in range(0,binHigh - binLow - order) :
        
        # Bins used range from bin to bin+degree inclusive.
        # Even orders use center of middle bin
        if order%2==0 :
          xval = usebinxvals[bin+int(float(order)/2.0)]
        # Odd orders use bin boundary between central two bins
        else :
          xval = usebinedges[bin+int(float(order)/2.0)+1]

        #print "in bin",bin,":"
        #print thisorderdict[bin]
        code = """myfunc = lambda pars : {0}""".format(thisorderdict[bin])
        exec code
        graphs[order].SetPoint(bin,xval,myfunc(pars))
  
        #print "At xval",xval,"assigning value",myfunc(pars)#thisorderdict[bin]

    for order in range(1,degree+1) :
      code = "self.der{0} = graphs[{0}]".format(order)
      exec code



