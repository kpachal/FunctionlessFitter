import ROOT
import numpy
from decimal import *
getcontext().prec = 28
import MathFunctions

class Dataset() :

  def __init__(self,inputProfile, seed=0, skipSetup=False) :
    "Initialise Dataset from a ROOT TProfile"

    # I know type checking is poor Python, but people *are* going
    # to feed this histograms, and we simply need to be prepared for that.
    if not type(inputProfile).__name__=="TProfile" :
      raise TypeError("""\tThis is not a TProfile.\n\t\tYou gave me a {0}.\n\t\tUse another constructor or a different input.""".format(type(inputProfile).__name__))

    if skipSetup :
      self.profile = inputProfile
      self.getHistFromTProfile()
      self.getBinXValsFromTProfile()
      
    self.initGeneral()
    return

  @classmethod
  def fromHistAndBinList(cls, inputHist, binSpecifier, seed=0) :
    "Initialise Dataset from an input hist and a list of bin x-values"

    if len(binSpecifier)!=inputHist.GetNbinsX() :
      raise ValueError("""Vector of x-values must be same length
                       as number of bins in histogram!""")

    self.histogram = inputHist
    self.binxvals = binSpecifier
    self.getTProfileFromHistAndVector()
    return cls(self.profile,seed=0,skipSetup=True)
    
  @classmethod
  def fromHistAndFunc(cls,inputHist,function,seed=0) :
    "Initialise Dataset from an input hist and a function fitted to it"

    self.histogram = inputHist
    self.getTProfileFromHistAndFunction(inputHist,function)
    return cls(self.profile,seed=0,skipSetup=True)

  def initGeneral(self) :
  
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
      self.binxvals.push_back(val)
    return

  def getTProfileFromHistAndVector(self) :
    name = self.histogram.GetName()+"_profile"
    self.profile = ROOT.TProfile(name, name, self.histogram.GetNbinsX(), self.histogram.GetXbins().GetArray())
    for xval, quantity in zip(self.binxvals,[self.histogram.GetBinContent(i) for i in range(1,self.histogram.GetNbinsX()+1)]) :
      self.profile.Fill(xval,xval,quantity)
    return

  def getTProfileFromHistAndFunction(self,function) :
    self.binxvals = []
    name = self.histogram.GetName()+"_profile"
    self.profile = ROOT.TProfile(name, name, self.historam.GetNbinsX(), self.histogram.GetXbins().GetArray())
    for bin in range(self.histogram.GetNbinsX()+1) :
      xlow = self.histogram.GetBinLowEdge(bin)
      xhigh = self.histogram.GetBinLowEdge(bin+1)
      integral = function.Integral(xlow,xhigh)
      for x in range(0,100) :
        xval = x*(xhigh-xlow)/100.0 + xlow
        thisint = function.Integral(xlow,xval)
        if thisint > 0.5*integral :
          break
      self.binxvals.append(xval)
      self.profile.Fill(xval,xval,self.histogram.GetBinContent(bin))
    return

#  def estimateBinXValsFromHist(self,type,slope=-1) :
#    self.binxvals = []
#    for bin in range(1, self.histogram.GetNbinsX()+1) :
#      if type=="central" :
#        self.binxvals.append(Decimal(self.histogram.GetBinCenter(bin)))
#      elif type=="exp" :
#        print "exp"
#        xlow = self.histogram.GetBinLowEdge(bin)
#        xhigh = self.histogram.GetBinWidth(bin)+xlow
#        ylow =
#        exp =
#      elif type=="linear" :
#        xlow = self.histogram.GetBinLowEdge(bin)
#        binwidth = self.histogram.GetBinWidth(bin)
#        if slope < 0 :
#          self.binxvals.append(Decimal(xlow+(1.0-1.0/numpy.sqrt(2.0))*binwidth))
#      else :
#        raise Exception( "Unrecognized bin center definition!" )

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
    print "usebincontents:",usebincontents
    print "usebinwidths:",usebinwidths
    pars = numpy.divide(usebincontents,usebinwidths)
    print "Parameters are",pars
    for order in range(1,degree+1) :
      print "\nBeginning order",order

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
  
        print "At xval",xval,"assigning value",myfunc(pars)#thisorderdict[bin]

    for order in range(1,degree+1) :
      code = "self.der{0} = graphs[{0}]".format(order)
      exec code



