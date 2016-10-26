import ROOT
import numpy
import scipy
from decimal import *

def computeDividedDifferences(degree,selectedbinxvals,selectedbinedges,scaleParsBy) :

    dividedDifferenceDatabase = {}
    jacobianDatabase = {}

    # 0th degree derivatives
    baseDict = {}
    for bin in range(len(selectedbinxvals)) :
      baseDict[bin] = "Decimal(pars[{0}]*{1})".format(bin,scaleParsBy[bin])
      #baseDict[bin] = "(pars[{0}]*{1})".format(bin,self.scaleParsBy[bin])
    dividedDifferenceDatabase[0] = baseDict

    # 0th degree Jacobian matrix
    baseJac = numpy.identity(len(selectedbinxvals),dtype=Decimal)*scaleParsBy

    # higher order derivatives and jacobians
    for order in range(1,degree+1) :
    
      lastorderdict = dividedDifferenceDatabase[int(order-1)]
      thisorderdict = {}
 
      # Matrix to multiply into current jacobian terms
      A = []

      for index in range(len(selectedbinxvals) - order) :
      
        # Difference for the constraint itself
        fxa = lastorderdict[index]
        fxb = lastorderdict[index+1]
        # Number of bins = order+1
        # So relevant bin = int(order/2)+1
        # Even orders use bin centers,
        # while odd orders use bin lower edge.
        # However, counting from bin = index, so don't need the 1.
        relevantbina = index+int(order/2)
        relevantbinb = index+int(order/2)+1
        if order%2==0 :
          xa = selectedbinxvals[relevantbina]
          xb = selectedbinxvals[relevantbinb]
        else :
          xa = selectedbinedges[relevantbina]
          xb = selectedbinedges[relevantbinb]
          
        diff = "({1}-{0})/({3}-{2})".format(fxa, fxb, "Decimal({0})".format(selectedbinxvals[index]),"Decimal({0})".format(selectedbinxvals[index+order]))
        #diff = "({1}-{0})/({3}-{2})".format(fxa, fxb, "Decimal({0})".format(xa),"Decimal({0})".format(xb))
        thisorderdict[index] = diff

        jacRow = []
        for column in range(len(selectedbinxvals) - order+1) :
          val = selectedbinxvals[index+order] - selectedbinxvals[index]
          if column == index :
            jacRow.append(Decimal(-1.0)/val)
          elif column == index+1 :
            jacRow.append(Decimal(1.0)/val)
          else :
            jacRow.append(Decimal(0.0))
        A.append(jacRow)

      dividedDifferenceDatabase[int(order)] = thisorderdict
      jacobianDatabase[int(order)] = numpy.dot(A,jacobianDatabase[int(order)-1])

    return dividedDifferenceDatabase,jacobianDatabase

def poissonPVal(data, bkg) :

  # Scipy has a function which is the sum of the first k terms
  # of the Poisson distribution: sum(exp(-m) * m**j / j!, j=0..k).
  
  answer = 1.0
  if (data < bkg) :
    # Sum downwards. Want the probability of getting d or fewer events given b.
    # That is, the first d terms of the Poisson distribution.
    answer = scipy.special.pdtr(data, bkg)
  else :
    # Sum upwards. Want the probability of getting d or more events.
    # That is, 1 minus the first (d-1) terms of the Poisson distribution.
    answer = 1.0 - scipy.special.pdtr(data-1,bkg)

  return answer

def poissonConvGammaPVal(data, bkg, bkgErr) :

  return PoissonPVal(data, bkg)

def probToSigma(prob) :

  assert(prob >=0 and prob <=1)
  
  # p = 0.5 - 0.5*erf(s/sqrt(2.0))
  # 2p = 1 - erf(s/sqrt(2.0))
  # erf(s/sqrt(2.0)) = 1-2p
  # s/sqrt(2.0) = erfInverse(1-2p)
  # s = sqrt(2.0) * erfInverse(1-2p)
  
  value = 1.0-2.0*prob
  if value>-1 and value<1 : return numpy.sqrt(2.0)*scipy.special.erfinv(value)
  elif (value==1) : return 1E10
  else : return -1E10

def sigmaToProb(sigma) :

  return 0.5*(1.0 - scipy.special.erf(sigma/sqrt(2.0)))

def makeHistFromVector(vector,alternateVal=-1) :

  nentries = len(vector)
  if alternateVal > 0 :
    nBins = int(float(nentries)/alternateVal)
  else :
    nBins = int(float(nentries)/10.0)

  maxVal = float(max(vector))
  minVal = float(min(vector))
  range = maxVal - minVal

  plotmin = minVal-0.05*range
  plotmax = maxVal+0.05*range

  statPlot = ROOT.TH1D("statPlot","",nBins,plotmin,plotmax)
  for item in vector : statPlot.Fill(item)

  return statPlot

def getPValFromVecAndStat(stat,vector) :

  vec = sorted(vector)
  nBelow = 0
  for item in vec :
    if item < stat :
      nBelow = nBelow+1
    else : break

  pVal = float(len(vector)-nBelow)/float(len(vector))
  return pVal

def getFlatVector(length, val) :
    return [Decimal(val)]*length
