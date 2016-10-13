import ROOT
import numpy
import scipy

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
