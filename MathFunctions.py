import numpy
import scipy

def PoissonPVal(data, bkg) :

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

def PoissonConvGammaPVal(data, bkg, bkgErr) :

  return PoissonPVal(data, bkg)

def ProbToSigma(prob) :

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

def SigmaToProb(sigma) :

  return 0.5*(1.0 - scipy.special.erf(sigma/sqrt(2.0)))