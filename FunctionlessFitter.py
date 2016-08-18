import ROOT
import numpy
import scipy
from scipy import optimize
from scipy import special
from dill.source import getsource
from HistWrapper import WrappedHist

class FunctionlessFitter :

  def __init__(self) :
    self.mode = "LogL"
    self.excludeWindow = False
    self.result = []
    self.dividedDifferenceDatabase = {}
    self.jacobianDatabase = {}
  
  def function(self,pars) :

    if self.mode == "LogL" :
      return self.computeLogL(self.selectedbincontents,pars)
    elif self.mode == "Chi2" :
      return self.computeChi2(self.selectedbincontents,pars)

  def function_der(self,pars) :
    
    if self.mode == "LogL" :
      return self.jacobianLogL(self.selectedbincontents,pars)
    elif self.mode == "Chi2" :
      return self.jacobianChi2(self.selectedbincontents,pars)

  def boundPositive(self) :
    lowers = []
    uppers = []
    for item in range(len(self.selectedbincontents)) :
      lowers.append(0.0)
      uppers.append(None)
    return zip(lowers,uppers)

  def getStartVals_flat(self) :
    if not self.selectedbincontents :
      print "Fill bins first!"
      return -1
  
    start_vals = []
    for bin in range(len(self.selectedbincontents)) :
      start_vals.append(1.0)
    return start_vals

  def getStartVals_linear(self) :
    if not self.selectedbincontents :
      print "Fill bins first!"
      return -1
    return []
    
  def getStartVals_exponential(self) :
    if not self.selectedbincontents :
      print "Fill bins first!"
      return -1
        
    # Fill with a falling exponential curve of form a e^(bx)
    # which passes through the first point and the last nonzero point.
    
    # First find that point.
    index = -1
    for content in self.selectedbincontents :
      index = index + 1
      if content > 0 :
        x1 = self.selectedbinxvals[index]
        y1 = content
        break
    index = len(self.selectedbincontents)
    for content in reversed(self.selectedbincontents) :
      index = index - 1
      if content > 0 :
        x2 = self.selectedbinxvals[index]
        y2 = content
        break

    # Now calculate a and b from x1 and x2:
    # b = (ln y1 - ln y2)/(x1 - x2)
    # a = y2 e^(- b x2)
    # y = y2 e^(b(x - x2))
    # So really we only need to compute b.
    b = (numpy.log(y1) - numpy.log(y2))/(x1 - x2)

    start_vals = []
    for bin in range(len(self.selectedbincontents)) :
      x = self.selectedbinxvals[bin]
      y = y2*numpy.exp(b*(x - x2))
      start_vals.append(y)
    return start_vals

  def computeConstraints(self, degree) :

    # 0th degree derivatives
    baseDict = {}
    for bin in range(len(self.selectedbinxvals)) :
      baseDict[bin] = "pars[{0}]".format(bin)
    self.dividedDifferenceDatabase[0] = baseDict
  
    # 0th degree Jacobian matrix
    baseJac = numpy.identity(len(self.selectedbinxvals))
    self.jacobianDatabase[0] = baseJac
  
    # higher order derivatives and jacobians
    for order in range(1,degree+1) :
    
      lastorderdict = self.dividedDifferenceDatabase[int(order-1)]
      thisorderdict = {}
 
      # Matrix to multiply into current jacobian terms
      A = []

      for index in range(len(self.selectedbinxvals) - order) :
      
        # Difference for the constraint itself
        fxa = lastorderdict[index]
        fxb = lastorderdict[index+1]
        diff = "({1}-{0})/({3}-{2})".format(fxa, fxb, self.selectedbinxvals[index],self.selectedbinxvals[index+order])
        thisorderdict[index] = diff

        jacRow = []
        for column in range(len(self.selectedbinxvals) - order+1) :
          val = self.selectedbinxvals[index+order] - self.selectedbinxvals[index]
          if column == index :
            jacRow.append(-1.0/val)
          elif column == index+1 :
            jacRow.append(1.0/val)
          else :
            jacRow.append(0)
        A.append(jacRow)

      self.dividedDifferenceDatabase[int(order)] = thisorderdict
      self.jacobianDatabase[int(order)] = numpy.dot(A,self.jacobianDatabase[int(order)-1])

  
  def getDerivativeConstraints(self, degree, slope) :

    if int(degree) not in self.dividedDifferenceDatabase.keys() :
      self.computeConstraints(degree)

    constraints = []
    nPars = len(self.selectedbincontents)
    
    for index in range(nPars-degree-1) :
    
      term1 = self.dividedDifferenceDatabase[degree][index]
      term2 = self.dividedDifferenceDatabase[degree][index+1]
    
      eqString = "["
      if slope < 0 : eqString = eqString + " - "
      eqString = eqString + term2
      if slope < 0 : eqString = eqString + " + "
      else : eqString = eqString + " - "
      eqString = eqString + term1
      eqString = eqString + "]"
    
      vec1 = numpy.array(self.jacobianDatabase[degree][index])
      vec2 = numpy.array(self.jacobianDatabase[degree][index+1])
      jacobian = slope*vec2 - slope*vec1
      jacString = "["
      for term in jacobian :
        jacString = jacString + "{0}, ".format(term)
      jacString = jacString + "]"

      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: numpy.array({1}),
        'jac' : lambda pars: numpy.array({2}){3}""".format("{",eqString,jacString,"}")
      exec code
      constraints.append(constraint)

    return constraints

  def fit(self,spectrum,firstBin=-1,lastBin=-1) :

    self.myHistogram = WrappedHist(spectrum)

    if firstBin < 0 or firstBin > spectrum.GetNbinsX() :
      self.rangeLow = self.myHistogram.firstBinWithData
    else : self.rangeLow = firstBin
    if lastBin < 0 or lastBin > spectrum.GetNbinsX() or lastBin < firstBin :
      self.rangeHigh = self.myHistogram.lastBinWithData
    else : self.rangeHigh = lastBin
    
    self.selectedbincontents, self.selectedbinxvals, self.selectedbinwidths = self.myHistogram.getSelectedBinInfo(self.rangeLow,self.rangeHigh)
    start_vals = self.getStartVals_exponential()
    #start_vals = self.getStartVals_flat()

    myBounds = self.boundPositive()
    
    # Calculate values to use for constraints: saves us doing it later
    self.computeConstraints(4)
    myConstraints = self.getDerivativeConstraints(0,-1)
    myConstraints = myConstraints + self.getDerivativeConstraints(1,1)
    #myConstraints = myConstraints + self.getDerivativeConstraints(2,-1)
    #myConstraints = myConstraints + self.getDerivativeConstraints(3,1)
    #myConstraints = myConstraints + self.getDerivativeConstraints(4,-1)

    print "Beginning fit to vals",self.selectedbincontents
    status = scipy.optimize.minimize(self.function, start_vals, method='SLSQP', jac=self.function_der, bounds=myBounds, constraints=myConstraints, options={'disp': True, 'maxiter':100000, })
    #print status

    self.result = status.x

    # Return a histogram with bin contents equal to fit results
    outputHist = spectrum.Clone("fitResult")
    outputHist.Reset()
    index = -1
    for bin in range(self.rangeLow,self.rangeHigh+1) :
      index = index+1
      outputHist.SetBinContent(bin,status.x[index]*self.selectedbinwidths[index])
    return outputHist
    
  def computeLogL(self,obs,exp) :

    # Likelihood for us is product of Poisson probability in each bin
    #   L = \Prod_{bins} (b^d * e^{-b})/d!
    # Log of this is:
    #   log L = \Sum_{bins} (log(b^2 * e^{-b}/d!))
    #         = \sum_{bins} (d ln b - b - ln(Gamma(d+1)))
    
    answer = 0
    for index in range(len(obs)) :

      #if self.excludeWindow and index > self.firstBinToExclude and index < self.lastBinToExclude :
      #  continue

      data = int(obs[index])
      bkg = exp[index]*self.selectedbinwidths[index]
      if data < 0.0 or bkg < 0.0 :
        thisterm = -1E10
      elif data == 0.0 :
        thisterm = -1.0*bkg
      else :
        thisterm = data * numpy.log(bkg) - bkg - scipy.special.gammaln(data+1)
      answer = answer - thisterm
    
    return answer


  def computeChi2(self,obs,exp) :

    return

  def jacobianLogL(self,obs,exp) :

    answer = []
    sum = 0
    for index in range(len(obs)) :

      if self.excludeWindow and index > self.firstBinToExclude and index < self.lastBinToExclude :
        continue

      data = int(obs[index])
      bkg = exp[index]*self.selectedbinwidths[index]
      if data < 0.0 or bkg < 0.0:
        thisterm = 0.0
      elif data == 0.0 :
        thisterm = 1.0
      else :
        thisterm = 1.0 - (float(data)/float(bkg))
      sum = sum+thisterm
      answer.append(thisterm)

    return numpy.array(answer)

  def jacobianChi2(self,obs,exp) :
  
    return

