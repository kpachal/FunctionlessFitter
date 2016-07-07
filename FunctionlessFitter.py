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

  def function(self,pars) :

    if self.mode == "LogL" :
      return self.computeLogL(self.bincontents,pars)
    elif self.mode == "Chi2" :
      return self.computeChi2(self.bincontents,pars)

  def function_der(self,pars) :
    
    if self.mode == "LogL" :
      return self.jacobianLogL(self.bincontents,pars)
    elif self.mode == "Chi2" :
      return self.jacobianChi2(self.bincontents,pars)

  def boundPositive(self) :
    lowers = []
    uppers = []
    for item in range(len(self.bincontents)) :
      lowers.append(0.0)
      uppers.append(None)
    return zip(lowers,uppers)

  def getStartVals_flat(self) :
    if not self.bincontents :
      print "Fill bins first!"
      return -1
  
    start_vals = []
    for bin in range(len(self.bincontents)) :
      start_vals.append(1.0)
    return start_vals

  def getStartVals_linear(self) :
    if not self.bincontents :
      print "Fill bins first!"
      return -1
    return []
    
  def getStartVals_exponential(self) :
    if not self.bincontents :
      print "Fill bins first!"
      return -1
        
    # Fill with a falling exponential curve of form a e^(bx)
    # which passes through the first point and the last nonzero point.
    
    # First find that point.
    index = -1
    for content in self.bincontents :
      index = index + 1
      if content > 0 :
        x1 = self.binxvals[index]
        y1 = content
        break
    index = len(self.bincontents)
    for content in reversed(self.bincontents) :
      index = index - 1
      if content > 0 :
        x2 = self.binxvals[index]
        y2 = content
        break

    # Now calculate a and b from x1 and x2:
    # b = (ln y1 - ln y2)/(x1 - x2)
    # a = y2 e^(- b x2)
    # y = y2 e^(b(x - x2))
    # So really we only need to compute b.
    print "x1,x2:",x1, x2
    b = (numpy.log(y1) - numpy.log(y2))/(x1 - x2)

    start_vals = []
    for bin in range(len(self.bincontents)) :
      x = self.binxvals[bin]
      y = y2*numpy.exp(b*(x - x2))
      start_vals.append(y)
    return start_vals

  def getConstraints_monotonicity(self,nPars,slope = -1) :
    constraints = []

    for index in range(nPars-1) :
    
      # Write jacobians for this
      jacString = "["
      for item in range(nPars) :
        if item == index :
          if slope < 0 :
            jacString = jacString + "1.0, "
          else :
            jacString = jacString + "-1.0, "
        elif item == index+1 :
          if slope < 0 :
            jacString = jacString + "-1.0, "
          else :
            jacString = jacString + "1.0, "
        else :
          jacString = jacString + "0.0, "
      jacString = jacString+"]"

      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: numpy.array([pars[{1}] - pars[{2}]]),
        'jac' : lambda pars: numpy.array({3}){4}"""
      if slope < 0 :
        code = code.format("{",index,index+1,jacString,"}")
      else :
        code = code.format("{",index+1,index,jacString,"}")
      exec code
      constraints.append(constraint)
  
    return constraints

  def getConstraints_1stDerSmooth(self,nPars,slope=1) :
    constraints = []

    for index in range(nPars-2) :
    
      # Find the two run-values for my slopes
      run1 = self.binxvals[index+1] - self.binxvals[index]
      run2 = self.binxvals[index+2] - self.binxvals[index+1]
    
      # Write jacobians for these:
      # they follow the pattern [1/run2, - (1/run1 + 1/run2), 1/run1]
      # for the three points of interest.
      jacString = "["
      for item in range(nPars) :
        if item == index :
          jacString = jacString + "{0}, ".format(float(slope)/run1)
        elif item == index+1 :
          jacString = jacString + "{0}, ".format(- 1.0 *float(slope) * (1.0/run1 + 1.0/run2))
        elif item == index+2 :
          jacString = jacString + "{0}, ".format(float(slope)/run2)
        else :
          jacString = jacString + "0.0, "
      jacString = jacString+"]"

      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: numpy.array([1.0/{1} * (pars[{5}]- pars[{4}]) - 1.0/{2} * (pars[{4}] - pars[{3}])]),
        'jac' : lambda pars: numpy.array({6}){7}""".format("{",run1,run2,index,index+1,index+2,jacString,"}")
      exec code
      constraints.append(constraint)
  
    return constraints

  def getConstraints_2ndDerSmooth(self,nPars,slope=1) :
    constraints = []

    for index in range(nPars-3) :
    
      # Have two consecutive bin pairs, each defining
      # a dSlopedX of rise/run between that pair of bin centers
      # Want to constrain such that dSlopedX1 > dSlopedX2
    
      # Find the various run-values for my slopes
      run1 = self.binxvals[index+1] - self.binxvals[index]
      run2 = self.binxvals[index+2] - self.binxvals[index+1]
      run3 = self.binxvals[index+3] - self.binxvals[index+2]
    
      run1prun2 = self.binxvals[index+2] - self.binxvals[index]
      run2prun3 = self.binxvals[index+3] - self.binxvals[index+1]
    
      # dSlopedX1 = (slope2-slope1)/(x2-x1)
      # dXlopedX2 = (slope3-slope2)/(x2-x1)
      # Constraint itself then follows format
      # (slope3 - slope2)/(x3-x2) - (slope2-slope1)/(x2-x1)
      # slope_a = (par[b] - par[a])/run
      # x for each slope is what? how shall I track it?
      # Change in slope over region x1 = over run1+run2
      # Change in slope over region x2 = over run2+run3
      # So constraint:
      # ((par[4]-par[3])/run3 - (par[3]-par[2])/run2)/(run2+run3) - ((par[3]-par[2])/run2 - (par[2]-par[1])/run1)/(run1+run2)
      # (par[4]-par[3])/(run3*(run2+run3)) - (par[3]-par[2])/(run2*(run2+run3)) - (par[3]-par[2])/(run2*(run1+run2)) + (par[2]-par[1])/(run1*(run1+run2))
    
      # Write jacobians for this:
      # they follow the pattern [-1/(run1*(run1+run2)), 1/(run2*(run2+run3)) + 1/(run2*(run1+run2)) + 1/(run1*(run1+run2)),
      #                          -1/(run3*(run2+run3)) - 1/(run2*(run2+run3)) - 1/(run2*(run1+run2)), 1/(run3*(run2+run3))]
      # for the four points of interest.
      jacString = "["
      for item in range(nPars) :
        if item == index :
          jacString = jacString + "{0}, ".format(-1.0*float(slope)/(run1*run1prun2))
        elif item == index+1 :
          jacString = jacString + "{0}, ".format(float(slope)/(run2*run2prun3) + float(slope)/(run2*(run1prun2)) + float(slope)/(run1*run1prun2))
        elif item == index+2 :
          jacString = jacString + "{0}, ".format(-1.0*float(slope)/(run3*run2prun3) - float(slope)/(run2*(run2prun3)) - float(slope)/(run2*run1prun2))
        elif item == index+3 :
          jacString = jacString + "{0}, ".format(float(slope)/(run3*run2prun3))
        else :
          jacString = jacString + "0.0, "
      jacString = jacString+"]"

      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: numpy.array([1.0/{1} * (pars[{5}]- pars[{4}]) - 1.0/{2} * (pars[{4}] - pars[{3}])]),
        'jac' : lambda pars: numpy.array({6}){7}""".format("{",run1,run2,index,index+1,index+2,jacString,"}")
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
    
    self.bincontents = []
    self.binxvals = []
    for bin in range(self.rangeLow,self.rangeHigh+1) :
      self.bincontents.append(spectrum.GetBinContent(bin))
      self.binxvals.append(spectrum.GetBinCenter(bin))
    #start_vals = self.getStartVals_exponential()
    start_vals = self.getStartVals_flat()
    print start_vals

    myBounds = self.boundPositive()
    myConstraints = self.getConstraints_monotonicity(len(start_vals))
    myConstraints = myConstraints + self.getConstraints_1stDerSmooth(len(start_vals))
    myConstraints = myConstraints + self.getConstraints_2ndDerSmooth(len(start_vals))

    print "Beginning fit to vals",self.bincontents
    #  jac=self.function_der,
    status = scipy.optimize.minimize(self.function, start_vals, method='SLSQP', jac=self.function_der, bounds=myBounds, constraints=myConstraints, options={'disp': True, 'maxiter':100000, })
    print status

    # Check that the 2nd derivative restriction worked correctly
    print "slopes:",
    for item in range(len(status.x)-1) :
      print round((status.x[item+1]-status.x[item])/(self.binxvals[item]-self.binxvals[item+1]),6),
    
    # Return a histogram with bin contents equal to fit results
    outputHist = spectrum.Clone("fitResult")
    outputHist.Reset()
    index = -1
    for bin in range(self.rangeLow,self.rangeHigh+1) :
      index = index+1
      outputHist.SetBinContent(bin,status.x[index])
    return outputHist
    
  def computeLogL(self,obs,exp) :

    # Likelihood for us is product of Poisson probability in each bin
    #   L = \Prod_{bins} (b^d * e^{-b})/d!
    # Log of this is:
    #   log L = \Sum_{bins} (log(b^2 * e^{-b}/d!))
    #         = \sum_{bins} (d ln b - b - ln(Gamma(d+1)))
    
    answer = 0
    for index in range(len(obs)) :

      if self.excludeWindow and index > self.firstBinToExclude and index < self.lastBinToExclude :
        continue

      data = int(obs[index])
      bkg = exp[index]
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
      bkg = exp[index]
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

