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

  def getConstraints_monotonicity(self,nPars,slope = -1) :
    constraints = []

    for index in range(nPars-1) :

      w0 = self.selectedbinwidths[index]
      w1 = self.selectedbinwidths[index+1]

      # slope * (p1/w1 - p0/w0) > 0
      eqString = "["
      if slope < 0 : eqString = eqString + " - "
      eqString = eqString+ "pars[{1}]/{3}"
      if slope < 0 : eqString = eqString + " + "
      else : eqString = eqString + " - "
      eqString = eqString+"pars[{0}]/{2}"
      eqString = eqString+"]"
      eqString = eqString.format(index,index+1,w0,w1)
      
      # Write jacobians for this
      # They follow the pattern: - slope/w0, slope/w1
      jacString = "["
      for item in range(nPars) :
        if item == index :
            jacString = jacString + "{0}, ".format(float(slope)*-1.0/w0)
        elif item == index+1 :
            jacString = jacString + "{0}, ".format(float(slope)/w1)
        else :
          jacString = jacString + "0.0, "
      jacString = jacString+"]"

      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: numpy.array({1}),
        'jac' : lambda pars: numpy.array({2}){3}"""
      code = code.format("{",eqString,jacString,"}")
      exec code
      constraints.append(constraint)
  
    return constraints

  def getConstraints_1stDerSmooth(self,nPars,slope=1) :
    constraints = []

    for index in range(nPars-2) :
    
      # Find the two run-values for my slopes
      w0 = self.selectedbinwidths[index]
      w1 = self.selectedbinwidths[index+1]
      w2 = self.selectedbinwidths[index+2]
      #run1 = w0 + w1
      #run2 = w1 + w2
      run1 = self.selectedbinxvals[index+1] - self.selectedbinxvals[index]
      run2 = self.selectedbinxvals[index+2] - self.selectedbinxvals[index+1]
      
      # Use finite difference formula to do this:
      # D = (val(up) - val(down))/(run)
      # will correspond to derivative at transition between bins.
      # So for us D2 - D1 > 0:
      # (p2/w2 - p1/w1)/run2 - (p1/w1 - p0/w0)/run1 > 0
    
      eqString = "["
      if slope < 0 : eqString = eqString + " - "
      eqString = eqString+ "(pars[{2}]/{7} - pars[{1}]/{6})/{4}"
      if slope < 0 : eqString = eqString + " + "
      else : eqString = eqString + " - "
      eqString = eqString+"(pars[{1}]/{6} - pars[{0}]/{5})/{3}"
      eqString = eqString+"]"
      eqString = eqString.format(index,index+1,index+2,run1,run2,w0,w1,w2)
    
      # Write jacobians for these:
      # they follow the pattern [slope*1/(w0*run1), - (1/(w1*run1) + 1/(w2*run2)), 1/(w2*run2)]
      # for the three points of interest.
      jacString = "["
      for item in range(nPars) :
        if item == index :
          jacString = jacString + "{0}, ".format(float(slope)/(w0*run1))
        elif item == index+1 :
          jacString = jacString + "{0}, ".format(- 1.0 *float(slope)/w1 * (1.0/run1 + 1.0/run2))
        elif item == index+2 :
          jacString = jacString + "{0}, ".format(float(slope)/(run2*w2))
        else :
          jacString = jacString + "0.0, "
      jacString = jacString+"]"

      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: numpy.array({1}),
        'jac' : lambda pars: numpy.array({2}){3}""".format("{",eqString,jacString,"}")
      exec code
      constraints.append(constraint)
  
    return constraints

#  def getConstraints_2ndDerSmooth(self,nPars,slope=1) :
#    constraints = []
#
#    for index in range(nPars-3) :
#    
#      # Have two consecutive bin pairs, each defining
#      # a second derivative between that pair of bin centers
#      # Want to constrain such that dSlopedX1 > dSlopedX2
#    
#      # Find the various run-values for my slopes
#      run1 = self.selectedbinxvals[index+1] - self.selectedbinxvals[index]
#      run2 = self.selectedbinxvals[index+2] - self.selectedbinxvals[index+1]
#      run3 = self.selectedbinxvals[index+3] - self.selectedbinxvals[index+2]
#      
#      run1prun2 = self.selectedbinxvals[index+2] - self.selectedbinxvals[index]
#      run2prun3 = self.selectedbinxvals[index+3] - self.selectedbinxvals[index+1]
#
#      w0 = self.selectedbinwidths[index]
#      w1 = self.selectedbinwidths[index+1]
#      w2 = self.selectedbinwidths[index+2]
#      w3 = self.selectedbinwidths[index+3]
#
#      # Use finite difference formulas to do this.
#      # D = (y2 - y1)/((x2 - x1)*(x2 - x0)) - (y1 - y0)/((x1 - x0)*(x2-x0)) > 0
#      # For us, 2nd derivative should be monotonically decreasing: slope = -1
#      # therefore, D2 - D1 > 0. Inequality is:
#      # slope * (y3 - y2)/((x3 - x2)*(x3 - x1)) - (y2 - y1)/((x2 - x1)*(x3-x1)) -
#      # slope * (y2 - y1)/((x2 - x1)*(x2 - x0)) - (y1 - y0)/((x1 - x0)*(x2-x0))
#      # keeping in mind that y3 = b3/w3 etc
#    
#      # Write an equation string for this.
#
#      rise1 = "pars[{1}]/{3} - pars[{0}]/{2}".format(index,index+1,w0,w1)
#      rise2 = "pars[{1}]/{3} - pars[{0}]/{2}".format(index+1,index+2,w1,w2)
#      rise3 = "pars[{1}]/{3} - pars[{0}]/{2}".format(index+2,index+3,w2,w3)
#    
#      dSlopedX1 = "(({1}*100.0/{3}) - ({0}*100.0/{2}))/({4})".format(rise1,rise2,run1,run2,run1prun2)
#      dSlopedX2 = "(({1}*100.0/{3}) - ({0}*100.0/{2}))/({4})".format(rise2,rise3,run2,run3,run2prun3)
#    
#      eqString = "["
#      if slope < 0 : eqString = eqString + " - "
#      eqString = eqString + dSlopedX1
#      if slope < 0 : eqString = eqString + " + "
#      else : eqString = eqString + " - "
#      eqString = eqString + dSlopedX2
#      eqString = eqString + "]"
#      
#      # Write jacobians for this:
#      # d/dy0 = -slope/(run1*(run1+run2))
#      # d/dy1 = 1/(run2*(run2+run3)) + 1/(run2*(run1+run2)) + 1/(run1*(run1+run2))
#      # d/dy2 =
#      # d/dy3 =
#      # they follow the pattern [+ 1/(run1*(run1+run2)), - 1/(run2*(run1+run2)) - 1/(run1*(run1+run2)) - 1/(run2*(run2+run3)),
#      #                          1/(run3*(run2+run3)) + 1/(run2*(run2+run3)) + 1/(run2*(run1+run2)), - 1/(run3*(run2+run3))]
#      # for the four points of interest.
#      jacString = "["
#      for item in range(nPars) :
#        if item == index :
#          jacString = jacString + "{0}, ".format(float(slope)/(run1*run1prun2))
#        elif item == index+1 :
#          jacString = jacString + "{0}, ".format(-1.0*float(slope)/(run1*run1prun2) + float(slope)/(run2*(run1prun2)) - float(slope)/(run3*run2prun3))
#        elif item == index+2 :
#          jacString = jacString + "{0}, ".format(float(slope)/(run2*run2prun3) - float(slope)/(run2*(run1prun2)) + float(slope)/(run3*run2prun3))
#        elif item == index+3 :
#          jacString = jacString + "{0}, ".format(-1.0*float(slope)/(run3*run2prun3))
#        else :
#          jacString = jacString + "0.0, "
#      jacString = jacString+"]"
#
##      code = """constraint = {0}'type': 'ineq',
##        'fun' : lambda pars: numpy.array([(pars[{8}]-pars[{7}])/({2}*{4}) + (pars[{7}]-pars[{6}])/({1}*{4}) - (pars[{9}]-pars[{8}])/({3}*{5}) - (pars[{8}]-pars[{7}])/({2}*{5})]),
##        'jac' : lambda pars: numpy.array({10}){11}""".format("{",run1,run2,run3,run1prun2,run2prun3,index,index+1,index+2,index+3,jacString,"}")
#      code = """constraint = {0}'type': 'ineq',
#        'fun' : lambda pars: numpy.array({1}){2}""".format("{",eqString,"}") #,
#        #'jac' : lambda pars: numpy.array({2}){3}""".format("{",eqString,jacString,"}")
#      exec code
#      constraints.append(constraint)
#
#      print constraints

#    return constraints

  def getConstraints_2ndDerSmooth(self,nPars,slope=-1) :
    constraints = []

    for index in range(nPars-3) :
    
      # Have two consecutive bin pairs, each defining
      # a second derivative between that pair of bin centers
      # Want to constrain such that dSlopedX1 > dSlopedX2
    
      # Find the various run-values for my slopes
      run1 = self.selectedbinxvals[index+1] - self.selectedbinxvals[index]
      run2 = self.selectedbinxvals[index+2] - self.selectedbinxvals[index+1]
      run3 = self.selectedbinxvals[index+3] - self.selectedbinxvals[index+2]
      
      run1prun2 = self.selectedbinxvals[index+2] - self.selectedbinxvals[index]
      run2prun3 = self.selectedbinxvals[index+3] - self.selectedbinxvals[index+1]

      w0 = self.selectedbinwidths[index]
      w1 = self.selectedbinwidths[index+1]
      w2 = self.selectedbinwidths[index+2]
      w3 = self.selectedbinwidths[index+3]

      # Use finite difference formulas to do this.
      # D = (y2 - y1)/((x2 - x1)*(x2 - x0)) - (y1 - y0)/((x1 - x0)*(x2-x0)) > 0
      # For us, 2nd derivative should be monotonically decreasing: slope = -1
      # therefore, D2 - D1 > 0. Inequality is:
      # slope * (y3 - y2)/((x3 - x2)*(x3 - x1)) - (y2 - y1)/((x2 - x1)*(x3-x1)) -
      # slope * (y2 - y1)/((x2 - x1)*(x2 - x0)) - (y1 - y0)/((x1 - x0)*(x2-x0))
      # keeping in mind that y3 = b3/w3 etc
    
      # Write an equation string for this.
      eqString = "["
      if slope < 0 : eqString = eqString + " - "
      eqString = eqString + "((pars[{3}]/{12} - pars[{2}]/{11})/{6} - (pars[{2}]/{11} - pars[{1}]/{10})/{5})/{8}"
      if slope < 0 : eqString = eqString + " + "
      else : eqString = eqString + " - "
      eqString = eqString + " ((pars[{2}]/{11} - pars[{1}]/{10})/{5} - (pars[{1}]/{10} - pars[{0}]/{9})/{4})/{7}"
      eqString = eqString + "]"
      eqString = eqString.format(index,index+1,index+2,index+3,run1,run2,run3,run1prun2,run2prun3,w0,w1,w2,w3)

      # Write jacobians for this:
      # d/dy0 = -slope/(run1*(run1+run2))
      # d/dy1 = 1/(run2*(run2+run3)) + 1/(run2*(run1+run2)) + 1/(run1*(run1+run2))
      # d/dy2 =
      # d/dy3 =
      # they follow the pattern [+ 1/(run1*(run1+run2)), - 1/(run2*(run1+run2)) - 1/(run1*(run1+run2)) - 1/(run2*(run2+run3)),
      #                          1/(run3*(run2+run3)) + 1/(run2*(run2+run3)) + 1/(run2*(run1+run2)), - 1/(run3*(run2+run3))]
      # for the four points of interest.
      jacString = "["
      for item in range(nPars) :
        if item == index :
          jacString = jacString + "{0}, ".format(float(slope)/(run1*run1prun2))
        elif item == index+1 :
          jacString = jacString + "{0}, ".format(-1.0*float(slope)/(run1*run1prun2) + float(slope)/(run2*(run1prun2)) - float(slope)/(run3*run2prun3))
        elif item == index+2 :
          jacString = jacString + "{0}, ".format(float(slope)/(run2*run2prun3) - float(slope)/(run2*(run1prun2)) + float(slope)/(run3*run2prun3))
        elif item == index+3 :
          jacString = jacString + "{0}, ".format(-1.0*float(slope)/(run3*run2prun3))
        else :
          jacString = jacString + "0.0, "
      jacString = jacString+"]"

#      code = """constraint = {0}'type': 'ineq',
#        'fun' : lambda pars: numpy.array([(pars[{8}]-pars[{7}])/({2}*{4}) + (pars[{7}]-pars[{6}])/({1}*{4}) - (pars[{9}]-pars[{8}])/({3}*{5}) - (pars[{8}]-pars[{7}])/({2}*{5})]),
#        'jac' : lambda pars: numpy.array({10}){11}""".format("{",run1,run2,run3,run1prun2,run2prun3,index,index+1,index+2,index+3,jacString,"}")
      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: numpy.array({1}){2}""".format("{",eqString,"}") #,
        #'jac' : lambda pars: numpy.array({2}){3}""".format("{",eqString,jacString,"}")
      code
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
    print start_vals

    myBounds = self.boundPositive()
    myConstraints = self.getConstraints_monotonicity(len(start_vals))
    myConstraints = myConstraints + self.getConstraints_1stDerSmooth(len(start_vals))
#    myConstraints = self.getConstraints_1stDerSmooth(len(start_vals))
    myConstraints = myConstraints + self.getConstraints_2ndDerSmooth(len(start_vals))

    print "Beginning fit to vals",self.selectedbincontents
    #  jac=self.function_der,
    status = scipy.optimize.minimize(self.function, start_vals, method='SLSQP', jac=self.function_der, bounds=myBounds, constraints=myConstraints, options={'disp': True, 'maxiter':100000, })
    print status

    # Check that the 1st derivative restriction worked correctly
    print "slopes:",
    for item in range(len(status.x)-1) :
      print round((status.x[item+1]-status.x[item])/(self.selectedbinxvals[item]-self.selectedbinxvals[item+1]),6),
    
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

