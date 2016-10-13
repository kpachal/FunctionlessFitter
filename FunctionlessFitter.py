import ROOT
import numpy
import scipy
from scipy import optimize
from scipy import special
from HistWrapper import WrappedHist
from fractions import *
from decimal import *
getcontext().prec = 28
from numbers import *

class FunctionlessFitter :

  def __init__(self) :
    self.mode = "LogL"
    self.excludeWindow = False
    self.result = []
    self.dividedDifferenceDatabase = {}
    self.jacobianDatabase = {}
    self.firstBinInWindow = -1
    self.lastBinInWindow = -1
    self.userStartVals = []
    self.nPEs = 10
    
    # Fitter internal parameters
    
    # Options: "SLSQP","COBYLA"
    #self.minAlg = "COBYLA"
    self.minAlg = "SLSQP"
    
    # Currently supported: "exp", "flat", "linear"
    self.startValFormat = "exp"
    self.flatStartVal = 1.0
  
    # Keys are order of derivative to restrict,
    # mapped values are sign of slope of derivative
    self.derivativeConstraints = {0:-1, 1:1, 2:-1, 3:1}
  
  def getOptionsDict(self,algorithm) :
    if algorithm == "SLSQP" :
      # 'disp': set verbosity
      # 'maxiter': keep high for diagnostics
      # 'ftol': stopping precision. Default is 1e-6
      # 'eps': step size for approximating the Jacobean (we don't have
      #        to do this because we know our Jacobean)
      options={'disp': True, 'maxiter':10000, 'ftol':1e-9} # 1e-6 is default
          # Testing in restricted TLA.
          # Does not converge with: 1e-10.
          # Does with 1e-9 and constraints up to 3.
          # Constraints up to 4: failed when precision 1e-9 required. 1e-8: no. 1e-7: no. 1e-6: yes.
          # Up to 5: worked at 1e-6.
          # Up to 6: worked at 1e-6 but the residuals were more narrow afterwards, although the function value was larger than for 5
          # Up to 7: works at 1e-6 but very slowly. width of dist:
          # Conclusion: default precision works for constraints from 4 through 6. Can go as tight as 1e-9 for constraints in 0-3. 7 and above require loosening beyond the default.
          # Did that change at all as a result of the re-binning?
          # deactivate rebinning and try up to 5 at 1e-6
          # Use 1e-9 to get satisfying identical results for ICHEP without parameter tricks.
    elif algorithm == "COBYLA" :
      # 'disp': set verbosity
      # 'maxiter': keep high for diagnostics
      # 'rhobeg': how far can we move from the start parameters?
      #           If you are using start parameters very far from
      #           the scale of the final out come (e.g. linear at 1)
      #           this needs to be set very high to allow them to change sufficiently.
      #           Right now the "true" final fit is done with scaled parameters
      #           so this doesn't need to be so big. Was 1e6 before.
      # 'tol': tolerance of final result. By default none specified
      # 'catol': tolerance of constraint violations. By default 0.0002
      options={'disp': True, 'maxiter':1000, 'rhobeg':50, 'tol':1e-5, 'catol':1e-5}
    else :
      raise ValueError("Unrecognized minimization algorithm!\nPlease use one of 'SLSQP','COBYLA'")

    return options

  def getFlatVector(self, length, val) :
    return [Decimal(val)]*length

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

  def getStartVals_flat(self,doScale=True) :
    if not self.selectedbincontents :
      print "Fill bins first!"
      return -1
  
    start_vals = []
    self.flatStartVal = Decimal(self.flatStartVal)
    for bin in range(len(self.selectedbincontents)) :
      if doScale :
        start_vals.append(self.flatStartVal/self.scaleParsBy[bin])
      else :
        start_vals.append(self.flatStartVal)
    return start_vals

  def getStartVals_linear(self,doScale=True) :
    if not self.selectedbincontents :
      print "Fill bins first!"
      return -1
    
    run = Decimal(self.selectedbinxvals[-1] - self.selectedbinxvals[0])
    slope = Decimal((self.selectedbincontents[-1]/self.selectedbinwidths[-1] - self.selectedbincontents[0]/self.selectedbinwidths[0])/run)
    start_vals = []
    for bin in range(len(self.selectedbincontents)) :
      val = self.selectedbincontents[0]/self.selectedbinwidths[0] + slope*(self.selectedbinxvals[bin]-self.selectedbinxvals[0])
      if doScale :
        start_vals.append(val/self.scaleParsBy[bin])
      else :
        start_vals.append(val)
    return start_vals
    
  def getStartVals_exponential(self,doScale=True) :
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
        y1 = content/self.selectedbinwidths[index]
        break
    index = len(self.selectedbincontents)
    for content in reversed(self.selectedbincontents) :
      index = index - 1
      if content > 0 :
        x2 = self.selectedbinxvals[index]
        y2 = content/self.selectedbinwidths[index]
        break

    # Now calculate a and b from x1 and x2:
    # b = (ln y1 - ln y2)/(x1 - x2)
    # a = y2 e^(- b x2)
    # y = y2 e^(b(x - x2))
    # So really we only need to compute b.
    b = (y1.ln() - y2.ln())/(x1 - x2)

    start_vals = []
    for bin in range(len(self.selectedbincontents)) :
      x = self.selectedbinxvals[bin]
      y = y2*numpy.exp(b*(x - x2))
      if doScale :
        start_vals.append(y/self.scaleParsBy[bin])
      else :
        start_vals.append(y)
    return start_vals
  
  def computeConstraints(self,degree) :

    # 0th degree derivatives
    baseDict = {}
    for bin in range(len(self.selectedbinxvals)) :
      baseDict[bin] = "Decimal(pars[{0}]*{1})".format(bin,self.scaleParsBy[bin])
      #baseDict[bin] = "(pars[{0}]*{1})".format(bin,self.scaleParsBy[bin])
    self.dividedDifferenceDatabase[0] = baseDict
  
    # 0th degree Jacobian matrix
    baseJac = numpy.identity(len(self.selectedbinxvals),dtype=Decimal)*self.scaleParsBy
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
        diff = "Decimal({1}-{0})/Decimal({3}-{2})".format(fxa, fxb, self.selectedbinxvals[index],self.selectedbinxvals[index+order])

        thisorderdict[index] = diff

        jacRow = []
        for column in range(len(self.selectedbinxvals) - order+1) :
          val = self.selectedbinxvals[index+order] - self.selectedbinxvals[index]
          if column == index :
            jacRow.append(Decimal(-1.0)/val)
          elif column == index+1 :
            jacRow.append(Decimal(1.0)/val)
          else :
            jacRow.append(Decimal(0.0))
        A.append(jacRow)

      self.dividedDifferenceDatabase[int(order)] = thisorderdict
      self.jacobianDatabase[int(order)] = numpy.dot(A,self.jacobianDatabase[int(order)-1])

  def getDerivativeConstraints(self, degree, slope) :

    if int(degree) not in self.dividedDifferenceDatabase.keys() :
      self.computeConstraints(degree)

    constraints = []
    self.eqDict = {}
    nPars = len(self.selectedbincontents)
    
    for index in range(nPars-degree-1) :
    
      term1 = self.dividedDifferenceDatabase[degree][index]
      term2 = self.dividedDifferenceDatabase[degree][index+1]
    
      # At this point it's possible to compare scalePars of terms that ended up relevant
      # and scale constraints and jacobians down by them. Let's see if that
      # stabilises the convergence in the scaled cases.
      # We know what indices were relevant in this spot:
#      relevantIndices = range(index,index+degree+2)
#      relevantScales = self.scaleParsBy[index:index+degree+2]
#      minScale = min(relevantScales)
#      newScales = [x/minScale for x in relevantScales]

      #eqString = "["
      eqString = ""
      if slope < 0 : eqString = eqString + " - "
      eqString = eqString + term2
      if slope < 0 : eqString = eqString + " + "
      else : eqString = eqString + " - "
      eqString = eqString + term1
      #eqString = eqString + "]"
    
      # Implement scaling of eqString to new lowest common denominators
#      for subindex in relevantIndices :
#        eqString = eqString.replace("pars[{0}]*{1}".format(subindex,relevantScales[relevantIndices.index(subindex)]),"pars[{0}]*{1}".format(subindex,newScales[relevantIndices.index(subindex)]))
#        eqString = eqString.replace("(pars[","{0}*(pars[".format(minScale*10.0))

      vec1 = numpy.array(self.jacobianDatabase[degree][index])
      vec2 = numpy.array(self.jacobianDatabase[degree][index+1])
      jacobian = slope*vec2 - slope*vec1
      jacString = "["
      for term in jacobian :
        jacString = jacString + "{0}, ".format(term)
      jacString = jacString + "]"

      self.eqDict[index] = {'eq' : eqString, 'jac' : eval(jacString)}

      code = """constraint = {0}'type': 'ineq',
        'fun' : lambda pars: {1},
        'jac' : lambda pars: numpy.array({2}){3}""".format("{",eqString,jacString,"}")
      exec code
      constraints.append(constraint)

    return constraints

  def fit(self,spectrum,firstBin=-1,lastBin=-1, errType = "None") :

    if firstBin < 0 or firstBin > spectrum.histogram.GetNbinsX() :
      self.rangeLow = spectrum.firstBinWithData
    else : self.rangeLow = firstBin
    if lastBin < 0 or lastBin > spectrum.histogram.GetNbinsX() or lastBin < firstBin :
      self.rangeHigh = spectrum.lastBinWithData
    else : self.rangeHigh = lastBin
    
    self.selectedbincontents, self.selectedbinxvals, self.selectedbinwidths, self.windowLow, self.windowHigh = spectrum.getSelectedBinInfo(self.rangeLow,self.rangeHigh,self.firstBinInWindow,self.lastBinInWindow)
    self.scaleParsBy = spectrum.scaleFactors
    
    if self.startValFormat == "exp" :
      print "Using exponential start values"
      start_vals = self.getStartVals_exponential()
      start_vals_unscaled = self.getStartVals_exponential(False)
    elif self.startValFormat == "flat" :
      print "Using flat start values."
      start_vals = self.getStartVals_flat()
      start_vals_unscaled = self.getStartVals_flat(False)
    elif self.startValFormat == "linear" :
      print "Using linear start values"
      start_vals = self.getStartVals_linear()
      start_vals_unscaled = self.getStartVals_linear(False)
    elif self.startValFormat == "user" :
      if len(self.userStartVals) > 0 :
        print "Using user-specified start values."
        temp_vals = tuple(self.userStartVals)
        start_vals = []
        start_vals_unscaled = []
        index = -1
        for v in temp_vals :
          index = index+1
          start_vals.append(Decimal(v)/self.scaleParsBy[index])
          start_vals_unscaled.append(Decimal(v))
      else :
        print "No start values specified by user!"
        print "Using default flat values."
        start_vals = self.getStartVals_flat()
        start_vals_unscaled = self.getStartVals_flat(False)
    else :
      raise ValueError("Requested start value format is unrecognized!")

    self.myBounds = self.boundPositive()
    
    # Calculate values to use for constraints: saves us doing it later
    orders = self.derivativeConstraints.keys()

    # Unscaled version
    self.scaleParsBy = self.getFlatVector(len(self.selectedbinxvals),1.0)
    self.computeConstraints(max(orders))
    self.unscaledConstraints = []
    for order in orders :
      slope = self.derivativeConstraints[order]
      self.unscaledConstraints = self.unscaledConstraints + self.getDerivativeConstraints(order,slope)

    # Scaled version (default)
    self.scaleParsBy = spectrum.scaleFactors
    self.computeConstraints(max(orders))
    self.myConstraints = []
    for order in orders :
      slope = self.derivativeConstraints[order]
      self.myConstraints = self.myConstraints + self.getDerivativeConstraints(order,slope)

    debug = False

    # Add smoothing routine for user-supplied start values
    # Use COBYLA with loosened tolerance on parameter constraint obedience
    # to smooth, for instance, data values so they obey further iterations
    # rhobeg was 1e4 for non-scaled system but is now lower because of parameter adjustments
    #looseOpts = {'disp': True, 'maxiter':1000, 'rhobeg':50, 'tol':1e-7, 'catol':1e-2} # 2 and 2 good.
    #looseOpts = {'disp': True, 'maxiter':1000, 'rhobeg':1e5, 'tol':1e-6, 'catol':1e-2} # 2 and 2 good.
    looseOpts = {'disp': True, 'maxiter':1000,'catol':1e-2} # 2 and 2 good.
    if self.startValFormat == "user" :
      # If we are using start values within an even vaguely reasonable distance of the
      # final result, we get our best chance by taking smaller steps and combing the
      # parameter space thoroughly. Use unscaled parameters.
      print "Beginning input value smoothing"
      print "We are starting from unscaled parameters."
      self.scaleParsBy = self.getFlatVector(len(self.selectedbinxvals),1.0)
      print "start_vals_unscaled are:",start_vals_unscaled
      if debug :
        status = scipy.optimize.minimize(self.function, start_vals_unscaled, method='COBYLA', options=looseOpts)
      else :
        status = scipy.optimize.minimize(self.function, start_vals_unscaled, method='COBYLA', constraints=self.unscaledConstraints, options=looseOpts)
      start_vals_unscaled = status.x
    elif self.startValFormat == "flat":
      # If we are very, very far from a minimum, oversized steps may be the only way to get anywhere close.
      print "Beginning input value smoothing"
      if debug :
        status = scipy.optimize.minimize(self.function, start_vals, method='COBYLA', options=looseOpts)
      else :
        status = scipy.optimize.minimize(self.function, start_vals, method='COBYLA', constraints=self.myConstraints, options=looseOpts)
      start_vals = status.x
      start_vals_unscaled = numpy.multiply([Decimal(val) for val in start_vals],spectrum.scaleFactors)

    #updated_start_vals = numpy.divide([Decimal(val) for val in start_vals_unscaled],spectrum.scaleFactors)

    self.scaleParsBy = self.getFlatVector(len(self.selectedbinxvals),1.0)
    print "Beginning simple fit."
    print "start_val_unscaled are:",start_vals_unscaled
    if debug :
      status = scipy.optimize.minimize(self.function, start_vals_unscaled, method='SLSQP', jac=self.function_der, bounds=self.myBounds, options={'disp': True, 'maxiter':10000})
    else :
      status = scipy.optimize.minimize(self.function, start_vals_unscaled, method='SLSQP', jac=self.function_der, bounds=self.myBounds, constraints=self.unscaledConstraints, options={'disp': True, 'maxiter':10000})

    updated_start_vals = status.x
    print "unscaled, updated_start_vals are:",updated_start_vals
    updated_start_vals = numpy.divide([Decimal(val) for val in updated_start_vals],spectrum.scaleFactors)

    # Work on tightening convergence criteria!
    print "Returning to scaled parameters and beginning robust fit"
    options_dict = self.getOptionsDict(self.minAlg)
    self.scaleParsBy = spectrum.scaleFactors
    print self.minAlg,options_dict
    print "scaled, updated_start_vals are:",updated_start_vals
    if debug :
      status = scipy.optimize.minimize(self.function, updated_start_vals, method=self.minAlg, jac=self.function_der, bounds=self.myBounds, options=options_dict)
    else :
      status = scipy.optimize.minimize(self.function, updated_start_vals, method=self.minAlg, jac=self.function_der, bounds=self.myBounds, constraints=self.myConstraints, options=options_dict)
    print "Final parameter values:"
    print status.x
    self.parameterVals = status.x

    self.result = numpy.multiply([Decimal(val) for val in status.x],self.scaleParsBy)


    # Return a histogram with bin contents equal to fit results
    outputHist = spectrum.histogram.Clone("fitResult")
    outputHist.Reset()
    outputHist.SetDirectory(0)
    index = -1
    for bin in range(self.rangeLow,self.rangeHigh+1) :
      index = index+1
      outputHist.SetBinContent(bin,self.result[index]*self.selectedbinwidths[index])

    # If we want uncertainties on this, handle it now.
    if errType == "None" :
      print "Fit returned without uncertainties."
    elif errType == "Hessian" :
      print "Doing covariance-based errors!"
      hessMatrix = self.hessianLogL(self.selectedbincontents, self.result)
      covMatrix = self.covLogL(self.selectedbincontents,self.result)
      # Is it positive definite?
      # Raises exception if not.
      decomp = numpy.linalg.cholesky(covMatrix)
      # Made it here: we have a positive definite covariance matrix.
      # Uncertainy for parameter i is just sqrt(m[ii]).
      for bin in range(0,outputHist.GetNbinsX()) :
        error = math.sqrt(covMatrix[bin][bin])
        outputHist.SetBinError(bin+1,error)
    elif errType == "Bootstrap" :
      # Now we do a number of pseudoexperiments and take the variance of the parameter values
      # as the uncertainty in each bin.
      print "Beginning error bar computation"
      variances = self.getManyPEErrors(outputHist,options_dict)
      for bin in range(1,outputHist.GetNbinsX()+1) :
        if bin < self.rangeLow or bin > self.rangeHigh :
          outputHist.SetBinError(bin,0.0)
        else :
          error = variances[bin - self.rangeLow]
          print "Adding error",error,"to bin",bin
          outputHist.SetBinError(bin,error)

    else :
      print "Unknown error type specified!"
      print "Fit result will be returned without uncertainties."
    return outputHist
    
  def computeLogL(self,obs,exp) :

    # Likelihood for us is product of Poisson probability in each bin
    #   L = \Prod_{bins} (b^d * e^{-b})/d!
    # Log of this is:
    #   log L = \Sum_{bins} (log(b^2 * e^{-b}/d!))
    #         = \sum_{bins} (d ln b - b - ln(Gamma(d+1)))
    
    # Expectation in each bin is parameter times scale factor to bring it
    # to the proper magnitude.
    
    answer = Decimal(0.0)
    #print "pars:"
    for index in range(len(obs)) :

      if self.excludeWindow and index > self.windowLow-1 and index < self.windowHigh+1 :
        continue

      data = int(obs[index])
      bkg = Decimal(exp[index])*self.selectedbinwidths[index]*self.scaleParsBy[index]
      #print "\tcompare",data,"to",bkg
      if data < 0.0 or bkg < 0.0 :
        thisterm = Decimal(-1E10)
      elif data == 0.0 :
        thisterm = Decimal(-1.0)*bkg
      else :
#        thisterm = data * numpy.log(bkg) - bkg - scipy.special.gammaln(data+1.0)
        thisterm = data * bkg.ln() - bkg - Decimal(scipy.special.gammaln(data+1.0))
      #print "\t",exp[index]
      answer = answer - thisterm
        
    return answer


  def computeChi2(self,obs,exp) :

    return


  # Derivative w.r.t. PARAMETER
  # not w.r.t. background prediction
  # Parameter is what shifts.
  def jacobianLogL(self,obs,exp) :

    answer = []
    for index in range(len(obs)) :

      if self.excludeWindow and index > self.windowLow-1 and index < self.windowHigh+1 :
        answer.append(0.0)
        continue

      data = int(obs[index])
      bkg = Decimal(exp[index])*self.selectedbinwidths[index]*self.scaleParsBy[index]
      if data < 0.0 or bkg < 0.0:
        thisterm = Decimal(0.0)
      elif data == 0.0 :
        thisterm = self.selectedbinwidths[index]*self.scaleParsBy[index]
      else :
        thisterm = (Decimal(1.0) - (Decimal(data)/bkg))*self.selectedbinwidths[index]*self.scaleParsBy[index]
      answer.append(thisterm)

    return numpy.array(answer,dtype=Decimal)

    
  def jacobianChi2(self,obs,exp) :
  
    return

  # This will be used for a hacky stand-in for uncertainty calculation
  # in the case that we do not have time for a boostrap method
  # (read: all my fit stability studies)
  def hessianLogL(self,obs,exp) :

    answer = []
    for index in range(len(obs)) :

      thisrow = []
      for inIndex in range(len(exp)) :

        if index != inIndex :
          thisrow.append(0)
          continue

        data = int(obs[index])
        bkg = exp[index]*self.selectedbinwidths[index]*self.scaleParsBy[index]

        if (self.excludeWindow and inIndex > self.windowLow-1 and inIndex < self.windowHigh+1) or\
           (data <= 0.0) or (bkg < 0.0) :
          thisrow.append(0)
          continue

        else :
          thisrow.append(-data/(bkg*bkg))
      
      answer.append(thisrow)

    return numpy.array(answer)

  # Inverting a diagonal matrix is easy: just inverse of every non-zero element.
  # Covariance matrix is inverse of negative Hessian.
  def covLogL(self,obs,exp) :

    answer = []
    for index in range(len(obs)) :

      thisrow = []
      for inIndex in range(len(exp)) :

        if index != inIndex :
          thisrow.append(0)
          continue

        data = int(obs[index])
        bkg = exp[index]*self.selectedbinwidths[index]*self.scaleParsBy[index]

        if (self.excludeWindow and inIndex > self.windowLow-1 and inIndex < self.windowHigh+1) or\
           (data <= 0.0) or (bkg < 0.0) :
          thisrow.append(0)
          continue

        else :
          thisrow.append((bkg*bkg)/float(obs[index]))
      
      answer.append(thisrow)

    return numpy.array(answer)

  def getManyPEErrors(self,nominalHist,options_dict) :

    # Only generate this once so the seed keeps
    # all following ones independent
    nominalHistWrapper = WrappedHist(nominalHist)

    binArrays = []
    for bin in range(len(self.result)) :
      binArrays.append([])

    self.lastUncertaintyFitStatuses = []
    for PE in range(self.nPEs) :

      thisPE = nominalHistWrapper.poissonFluctuateBinByBin()
      PEWrapper = WrappedHist(thisPE)
      # Need to reset thing we run on to be contents of thisPE
      self.selectedbincontents, self.selectedbinxvals, self.selectedbinwidths, self.windowLow, self.windowHigh = PEWrapper.getSelectedBinInfo(self.rangeLow,self.rangeHigh,self.firstBinInWindow,self.lastBinInWindow)
      thisStatus = scipy.optimize.minimize(self.function, self.parameterVals, method=self.minAlg, jac=self.function_der, bounds=self.myBounds, constraints=self.myConstraints, options=options_dict)
      paramResults = thisStatus.x
      # Multiply parameter values by scales to get equivalent of function val
      binResults = numpy.multiply([Decimal(val) for val in paramResults],self.scaleParsBy)
      # ... and now by bin width to get equivalent of bin content.
      binResults = numpy.multiply(binResults,self.selectedbinwidths)
      nominalBinResults = numpy.multiply(self.result,self.selectedbinwidths)
      self.lastUncertaintyFitStatuses.append(thisStatus.status)
      index = -1
      for newval, originalval in zip(binResults,nominalBinResults) :
        index = index+1
        binArrays[index].append(newval - originalval)

    variances = []
    self.lastUncertaintyBinContents = []
    for bin in range(len(self.result)) :
      vec = numpy.array(binArrays[bin])
      self.lastUncertaintyBinContents.append(vec)

      mean = numpy.mean(vec)
      x2 = 0
      for input in vec :
        x2 = x2 + (input-mean)*(input-mean)
      variances.append(numpy.sqrt(x2/len(vec)))

    return variances


