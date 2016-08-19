import ROOT
import StatisticalTest

class BumpHunter(StatisticalTest) :

  def __init__(self) :
    self.allowDeficit = False
    self.useSidebands = False
    self.minBinsInBump = 1
    self.maxBinsInBump = 1e5
    self.nBinsInSideband = 1
    self.doErr = False
    self.lowEdge = 0.0
    self.highEdge = 0.0
    self.lowEdgesAllBumps = []
    self.highEdgesAllBumps = []
    self.probAllBumps = []

  def doTest(self, dataHist, bkgHist, firstBinToUse, lastBinToUse) :


  def doCalculationCore(self, dataHist, bkgHist, minWidth, maxWidth, firstBin, lastBin) :

    self.bumpInfoList = []
    self.mostInterestingDict = {"binlow" : 0, "binhigh" : 0, "prob" : 1.0}

    for width in range(minWidth,maxWidth+1) :

      # Sideband width needs to be something sensible
      if self.nBinsInSideband > 1 :
        sidebandWidth = self.nBinsInSideband
      else : sidebandWidth = max(1,int(width/2.0))
      
      smallestPForWidth = 1.0
      
      if self.useSidebands :
        minBinL = firstBin + sidebandWidth
        maxBinL = lastBin - width - sidebandWidth + 1
      else :
        minBinL = firstBin
        maxBinL = lastBin - width + 1

      # Loop over left edges possible with this bin width
      for windowLeft in range(minBinL, maxBinL+1) :

        # Other limits on window & sidebands
        windowRight = windowLeft + width - 1
        sidebandLeft = windowLeft - sidebandWidth
        sidebandRight = windowRight + sidebandWidth

        data, dataErr, bkg, bkgErr = getEffectiveBandContentsWithError(dataHist, bkgHist, windowLeft, windowRight)
        
        # Don't care about deficits unless otherwise specified
        if not self.allowDeficit and data < bkg :
          continue
          
        # Use uncertainty convolution if specified
        if self.doErr :
          probability = poissonConvGammaPVal(data, bkg, bkgErr)
        else :
          probability = poissonPVal(data, bkg)

        # Get probabilities for sidebands if desired.
        # If we have a big discrepancy in the sidebands then we do not
        # keep considering this window.
        if self.useSidebands :
          LSdata, LSdataErr, LSbkg, LSbkgErr = getEffectiveBandContentsWithError(dataHist, bkgHist, sidebandLeft, windowLeft - 1)
          RSdata, RSdataErr, RSbkg, RSbkgErr = getEffectiveBandContentsWithError(dataHist, bkgHist, windowRight + 1, sidebandRight)
          if self.doErr :
            probLeftSideband = poissonConvGammaPVal(LSdata, LSbkg, LSbkgErr)
            probRightSideband = poissonConvGammaPVal(RSdata, RSbkg, RSbkgErr)
          else :
            probLeftSideband = poissonPVal(LSdata, LSbkg, LSbkgErr)
            probRightSideband = poissonPVal(RSdata, RSbkg, RSbkgErr)

          if probLeftSideband < 1E-3 or probRightSideband < 1E-3 :
            continue
            
        # Save information on this window for the tomography plot
        windowDict = {"binlow" : windowLeft, "binhigh" : windowRight, "prob" : probability}
        self.bumpInfoList.append(windowDict)
        if probability < mostInterestingDict["prob"] :
          self.mostInterestingDict = windowDict


  def getEffectiveBandContentsWithError(self,data, bkg, firstBin, lastBin) :
  
    dataInt = dataErr = bkgInt = bkgErr = 0.0
    for bin in range(firstBin, lastBin+1) :
      dataInt = dataInt + data.GetBinContent(bin)
      dataErr = dataErr + data.GetBinError(bin)
      bkgInt = bkgInt + bkg.GetBinContent(bin)
      bkgErr = bkgErr + bkg.GetBinError(bin)
    
    return dataInt, dataErr, bkgInt, bkgErr

  def findBumpInCaseOfIncalculable(self, data, bkg, firstBin, lastBin) :

    lastWasInf = False
    allConsecutive = True
    singlebinsinf = []
    for bin in range(firstBin, lastBin+1) :
      D = data.GetBinContent(bin)
      B = bkg.GetBinContent(bin)
      thisbinpval = poissonPval(D,B)
      if thisbinpval==0 and D>B :
        if len(singlebinsinf)>0 and lastWasInf==False :
          allInfsConsecutive = False
        singlebinsinf.append(bin)
        lastWasInf = True
      else : lastWasInf = False
  
    if len(singlebinsinf) > 0 and allInfsConsecutive :
      windowDict = {"binlow" : data.GetBinLowEdge(singlebinsinf[0]),\
                    "binhigh" : data.GetBinLowEdge(singlebinsinf[-1]+1),
                    "prob" : 0.0}
      self.mostInterestingDict = windowDict

  def getFurtherInformation(self) :
    return self.bumpLowEdge, self.bumpHighEdge