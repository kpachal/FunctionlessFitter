import ROOT
import StatisticalTest
import numpy
from MathFunctions import poissonPVal, poissonConvGammaPVal
import HistWrapper

class BumpHunter(StatisticalTest.StatisticalTest) :

  def __init__(self) :
    StatisticalTest.StatisticalTest.__init__(self)
    self.allowDeficit = False
    self.useSidebands = False
    self.minBinsInBump = 2
    self.maxBinsInBump = 1e5
    self.nBinsInSideband = 1
    self.doErr = False
    self.tomography = None
    self.excludeWindow = False
    self.firstBinToExclude = -1
    self.lastBinToExclude = -1

  def doTest(self, dataHist, bkgHist, firstBinToUse, lastBinToUse) :

    dataCore = dataHist.histogram
    bkgCore = bkgHist.histogram
  
    assert dataCore.GetNbinsX() == bkgCore.GetNbinsX()
  
    # Find first and last bins with data
    # If reasonable, overwrite with user's choice
    firstBin = dataHist.firstBinWithData
    lastBin = dataHist.lastBinWithData
    if firstBinToUse>0 and firstBinToUse > firstBin and firstBinToUse < lastBin : firstBin = firstBinToUse
    if lastBinToUse > firstBinToUse and lastBinToUse>0 and lastBinToUse > firstBin and lastBinToUse < lastBin :
      lastBin = lastBinToUse
      
    regionsDef = []
    if self.excludeWindow :
      regionsDef.append([firstBin,self.firstBinToExclude-1])
      regionsDef.append([self.lastBinToExclude+1,lastBin])
    else :
      regionsDef.append([firstBin,lastBin])

    self.mostInterestingDict = {"binlow" : 0, "binhigh" : 0, "prob" : 1.0}

    for region in regionsDef :
      
      nBins = region[1] - region[0] + 1
      minWidth = max(self.minBinsInBump,1)
      maxWidth = min(self.maxBinsInBump,int(nBins/2.0))

      self.doCalculationCore(dataCore,bkgCore,minWidth,maxWidth,region[0],region[1])
      
    self.tomography = ROOT.TGraphErrors()
    index = -1
    for windowDict in self.bumpInfoList :
      index = index+1
      self.tomography.SetPoint(index,(windowDict["binlow"]+windowDict["binhigh"])/2.0,windowDict["prob"])
      self.tomography.SetPointError(index,(windowDict["binhigh"]-windowDict["binlow"])/2.0,0)
    
    if self.mostInterestingDict["prob"] == 0 :
      self.findBumpInCaseOfIncalculable(dataCore,bkgCore,firstBin,lastBin)

    return -numpy.log(self.mostInterestingDict["prob"])

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

        data, dataErr, bkg, bkgErr = self.getEffectiveBandContentsWithError(dataHist, bkgHist, windowLeft, windowRight)
        
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
          LSdata, LSdataErr, LSbkg, LSbkgErr = self.getEffectiveBandContentsWithError(dataHist, bkgHist, sidebandLeft, windowLeft - 1)
          RSdata, RSdataErr, RSbkg, RSbkgErr = self.getEffectiveBandContentsWithError(dataHist, bkgHist, windowRight + 1, sidebandRight)
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
        if probability < self.mostInterestingDict["prob"] :
          self.mostInterestingDict = windowDict


  def getEffectiveBandContentsWithError(self,data, bkg, firstBin, lastBin) :
  
#    print "With firstBin, lastBin =",firstBin,lastBin
#    data.Print("all")
#    bkg.Print("all")

    dataInt = dataErr = bkgInt = bkgErr = 0.0
    for bin in range(firstBin, lastBin+1) :
      dataInt = dataInt + data.GetBinContent(bin)
      dataErr = dataErr + data.GetBinError(bin)
      bkgInt = bkgInt + bkg.GetBinContent(bin)
      bkgErr = bkgErr + bkg.GetBinError(bin)
    
    return dataInt, dataErr, bkgInt, bkgErr

  def findBumpInCaseOfIncalculable(self, data, bkg, firstBin, lastBin) :

    lastWasInf = False
    allInfsConsecutive = True
    singlebinsinf = []
    for bin in range(firstBin, lastBin+1) :
      D = data.GetBinContent(bin)
      B = bkg.GetBinContent(bin)
      thisbinpval = poissonPVal(D,B)
      if thisbinpval==0 and D>B :
        if len(singlebinsinf)>0 and lastWasInf==False :
          allInfsConsecutive = False
        singlebinsinf.append(bin)
        lastWasInf = True
      else : lastWasInf = False
  
    if len(singlebinsinf) > 0 and allInfsConsecutive :
      windowDict = {"binlow" : singlebinsinf[0],\
                    "binhigh" : singlebinsinf[-1],
                    "prob" : 0.0}
      self.mostInterestingDict = windowDict

  def getFurtherInformation(self) :
    return self.mostInterestingDict["binlow"], self.mostInterestingDict["binhigh"], self.tomography

