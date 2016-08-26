import ROOT
import StatisticalTest
import numpy
import scipy

class Chi2Test(StatisticalTest.StatisticalTest) :

  def __init__(self) :
    StatisticalTest.StatisticalTest.__init__(self)

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

    answer = 0.0
    for bin in range (firstBin, lastBin+1) :

      if self.excludeWindow :
        if (bin > self.firstBinToExclude - 1 and bin < self.lastBinToExclude+1) : continue

      data = dataCore.GetBinContent(bin)
      if (data==0) : continue
      bkg = bkgCore.GetBinContent(bin)
      deltaB = bkgCore.GetBinError(bin)

      # Definitions of chi2 itself
      # term = (d - b) / sqrt(b) //Pearson's Chi is defined with sqrt(b), and only for populated bins.
      term = (data - bkg) / numpy.sqrt(bkg+deltaB*deltaB) #give it an extra uncertainty in the background if it's there
      answer = answer + (term*term)

    return answer
