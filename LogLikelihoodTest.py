import ROOT
import numpy
import scipy
import StatisticalTest

class LogLikelihoodTest(StatisticalTest.StatisticalTest) :

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

      if data < 0.0 or bkg < 0.0 :
        thisterm = -1E10
      elif data == 0.0 :
        thisterm = -1.0*bkg
      else :
        thisterm = data * numpy.log(bkg) - bkg - scipy.special.gammaln(data+1)
      answer = answer + thisterm

    return answer