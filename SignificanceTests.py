from MathFunctions import PoissonPVal, ProbToSigma, PoissonConvGammaPVal
import numpy
import scipy
import ROOT
from HistWrapper import WrappedHist

def getResidual(dataHist, bkgHist, firstBinToUse=-1, lastBinToUse=-1, errHist=None) :

  if type(dataHist) is WrappedHist :
    normData = dataHist.histogram
  else :
    normData = dataHist
    dataHist = WrappedHist(normData)
  if type(bkgHist) is WrappedHist :
    normBkg = bkgHist.histogram
  else :
    normBkg = bkgHist
    bkgHist = WrappedHist(normBkg)

  assert(normData.GetNbinsX()==normBkg.GetNbinsX())

  result = ROOT.TH1D(normBkg)
  resultname = "result_residual_{0}".format(normBkg.GetName())
  result.SetName(resultname)
  result.Reset()

  if firstBinToUse < 0 or firstBinToUse > normData.GetNbinsX()+1 :
    firstBin = dataHist.firstBinWithData
  else : firstBin = firstBinToUse
  if lastBinToUse < 0 or lastBinToUse > normData.GetNbinsX()+1 or lastBinToUse < firstBinToUse :
    lastBin = dataHist.lastBinWithData
  else : lastBin = lastBinToUse

  for bin in range(result.GetNbinsX()+2) :
    if bin>firstBin-1 and bin<lastBin+1 :

      D = normData.GetBinContent(bin)
      B = normBkg.GetBinContent(bin)
      bErr = normBkg.GetBinError(bin)

      if D != 0 : # Don't want residual Getting plotted outside of bins that actually contain something
        if (errHist) : PVal = PoissonConvGammaPVal(D,B,bErr)
        else : PVal = PoissonPVal(D,B)
        frac = ProbToSigma(PVal)
        # Trim to a reasonable size for display
        if (frac > 100) : frac = 20
        # use poissonPvalNonNegative
        if (frac<0.0) : frac = 0.0 #if it's negative (very insignificant), then set it to zero.
        if (D < B) : frac = frac * -1.0 #make it signed.
  
        result.SetBinContent(bin,frac)
        result.SetBinError(bin,0)

      # Data is zero: no significant result
      else :
        result.SetBinContent(bin,0)
        result.SetBinError(bin,0)

    # Outside our specified bin range: no significant result
    else :
      result.SetBinContent(bin,0)
      result.SetBinError(bin,0)

  return result
