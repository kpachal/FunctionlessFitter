import ROOT

class StatisticalTest(object) :

  def __init__(self) :
    self.excludeWindow = False
    self.firstBinToExclude = 0
    self.lastBinToExclude = 0

  def setWindowExclusion(self, lowBin, highBin) :
    self.firstBinToExclude = lowBin
    self.lastBinToExclude = highBin

  def doTest(self) :
    return 0

  def getFurtherInformation(self) :
    return []