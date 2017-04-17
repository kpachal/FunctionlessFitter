import ROOT
import numpy
import StatisticalTest
import FunctionlessFitter
from HistWrapper import Dataset
from MathFunctions import makeHistFromVector, getPValFromVecAndStat

class PseudoExperimenter(object) :

  def __init__(self) :
    self.fitSuccessRate = 0.0

  def getPseudoexperiments(self,dataHist,bkgHist,statTests,firstBinToUse,lastBinToUse,nPEs=1000,fitter=None) :

    if not isinstance(statTests, (list, tuple)) :
      statTests = [statTests]

    results = {}

    # Get statistics for true data
    for statTest in statTests :
    
      results[statTests.index(statTest)] = {}
    
      originalStat = statTest.doTest(dataHist, bkgHist, firstBinToUse, lastBinToUse)
      originalFurtherInformation = statTest.getFurtherInformation()

      results[statTests.index(statTest)]["stat"] = originalStat
      results[statTests.index(statTest)]["furtherInformation"] = originalFurtherInformation
      results[statTests.index(statTest)]["PEStats"] = []
      results[statTests.index(statTest)]["PEFurtherInformation"] = []

    # Calculate each stat for each PE
    PE = -1
    while PE < nPEs+1 :
      PE = PE + 1
    
      if PE%50 == 0 : print "on PE",PE

      PEHist = Dataset(bkgHist.poissonFluctuateBinByBin(),binSpecifier=bkgHist.binxvals,baseName="toy{0}".format(PE))

      # Default: we throw PEs and compare directly to the nominal background
      # prediction rather than refitting
      if not fitter :
        useBkg = bkgHist
      else :
        result = fitter.fit(PEHist,firstBinToUse,lastBinToUse)
        result.SetDirectory(0)
        useBkg = Dataset(result)
      
      statVector = []
      FUVector = []
      skip = False
      for statTest in statTests :

         thisStat = statTest.doTest(PEHist, useBkg, firstBinToUse, lastBinToUse)
         thisFurtherInformation = statTest.getFurtherInformation()
         
         # Protect against unphysical values
         if numpy.isinf(thisStat) or numpy.isnan(thisStat) : skip = True
         
         # If physical, save it
         statVector.append(thisStat)
         FUVector.append(thisFurtherInformation)
        
      if skip :
        PE = PE - 1
        continue
        
      # Now we know all values are physical, store this PE.
      for statTest in statTests :
        results[statTests.index(statTest)]["PEStats"].append(statVector[statTests.index(statTest)])
        results[statTests.index(statTest)]["PEFurtherInformation"].append(FUVector[statTests.index(statTest)])

    # Work out a few more quantities we'll want.
    for statTest in statTests :
    
      pVal = getPValFromVecAndStat(results[statTests.index(statTest)]["stat"],results[statTests.index(statTest)]["PEStats"])
      results[statTests.index(statTest)]["pValue"] = pVal
    
      histogram = makeHistFromVector(results[statTests.index(statTest)]["PEStats"])
      histogram.SetName(histogram.GetName()+"_{0}".format(statTests.index(statTest)))
      results[statTests.index(statTest)]["statHist"] = histogram


    # Now return a dictionary of results
    return results

