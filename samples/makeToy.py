import ROOT
import numpy

infile = ROOT.TFile("dataLikeHistograms_IBLOff.2015.root","READ")
inhist = infile.Get("Nominal/mjj_Data_2015_3p57fb")

binEdgeVec = []
for bin in range(100,106) :
  binEdgeVec.append(inhist.GetBinLowEdge(bin))

toyHist = ROOT.TH1D("newHist_varyingBins","newHist_varyingBins",len(binEdgeVec)-1,numpy.array(binEdgeVec))
toyHist.SetDirectory(0)
toyHist2 = ROOT.TH1D("newHist_equalBins","newHist_equalBins",len(binEdgeVec)-1,binEdgeVec[0],binEdgeVec[-1])
toyHist2.SetDirectory(0)

index = 0
for bin in range(100,106) :
  index = index+1
  toyHist.SetBinContent(index,inhist.GetBinContent(bin))
  toyHist2.SetBinContent(index,inhist.GetBinContent(bin))

infile.Close()

nominalFile = ROOT.TFile("Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","READ")
nomhist = nominalFile.Get("basicBkgFrom4ParamFit")
nomhist.SetDirectory(0)
nomhist.SetName("boop")
nomFitResults = nominalFile.Get("fittedParameters")
nominalFile.Close()

nomHist = ROOT.TH1D("basicBkgFrom4ParamFit","basicBkgFrom4ParamFit",len(binEdgeVec)-1,numpy.array(binEdgeVec))
nomHist.SetDirectory(0)
nomHist2 = ROOT.TH1D("basicBkgFrom4ParamFit2","basicBkgFrom4ParamFit2",len(binEdgeVec)-1,binEdgeVec[0],binEdgeVec[-1])
nomHist2.SetDirectory(0)

index = 0
for bin in range(100,106) :
  index = index+1
  nomHist.SetBinContent(index,nomhist.GetBinContent(bin))
  nomHist2.SetBinContent(index,nomhist.GetBinContent(bin))

outfile = ROOT.TFile("toyTinyHists.root","RECREATE")
outfile.cd()
toyHist.Write()
toyHist2.Write()
nomHist.Write()
nomHist2.Write()
nomFitResults.Write("fittedParameters")
outfile.Close()
