import ROOT
import numpy
from sympy import *

valuesList = [
{"inname":"Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","outname":"newNominal_EOYE.root","vallow": 1100.00, "valhigh": 7052.00 }
]
for valdict in valuesList :

  outname = valdict["outname"]
  infile = ROOT.TFile(valdict["inname"],"READ")
  inhist = infile.Get("basicBkgFrom4ParamFit")
  inhist.SetName("nominal_original")
  inhist.SetDirectory(0)

  parvec = infile.Get("fittedParameters")
  par0 = parvec[0]
  par1 = parvec[1]
  par2 = parvec[2]
  x = Symbol('x')
  y = par0*((1-x/13000.0)**par1)*((x/13000.0)**par2)
  funcstring = "{0}".format(y).replace("**","^")
  myFunc = ROOT.TF1("function",funcstring, valdict["vallow"], valdict["valhigh"])

  infile.Close()

  remakeNominal = inhist.Clone()
  remakeNominal.SetName("nominal_recreated")
  remakeNominal.SetDirectory(0)
  remakeNominal.Reset()

  newNominal = inhist.Clone()
  newNominal.SetName("nominal_simple")
  newNominal.SetDirectory(0)
  newNominal.Reset()

  for bin in range(1,inhist.GetNbinsX()+1) :

    if inhist.GetBinContent(bin) == 0 :
      remakeNominal.SetBinContent(bin,0)
      remakeNominal.SetBinError(bin,0)
      newNominal.SetBinContent(bin,0)
      newNominal.SetBinError(bin,0)
      continue

    x1 = inhist.GetBinLowEdge(bin)
    x2 = x1 + inhist.GetBinWidth(bin)
    simple = myFunc.Eval(inhist.GetBinCenter(bin))*inhist.GetBinWidth(bin)
    complex = myFunc.Integral(x1,x2)
    newNominal.SetBinContent(bin,simple)
    newNominal.SetBinError(bin,sqrt(simple))
    remakeNominal.SetBinContent(bin,complex)
    remakeNominal.SetBinError(bin,sqrt(complex))

  outfile = ROOT.TFile(outname,"RECREATE")
  outfile.cd()
  inhist.Write()
  remakeNominal.Write()
  newNominal.Write()
  parvec.Write("fittedParameters")
  outfile.Close()
