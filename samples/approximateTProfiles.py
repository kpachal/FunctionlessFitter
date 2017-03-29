import ROOT
import numpy
from sympy import *
from array import array

valuesList = [
{"inname":"Step1_SearchPhase_mjj_Data_2015_3p57fb_fluctOnData.root","outname":"tProfile_EOYE.root","vallow": 1100.00, "valhigh": 7052.00 },
{"inname":"SearchResultData_UA2_fullDataset_from443.root","outname":"tProfile_TLA.root","vallow":444,"valhigh":1260},
{"inname":"Step1_SearchPhase_mjj_Data_2016_15p7fb.root","outname":"tProfile_ICHEP.root","vallow":1100,"valhigh":7000}
]

for valdict in valuesList :

  outname = valdict["outname"]
  print "Beginning",outname
  infile = ROOT.TFile(valdict["inname"],"READ")
  inhist = infile.Get("basicBkgFrom4ParamFit")
  inhist.SetName("nominal_original")
  inhist.SetDirectory(0)

  profile = ROOT.TProfile("data_profile","data_profile",inhist.GetNbinsX(),array('d',inhist.GetXaxis().GetXbins()))

  myFunc = infile.Get("theFitFunction")

  for bin in range(inhist.FindBin(valdict["vallow"]),inhist.FindBin(valdict["valhigh"])+1) : # 1,inhist.GetNbinsX()+1

      xlow = inhist.GetBinLowEdge(bin)
      xhigh = inhist.GetBinLowEdge(bin+1)
      integral = myFunc.Integral(xlow,xhigh)
      #print "Func at low (",xlow,") and high (",xhigh,") edges of bin",bin,"is",myFunc.Eval(xlow),myFunc.Eval(xhigh),". Integral",integral
      for x in range(300,600) :
        #print x, x*(xhigh-xlow)/100.0
        xval = x*(xhigh-xlow)/1000.0 + xlow
        thisint = myFunc.Integral(xlow,xval)
        #print xval, thisint
        if not thisint > 0 :
          xval = inhist.GetBinCenter(bin)
          break
        #print "testing xval",xval,": integral",xlow,xval,"=",thisint,". Compare to",0.5*integral
        if thisint > 0.5*integral :
          break
#      print "for bin",bin,"kept xval",xval,"compared to nominal",inhist.GetBinCenter(bin), "( diff", inhist.GetBinCenter(bin) - xval, ")"
      profile.Fill(xval,xval,inhist.GetBinContent(bin))

  outfile = ROOT.TFile(outname,"RECREATE")
  outfile.cd()
  inhist.Write()
  profile.Write()
  outfile.Close()
