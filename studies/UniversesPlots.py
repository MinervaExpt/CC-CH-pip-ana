

from ROOT import *
from ROOT import PlotUtils

try:
    import array as np
except:
    exit()
gROOT.SetBatch()  # Don't render histograms to a window.  Also gets filled areas correct.

TH1.AddDirectory(False)
variables = [
#    "mixtpi",
#    "mixthetapi_deg",
#    "enu",
#    "pmu",
#    "ptmu",
#    "pzmu",
#    "q2",
#    "thetamu_deg",
    "mixtpi_true",
    "mixthetapi_deg_true",
    "enu_true",
    "pmu_true",
    "ptmu_true",
    "pzmu_true",
    "q2_true",
    "thetamu_deg_true",
]

#File = TFile.Open("MCXSecInputs_1110_ME1A_0_2024-12-02.root")
File = TFile.Open("DataXSecInputs_20241204_ALL_mix_CCPiCVandSys_p4.root")
Max = 0
Min = 999999 
MaxVar = ""
MinVar = "" 
h = File.Get("selection_mc_q2")
study = "cross_section"
#study = "selection_mc"
#study = "effnum"
#ErrNamesVec = h.GetVertErrorBandNames() 
ErrNamesVec = ["CCPi+ Tune"]
print (ErrNamesVec)
for err in ErrNamesVec: 
  print (err)
  for v in variables:
    errname = ""
    if err == "CCPi+ Tune":
        errname = "CCPi+_Tune"
    else:
        errname = err
    c1 = TCanvas("Warping studies for")
    step = "{STU}_{VAR}".format(STU=study,VAR=v)
    title = "{STU}_{VAR}_{ERR}".format(STU=study,VAR=v,ERR=err)
    hist = File.Get(step)
    hErr=hist.GetVertErrorBand(err).GetHist(0)
    
    Max=hErr.GetBinContent(hErr.GetMaximumBin())*1.4
    Min=hErr.GetBinContent(hErr.GetMinimumBin())
    
    hErr.GetYaxis().SetRangeUser(0,Max)
    gStyle.SetOptStat(0)


    hErr.Draw("HIST")
    hist.Draw("SAME HIST")
    legend = TLegend(0.7, 0.75, 0.9, 0.9)
    legend.AddEntry(hist,"CV","l")
    legend.AddEntry(hErr,err,"l")
    legend.Draw()
    c1.Print("Universes_CV_{ERR}_{VAR}_{ST}.png".format(ERR=errname, VAR=v, ST=study))
     

print("DONE")
