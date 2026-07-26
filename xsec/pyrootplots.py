
#!/usr/bin/python

from ROOT import *
from ROOT import PlotUtils

try:
    import array as np
except:
    exit()
gROOT.SetBatch()  # Don't render histograms to a window.  Also gets filled areas correct.

TH1.AddDirectory(False)
variables = [
   "mixtpi",
#   "mixthetapi_deg",
#   "enu",
   "pmu",
    "ptmu",
    "pzmu",
    "q2",
    "thetamu_deg",
    "wexp"
]
# variables = ["mixtpi"]
date = "20241122"
warp = "WARP1"
plist = "ALL"

fdenom = TFile.Open(
        "DataXSecInputs_{DATE}_{PL}_mixed_CCPiSysAdded_sys_p4.root".format(
#        "DataXSecInputs_{DATE}_{PL}_mix_CCPiWeightSys_p4.root".format(
            PL=plist, DATE=date))

fnum = TFile.Open(
        "MCXSecInputs_20241217_ALL_mix_CCPiWeightNoSys_p4_{WARP}.root".format(
            PL=plist, DATE=date, WARP=warp))
fdata = TFile.Open(
#        "DataXSecInputs_{DATE}_{PL}_mix_CCPiWeightSys_p4.root".format(
        "DataXSecInputs_{DATE}_{PL}_mixed_CCPiSysAdded_sys_p4.root".format(
            PL=plist, DATE=date, WARP=warp))


for var in variables:
    lineWidth = 3

    #num = fnum.Get("effnum_{V}".format(V=var))
    #denom = fdenom.Get("effnum_{V}".format(V=var))
    #data = fdata.Get("bg_subbed_data_{V}".format(V=var))

    num = fdenom.Get("selection_mc_{V}".format(V=var))
    denom = fdenom.Get("selection_mc_untracked_no_tpi_weight_{V}".format(V=var))
    data = fdata.Get("selection_data_{V}".format(V=var))

    datapot = fdata.Get("data_pot").GetBinContent(1)
    mcpot = fdenom.Get("mc_pot").GetBinContent(1)
    potScale = datapot/mcpot
    nxbins = num.GetXaxis().GetNbins()
    denom.Scale(potScale)
    denom1 = denom.Clone()
    num.Scale(potScale)
    ratio = num.Clone()
    oneline = num.Clone()
    for i in range(nxbins):
        oneline.SetBinContent(i + 1, 1)

    titleName = ""
    units = ""
    if var == "mixtpi":
        titleName = "T_{#pi} (MeV)"
        units="(MeV)"
    if var == "pzmu":
        titleName = "p^{z}_{#mu} (MeV)"
        units="(MeV)"
    if var == "ptmu":
        titleName = "p^{t}_{#mu} (MeV)"
        units="(MeV)"
    if var == "thetapi_deg":
        titleName = "#theta_{#pi} (deg)"
        units="(deg)"
    if var == "pmu":
        titleName = "p_{#mu} (MeV)"
        units="(MeV)"
    if var == "q2":
        titleName = "Q^{2} (MeV^{2})"
        units="(MeV^{2})"
    if var == "thetamu_deg":
        titleName = "#theta_{#mu} (deg)"
        units="(deg)"
    if var == "wexp":
        titleName = "W_{exp} (MeV)"
        units="(MeV)"
    if var == "mixthetapi_deg":
        titleName = "#theta_{#pi} (deg)"
        units="(deg)"
    if var == "enu":
        titleName = "E_{#nu} (MeV)"
        units="(MeV)"

    warptitle = ""
  
    if warp == "WARP4":
      warptitle = "Warp = +20% M^{RES}_{A}"
    if warp == "WARP2":
      warptitle = "Warp = Anisotropic #Delta Decay"
    if warp == "WARP3":
      warptitle = "Warp = MK"
    if warp == "WARP5":
      warptitle = "Warp =T_{#pi} reweight"
    if warp == "NOMINAL":
      warptitle = "Closure test"
    c1 = TCanvas("Comparing T_{#pi} Weight aplication")
    Title = TPaveText(0, 0, 0.1, 0.1)
    Title.Clear()
    Title.AddText("Nominal vs " + warp + " " + titleName)
#    Title.AddText("Comparing T_{#pi} weight")
    Title.SetShadowColor(0)
    Title.SetLineColor(0)
    Title.SetFillColor(0)
    legend = TLegend(0.6, 0.75, 0.9, 0.9)
    c1.UseCurrentStyle()

    num.SetLineWidth(lineWidth)
    num.SetLineColor(kRed - 4)

    denom.SetLineWidth(lineWidth)
    denom.SetLineColor(kMagenta - 7)
    denom.SetTitle("Nominal vs " + warptitle)
    denom.GetXaxis().SetTitle(titleName)
    denom.GetXaxis().CenterTitle()
    denom.GetXaxis().SetTitleOffset(.9)
    denom.GetXaxis().SetTitleSize(0.04)

    denom.GetYaxis().SetTitle("Events/" + units)
    denom.GetYaxis().CenterTitle()
    denom.GetYaxis().SetTitleOffset(1.)
    denom.GetYaxis().SetTitleSize(0.05)
    data.SetMarkerStyle(20)
    data.SetLineWidth(1)
    data.SetMarkerColor(1)
    data.SetLineColor(1)
    denom.Scale(1., "width")
    num.Scale(1., "width")
    data.Scale(1., "width")
    ymax = denom.GetBinContent(denom.GetMaximumBin())*1.6
    denom.SetMaximum(ymax)
    denom.Draw("HIST E2")
    gStyle.SetOptStat(0)
    
    #c1.Update()

    num.Draw("SAME HIST E2")
    data.Draw("SAME")
    #c1.SetLogx()
    if var == "q2":
        c1.SetLogx()
    #c1.SetLogy()
    # c1.BuildLegend(0.7, 0.6, 0.9, 0.9)
    #legend.AddEntry(num, "T_{#pi} weight", "l")
    #legend.AddEntry(denom, " No T_{#pi} weight", "l")
    legend.AddEntry(num, warptitle, "l")
    legend.AddEntry(denom, "Nominal", "l")
    legend.AddEntry(data, "Data")
    legend.Draw()
    #Title.Draw()
    c1.Print(
        "Comparing_{PL}_newtpibinning_{DATE}_{VAR}_{WARP}.png".format(
            PL=plist, DATE=date, VAR=var, WARP=warp
        )
    )
    c1.Clear()
    
    c2 = TCanvas("Ratio Tpi weight")

    Title1 = TPaveText(0., 0., 0.1, 0.1)
#    Title1.Clear()
#    Title1.AddText("Ratio Nominal vs " + warp + " " + titleName)
    Title1.AddText("Ratio T_{#pi} weight")
    Title1.SetShadowColor(0)
    Title1.SetLineColor(0)
    Title1.SetFillColor(0)
    legend1 = TLegend(0.6, 0.75, 0.9, 0.9)
    ratiodata = data.Clone()
    ratiodata.SetMarkerStyle(20)
    ratiodata.SetLineWidth(1)
    ratiodata.SetMarkerColor(1)
    ratiodata.SetLineColor(1)
    ratio.Divide(ratio, denom1)
    ratiodata.Divide(ratiodata,denom)
    ratio.SetLineWidth(lineWidth)
    ratio.SetLineColor(kRed - 4)
    c2.UseCurrentStyle()

    if var == "q2":
        c2.SetLogx()

    oneline.SetTitle("Ratio Nominal vs " + warptitle )
    oneline.SetLineWidth(lineWidth)
    oneline.SetLineStyle(2)
    oneline.SetLineColor(kGray)

    oneline.GetXaxis().SetTitle(titleName)
    oneline.GetXaxis().CenterTitle()
    oneline.GetXaxis().SetTitleOffset(1.1)
    oneline.GetXaxis().SetTitleSize(0.04)

#    oneline.GetYaxis().SetTitle("T_{#pi}/NoT_{#pi}W")
    oneline.GetYaxis().CenterTitle()
    oneline.GetYaxis().SetTitleOffset(1.3)
    oneline.GetYaxis().SetTitleSize(0.04)

    ymax = ratiodata.GetBinContent(ratiodata.GetMaximumBin())*1.6
    oneline.SetMaximum(ymax)
    oneline.SetMinimum(0.0)

    oneline.Draw("HIST")
#    gStyle.SetOptStat(0)
    c2.SetTitle("Title")
    c2.Update()

    ratio.Draw("SAME HIST")
    ratiodata.Draw("SAME")
    #c1.SetLogx()
    #c1.SetLogy()
    # c1.BuildLegend(0.7, 0.6, 0.9, 0.9)
    #legend1.AddEntry(ratio, "T_{#pi}W/NoT_{#pi}W", "l")
    legend1.AddEntry(ratio, warptitle, "l")
    legend1.AddEntry(ratiodata, "data", "p")
    legend1.Draw()
    #Title1.Draw()
    c2.Print(
        "Ratio_{PL}_newtpibinning_{DATE}_{VAR}_{WARP}.png".format( PL=plist, DATE=date, VAR=var, WARP=warp
        )
    )

    c2.Clear()


