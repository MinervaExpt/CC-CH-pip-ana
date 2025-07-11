# TH1D * hErrtpi= dynamic_cast<TH1D*>(cross_section_mixtpi->GetVertErrorBand("GENIE_CCQEPauliSupViaKF")->GetErrorBand( true, false ).Clone("Test1"))

import math
from ROOT import *
from ROOT import PlotUtils

try:
    import array as np
except:
    exit()
gROOT.SetBatch()  # Don't render histograms to a window.  Also gets filled areas correct.


def round_to_sig_figs(num, sig_figs):
    if num == 0:
        return 0  # Edge case: If the number is zero, just return zero
    else:
        return round(num, sig_figs - int(math.floor(math.log10(abs(num)))) - 1)


TH1.AddDirectory(False)
variables = [
    "mixtpi",
    #   "mixthetapi_deg",
    "enu",
    "pmu",
    "ptmu",
    "pzmu",
    "q2",
    "thetamu_deg",
    "wexp",
]

File = TFile.Open("DataXSecInputs_20240622_ALL_mixed_newtpibinning_noSys_p4.root")
# File = TFile.Open("DataXSecInputs_20240610_ALL_thetapisig_NewEstimatorptmucut_Sys_p4.root")
Max = 0
Min = 999999
MaxVar = ""
MinVar = ""
h = File.Get("cross_section_q2")
study = "selection_mc"
ErrNamesVec = h.GetVertErrorBandNames()

for v in variables:
    h_cv_data = File.Get("cross_section_{VAR}".format(VAR=v))
    h_sel_data_tracked = File.Get("selection_data_tracked_{VAR}".format(VAR=v))
    h_sel_data_untracked = File.Get("selection_data_untracked_{VAR}".format(VAR=v))
    h_sel_data_mixed = File.Get("selection_data_mixed_{VAR}".format(VAR=v))
    h_sel_data = File.Get("selection_data_{VAR}".format(VAR=v))
    h_sel_mc = File.Get("selection_mc_{VAR}".format(VAR=v))
    #  h_cv_data.Scale(1.e42,"width")

    unfact = 0
    unitsCorr = 1
    if v == "mixtpi":
        unfact = 7.8
        h_cv_data.ModifyStatisticalUnc(unfact, "unfolding_cov_matrix_{V}".format(V=v))
    if v == "q2":
        unfact = 6.9
        h_cv_data.ModifyStatisticalUnc(unfact, "unfolding_cov_matrix_{V}".format(V=v))
        unitsCorr = 1000000
    if v == "mixthetapi_deg":
        unfact = 10.2
        h_cv_data.ModifyStatisticalUnc(unfact, "unfolding_cov_matrix_{V}".format(V=v))
    if v == "ptmu":
        unfact = 7.9
        h_cv_data.ModifyStatisticalUnc(unfact, "unfolding_cov_matrix_{V}".format(V=v))
        unitsCorr = 1000
    if v == "pmu":
        unitsCorr = 1000
    if v == "pzmu":
        unitsCorr = 1000
    if v == "enu":
        unitsCorr = 1000

    h_cv_mc = File.Get("mc_cross_section_{VAR}".format(VAR=v))
    h_cv_mc.Scale(1.0e42 * unitsCorr, "width")
    h_cv_data.Scale(1.0e42 * unitsCorr, "width")
    h_cv_data_flux = (
        h_cv_data.GetVertErrorBand("Flux").GetErrorBand(True, False).Clone("Flux")
    )
    h_total_error_abs = h_cv_data.GetTotalError(True, False)
    h_total_error_frac = h_cv_data.GetTotalError(True, True)
    h_sel_stat_frac = h_sel_data.GetStatError(True)
    h_sel_sys_frac = h_sel_mc.GetTotalError(False, True)
    h_staterr_abs = h_cv_data.GetStatError(False)
    h_staterr_frac = h_cv_data.GetStatError(True)

    nbins = h_cv_mc.GetNbinsX()
    n = 6
    ne = 5
    nf = 2
    factor = 1
    print("Cross section Table")
    print("Variable {VAR}".format(VAR=v))
    for i in range(1, nbins + 1):
        lowedgebin = h_cv_mc.GetBinLowEdge(i) / unitsCorr
        upedgebin = h_cv_mc.GetBinLowEdge(i + 1) / unitsCorr
        cv_data = round_to_sig_figs(h_cv_data.GetBinContent(i) * factor, n)
        cv_mc = round_to_sig_figs(h_cv_mc.GetBinContent(i) * factor, n)
        total_error_abs = round_to_sig_figs(
            h_total_error_abs.GetBinContent(i) * factor, ne
        )
        total_error_frac = round_to_sig_figs(h_total_error_frac.GetBinContent(i), nf)
        staterr_abs = round_to_sig_figs(h_staterr_abs.GetBinContent(i) * factor, ne)
        staterr_frac = round_to_sig_figs(h_staterr_frac.GetBinContent(i), nf)
        syserr_abs = round_to_sig_figs(
            math.sqrt(total_error_abs * total_error_abs - staterr_abs * staterr_abs), ne
        )
        flux_frac = round_to_sig_figs(h_cv_data_flux.GetBinContent(i), nf)

        syserr_frac = 0
        if cv_data != 0:
            syserr_frac = round_to_sig_figs(syserr_abs / cv_data, nf)
        # print(r"{LBE} - {UBE} & {CVD}e-02 & {MC}e-02 & {TEA}e-02 & {TEF} & {SEA}e-02 & {SEF} & {SYSEA}e-02 & {SYSEF}\\ \hline".format(
        print(
            r"{LBE} - {UBE} & {CVD} & {TEA} & {SEF} & {SYSEF} & {FLUXEF}\\".format(
                LBE=lowedgebin,
                UBE=upedgebin,
                CVD=cv_data,
                MC=cv_mc,
                TEA=total_error_abs,
                TEF=total_error_frac,
                SEA=staterr_abs,
                SEF=staterr_frac,
                SYSEA=syserr_abs,
                SYSEF=syserr_frac,
                FLUXEF=flux_frac,
            )
        )

    print("--------------------------------------------------------------------")

    print("Selection table")
    for i in range(1, nbins + 1):
        lowedgebin = h_cv_mc.GetBinLowEdge(i) / unitsCorr
        upedgebin = h_cv_mc.GetBinLowEdge(i + 1) / unitsCorr
        sel_data_total = round_to_sig_figs(h_sel_data.GetBinContent(i), n)
        sel_data_tracked = round_to_sig_figs(h_sel_data_tracked.GetBinContent(i), n)
        sel_data_untracked = round_to_sig_figs(h_sel_data_untracked.GetBinContent(i), n)
        sel_data_mixed = round_to_sig_figs(h_sel_data_mixed.GetBinContent(i), n)
        sel_stat_err = round_to_sig_figs(h_sel_stat_frac.GetBinContent(i), nf)
        sel_sys_err = round_to_sig_figs(h_sel_sys_frac.GetBinContent(i), nf)

        syserr_frac = 0
        if cv_data != 0:
            syserr_frac = round_to_sig_figs(syserr_abs / cv_data, nf)
        # print(r"{LBE} - {UBE} & {CVD}e-02 & {MC}e-02 & {TEA}e-02 & {TEF} & {SEA}e-02 & {SEF} & {SYSEA}e-02 & {SYSEF}\\ \hline".format(
        print(
            r"{LBE} - {UBE} & {SD} & {ST} & {SU} & {SM} & {SSTE} & {SSYE} \\".format(
                LBE=lowedgebin,
                UBE=upedgebin,
                SD=sel_data_total,
                ST=sel_data_tracked,
                SU=sel_data_untracked,
                SM=sel_data_mixed,
                SSTE=sel_stat_err,
                SSYE=sel_sys_err,
            )
        )

    print("--------------------------------------------------------------------")

print("DONE")
