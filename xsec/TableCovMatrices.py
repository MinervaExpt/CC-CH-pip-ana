# MakeCovarianceLatexTables.py

from collections import OrderedDict, namedtuple
import ROOT
import math
import sys, os

# LaTeX document structure
latex_prefix = [
    r"\documentclass[11pt]{article}",
    r"\usepackage{graphicx}",
    r"\usepackage{adjustbox}",
    r"\usepackage{lscape}",
    r"\usepackage{amsmath}",
    r"\usepackage{xspace}",
    r"\begin{document}",
    r"\setlength\tabcolsep{0.5em}",
    r"\include{minerva_preamble}",
]

latex_suffix = [r"\end{document}"]

# Data structures
varprops = namedtuple("varprops", "title unit binformat limits land")
covprops = namedtuple("covprops", "bins stat sys")

# Variable definitions from your script
variables = OrderedDict()
variables["mixtpi"] = varprops(r"$T_{\pi}$", r"\MeV", "{0:.0f}-{1:.0f}", [], False)
#variables["enu"] = varprops(r"$E_{\nu}$", r"GeV", "{0:.1f}-{1:.1f}", [], False)
variables["pmu"] = varprops(r"$p_{\mu}$", r"\GeVc", "{0:.1f}-{1:.1f}", [], False)
variables["ptmu"] = varprops(r"$p_{\mu,T}$", r"\GeVc", "{0:.2f}-{1:.2f}", [], True)
variables["pzmu"] = varprops(r"$p_{\mu,||}$", r"\GeVc", "{0:.0f}-{1:.0f}", [], False)
variables["q2"] = varprops(r"$Q^2$", r"\GeVsqcsq", "{0:.3f}-{1:.3f}", [], True)
variables["thetamu_deg"] = varprops(
    r"$\theta_{\mu}$", r"degree", "{0:.0f}-{1:.0f}", [], False
)
#variables["wexp"] = varprops(r"$W_{exp}$", r"GeV/c$^2$", "{0:.1f}-{1:.1f}", [], False)

# Units correction factors from your script
units_corrections = {
    "mixtpi": 1,
    "q2": 1e6,
    "ptmu": 1e3,
    "pmu": 1e3,
    "pzmu": 1e3,
    "enu": 1e3,
    "thetamu_deg": 1,
    "wexp": 1,
}

# Unfolding factors
unfolding_factors = {
    "mixtpi": 7.8,
    "q2": 6.9,
    "mixthetapi_deg": 10.2,
    "ptmu": 7.9,
}


def get_covariance_info(xsec_hist, var_name, var_info, units_corr):
    """Extract covariance matrix information from MnvH1D"""
    # Get covariance matrices
    stat_cov = xsec_hist.GetStatErrorMatrix()
    unfold_cov = xsec_hist.GetSysErrorMatrix("unfolding_cov_matrix_{0}".format(var_name))
    sys_cov = xsec_hist.GetTotalErrorMatrix(False)

    n_bins = xsec_hist.GetNbinsX()

    # Build bin edges
    bin_edges = []
    for i_bin in range(1, n_bins + 1):
        if len(var_info.limits) != 0:
            if i_bin < var_info.limits[0]:
                continue
            if var_info.limits[1] != -1 and i_bin > var_info.limits[1]:
                continue

        low = xsec_hist.GetBinLowEdge(i_bin) / units_corr
        high = xsec_hist.GetBinLowEdge(i_bin + 1) / units_corr
        bin_edges.append(var_info.binformat.format(low, high))

    # Build covariance matrices
    stat_cov_matrix = []
    sys_cov_matrix = []

    for i_bin_x in range(1, n_bins + 1):
        if len(var_info.limits) != 0:
            if i_bin_x < var_info.limits[0]:
                continue
            if var_info.limits[1] != -1 and i_bin_x > var_info.limits[1]:
                continue

        stat_cov_row = []
        sys_cov_row = []

        for i_bin_y in range(1, n_bins + 1):
            if len(var_info.limits) != 0:
                if i_bin_y < var_info.limits[0]:
                    continue
                if var_info.limits[1] != -1 and i_bin_y > var_info.limits[1]:
                    continue

            stat_val = (
                stat_cov[i_bin_x][i_bin_y] + unfold_cov[i_bin_x][i_bin_y]
            )
            sys_val = (sys_cov[i_bin_x][i_bin_y] - unfold_cov[i_bin_x][i_bin_y])

            stat_cov_row.append(stat_val)
            sys_cov_row.append(sys_val)

        stat_cov_matrix.append(stat_cov_row)
        sys_cov_matrix.append(sys_cov_row)

    return covprops(bin_edges, stat_cov_matrix, sys_cov_matrix)


def write_covariance_table(tex_table, caption, label, cov_info, var_info, cov_type):
    """Write a covariance matrix table in LaTeX format"""

    # Landscape mode if needed
    if var_info.land:
        tex_table.append(r"\pagebreak")
        tex_table.append(r"\begin{landscape}")

    # Table header
    tex_table.append(r"\begin{table}")
    tex_table.append(r"  \centering")
    tex_table.append(r"  \renewcommand{\arraystretch}{1.15}")
    tex_table.append(r"  \caption{" + caption + "}")

    if var_info.land:
        tex_table.append(r"  \begin{adjustbox}{max width=1.3\textheight}")
    else:
        tex_table.append(r"  \begin{adjustbox}{max width=\textwidth}")

    # Build column format
    colformat = "c||"
    header = ["Bin edges ({0})".format(var_info.unit)]
    for bin_range in cov_info.bins:
        colformat += "c"
        header.append("{0}".format(bin_range))

    tex_table.append(r"    \begin{tabular}{" + colformat + "}")
    tex_table.append("      " + " & ".join(header) + r"\\")
    tex_table.append(r"      \hline")
    tex_table.append(r"      \hline")

    # Matrix data
    if cov_type == "stat":
        matrix_data = cov_info.stat
    else:
        matrix_data = cov_info.sys

    for i, row in enumerate(matrix_data):
        line = cov_info.bins[i]
        for element in row:
            formatted = "{0:.3f}".format(element)
            if formatted.startswith("-0.") and float(formatted) == 0:
                formatted = formatted[1:]  # Remove the minus sign
            line += " & " + formatted
        if i != len(matrix_data) - 1:
            tex_table.append("      " + line + r"\\")
        else:
            tex_table.append("      " + line)

    # Table footer
    tex_table.append(r"    \end{tabular}")
    tex_table.append(r"  \end{adjustbox}")
    tex_table.append(r"  \label{" + label + "}")
    tex_table.append(r"\end{table}")

    if var_info.land:
        tex_table.append(r"\end{landscape}")
        tex_table.append(r"\clearpage")


def main():
    # Configuration
    input_file = "../rootfiles/DataXSecInputs_20241204_ALL_mix_CCPiWeightSys_p4.root"
    target_name = "scintillator"
    norm_type = "nucleon"

    # Open ROOT file
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        print("Error: Cannot open file", input_file)
        return

    # Storage for all tables
    all_tex_tables = []

    # Process each variable
    for var_name, var_info in variables.items():
        print("Processing variable:", var_name)

        # Get cross section histogram
        h_xsec = f.Get("cross_section_{0}".format(var_name))
        if not h_xsec:
            print("Warning: Cannot find histogram for", var_name)
            continue

        # Clone and prepare histogram
        xsec = h_xsec.Clone("xsec_{0}_scaled".format(var_name))

        # Apply unfolding factor if needed
        if var_name in unfolding_factors:
            xsec.ModifyStatisticalUnc(
                unfolding_factors[var_name], "unfolding_cov_matrix_{0}".format(var_name)
            )

        # Get units correction
        units_corr = units_corrections.get(var_name, 1)

        # Scale by units and bin width
        xsec.Scale(1e42 * units_corr, "width")

        # Get covariance information
        cov_info = get_covariance_info(xsec, var_name, var_info, units_corr)

        # Create statistical covariance table
        stat_caption = (
            "Statistical covariance matrix of measured cross section as function of {var}, "
            "in units of $10^{{-84}}$ $(\cmsq/{unit}/\\mathrm{{{norm}}})^2$"
        ).format(
            var=var_info.title, tar=target_name, unit=var_info.unit, norm=norm_type
        )
        stat_label = "tbl:statcov_{0}_{1}".format(var_name, target_name)

        stat_tex_table = []
        write_covariance_table(
            stat_tex_table, stat_caption, stat_label, cov_info, var_info, "stat"
        )
        all_tex_tables.extend(stat_tex_table)
        all_tex_tables.append("")  # Empty line between tables

        # Create systematic covariance table
        sys_caption = (
            "Systematic covariance matrix of measured cross section as function of {var}, "
            "in units of $10^{{-84}}$ $(\cmsq/{unit}/\\mathrm{{{norm}}})^2$"
        ).format(
            var=var_info.title, tar=target_name, unit=var_info.unit, norm=norm_type
        )
        sys_label = "tbl:syscov_{0}_{1}".format(var_name, target_name)

        sys_tex_table = []
        write_covariance_table(
            sys_tex_table, sys_caption, sys_label, cov_info, var_info, "sys"
        )
        all_tex_tables.extend(sys_tex_table)
        all_tex_tables.append("")  # Empty line between tables

        # Clean up
        del xsec

    # Write output file
    with open("CovarianceTables.tex", "w") as output:
        for line in latex_prefix:
            output.write(line + "\n")

        for line in all_tex_tables:
            output.write(line + "\n")

        for line in latex_suffix:
            output.write(line + "\n")

    print("LaTeX tables written to CovarianceTables.tex")

    # Close ROOT file
    f.Close()


if __name__ == "__main__":
    main()
