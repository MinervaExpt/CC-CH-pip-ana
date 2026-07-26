#ifndef plotNuisance_C
#define plotNuisance_C

#include <cassert>  // !! must compile in debug mode root script.C++g
#include <iostream>
#include <sstream>

#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif

#include "Constants.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MnvColors.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvVertErrorBand.h"
#include "Plotter.h"
#include "SignalDefinition.h"
#include "Systematics.h"
#include "TArrayD.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TCollection.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TText.h"
#include "Variable.h"
#include "includes/MacroUtil.h"
#include "includes/Plotter.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "includes/myPlotStyle.h"
#include "xsec/makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "xsec/plotting_functions.h"        // RebinQ2Plot

void plot_all_models(Plotter p, MnvH1D* data,
                     std::map<std::string, PlotUtils::MnvH1D*> mc,
                     PlotUtils::MnvH1D* Mnv_mc_xsec, std::string outdir = ".",
                     double ymax = -1, bool do_log_scale = false,
                     bool do_bin_width_norm = true, bool plotRatios = false) {
  std::cout << "Plotting CrossSection " << p.m_variable->Name() << "\n";

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  for (const auto& i : mc) assert(i.second);
  //  p.m_mnv_plotter.ApplyStyle(PlotUtils::kCCCohStyle);

  // Clone the input hists. Don't modify originals.
  PlotUtils::MnvH1D* data_xsec = nullptr;
  PlotUtils::MnvH1D* mnv_mc_xsec = nullptr;
  std::map<std::string, PlotUtils::MnvH1D*> mc_xsec;
  std::map<std::string, PlotUtils::MnvH1D*> ratio_mc_xsec;

  double uncfactor;
  std::string name = p.m_variable->Name();
  if (name == "mixtpi") {
    // uncfactor = 5.7;// Without CCPi weight
    uncfactor = 6.5;  // Wight CCPi Weight
    data->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "mixthetapi_deg") {
    // uncfactor = 10.2;// Without CCPi weight
    uncfactor = 6.2;  // Wight CCPi Weight
    data->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "q2") {
    // uncfactor = 6.9;// Without CCPi weight
    uncfactor = 7.5;  // Wight CCPi Weight
    data->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "ptmu") {
    // uncfactor = 7.9;// Without CCPi weight
    uncfactor = 7.3;  // Wight CCPi Weight
    data->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "thetamu_deg") {
    uncfactor = 15.2;  // Wight CCPi Weight
    data->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  p.m_mnv_plotter.mc_line_width = 3;
  p.m_mnv_plotter.good_colors = {1,
                                 1,
                                 TColor::GetColor("#E69F00"),
                                 TColor::GetColor("#56B4E9"),
                                 TColor::GetColor("#009E73"),
                                 TColor::GetColor("#0072B2"),
                                 TColor::GetColor("#D55E00"),
                                 TColor::GetColor("#CC79A7")};
  // rebin q2
  const bool do_q2_rebin = true;
  const bool do_rebin_GeV = true;
  if (do_q2_rebin and p.m_variable->Name() == "q2") {
    data_xsec = RebinQ2Plot(*data);
    mnv_mc_xsec = RebinQ2Plot(*Mnv_mc_xsec);
    for (auto i : mc) {
      mc_xsec[i.first] = RebinQ2Plot(*i.second);
    }
  } else if (do_q2_rebin and
             (p.m_variable->Name() == "pmu" || p.m_variable->Name() == "enu" ||
              p.m_variable->Name() == "ptmu" ||
              p.m_variable->Name() == "pzmu")) {
    data_xsec = RebinningtoGeV(*data, p.m_variable->Name());
    mnv_mc_xsec = RebinningtoGeV(*Mnv_mc_xsec, p.m_variable->Name());
    std::cout << "data_xsec->GetXaxis()->GetNbins() = "
              << data_xsec->GetXaxis()->GetNbins() << "\n";
    for (auto i : mc)
      mc_xsec[i.first] = RebinningtoGeV(*i.second, p.m_variable->Name());
  } else {
    data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
    mnv_mc_xsec = (PlotUtils::MnvH1D*)Mnv_mc_xsec->Clone("Mnv_mc_sec");
    for (const auto& i : mc)
      mc_xsec[i.first] = (PlotUtils::MnvH1D*)i.second->Clone(i.first.c_str());
  }

  PlotUtils::MnvH1D* h_flux = nullptr;
  PlotUtils::MnvH1D* h_flux_integrated = nullptr;
  PlotUtils::MnvH1D* flux = nullptr;

  if (p.m_variable->Name() == "enu") {
    h_flux = (PlotUtils::MnvH1D*)data_xsec->Clone("enu_clone");
    h_flux->Reset();

    // Get the flux histo, to be integrated
    static PlotUtils::FluxReweighter* frw = new PlotUtils::FluxReweighter(
        14, CCNuPionIncConsts::kUseNueConstraint, "minervame1D1M1NWeightedAve",
        PlotUtils::FluxReweighter::gen2thin,
        PlotUtils::FluxReweighter::g4numiv6,
        CCNuPionIncConsts::kNFluxUniverses);

    h_flux_integrated = frw->GetIntegratedFluxReweighted(14, h_flux, 0., 100.);
    h_flux_integrated->Scale(1.0e-4);

    h_flux = UndoBWN(frw->GetRebinnedFluxReweighted(14, h_flux));
    h_flux->Scale(1.0e-4);

    TH1* h_flux_aux = (TH1*)h_flux->Clone("TH1Flux");

    for (int i = 1; i <= h_flux->GetNbinsX(); i++) {
      std::cout << "Bin = " << i << " Low Edge = " << h_flux->GetBinLowEdge(i)
                << "  " << h_flux->GetBinContent(i) << "  "
                << h_flux_integrated->GetBinContent(i) << "\n";
    }

    data_xsec->Scale(h_flux_integrated->GetBinContent(1));
    data_xsec->DivideSingle(data_xsec, h_flux);
    for (const auto& i : mc) {
      mc_xsec[i.first]->Scale(h_flux_integrated->GetBinContent(1));
      mc_xsec[i.first]->DivideSingle(mc_xsec[i.first], h_flux);
    }
  }

  data_xsec->SetTitle("Data");
  TCanvas canvas("c1", "c1");

  // Get Hists
  TH1D* data_xsec_w_tot_error = new TH1D(data_xsec->GetCVHistoWithError());
  TH1D* data_xsec_w_stat_error = new TH1D(data_xsec->GetCVHistoWithStatError());
  TH1* num_data = (TH1*)data_xsec_w_tot_error->Clone("num_Data");
  TH1* num_data_stat_error =
      (TH1*)data_xsec_w_stat_error->Clone("num_Data_stat_error");
  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
  if (p.m_variable->Name() == "q2") {
    canvas.SetLogx();
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;
  if (p.m_variable->Name() == "mixthetapi_deg")
    p.m_mnv_plotter.axis_maximum = 48;
  if (p.m_variable->Name() == "q2") p.m_mnv_plotter.axis_maximum = 8.5;
  if (p.m_variable->Name() == "ptmu") p.m_mnv_plotter.axis_maximum = 6.5;

  // Y-label
  p.m_mnv_plotter.axis_title_offset_y = 0.;

  // X label
  p.SetXLabel(data_xsec_w_tot_error);
  p.SetXLabel(data_xsec_w_stat_error);
  p.SetXLabel(data_xsec);
  p.SetXLabel(num_data);
  p.SetXLabel(mnv_mc_xsec);
  for (auto i : mc_xsec) p.SetXLabel(i.second);

  for (auto i : mc_xsec) {
    ratio_mc_xsec[i.first] =
        (PlotUtils::MnvH1D*)i.second->Clone(i.first.c_str());
  }

  for (auto i : ratio_mc_xsec) p.SetXLabel(i.second);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = 1.;
  }
  double yscale = 1.e42;
  std::string stryscale = "-42";
  // Bin Width Normalization, Y-axis label, and 10^-42 shift
  if (do_bin_width_norm) {
    if (p.m_variable->Name() == "q2" || p.m_variable->Name() == "enu" ||
        p.m_variable->Name() == "ptmu") {
      yscale = 1.e39;
      stryscale = "-39";
    }
    data_xsec_w_tot_error->Scale(yscale, "width");
    data_xsec_w_stat_error->Scale(yscale, "width");
    data_xsec->Scale(yscale, "width");
    for (auto i : mc_xsec) i.second->Scale(yscale, "width");
    // Fix q2 first bin
    if (do_q2_rebin and p.m_variable->Name() == "q2") {
      // already divided by 0.25 - 0.006
      // but we really want to divide by 0.25
      // TODO sketchy AF
      double scale = (0.025 - 0.006) / 0.025;

      data_xsec_w_tot_error->SetBinContent(
          2, data_xsec_w_tot_error->GetBinContent(2) * scale);
      data_xsec_w_tot_error->SetBinError(
          2, data_xsec_w_tot_error->GetBinError(2) * scale);

      data_xsec_w_stat_error->SetBinContent(
          2, data_xsec_w_stat_error->GetBinContent(2) * scale);
      data_xsec_w_stat_error->SetBinError(
          2, data_xsec_w_stat_error->GetBinError(2) * scale);

      data_xsec->SetBinContent(2, data_xsec->GetBinContent(2) * scale);
      data_xsec->SetBinError(2, data_xsec->GetBinError(2) * scale);

      for (auto i : mc_xsec) {
        i.second->SetBinContent(2, i.second->GetBinContent(2) * scale);
        i.second->SetBinError(2, i.second->GetBinError(2) * scale);
      }
    }

    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel +
                        Form(" (10^{%s} cm^{2}/", stryscale.c_str()) +
                        p.m_variable->m_units + "/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    for (auto i : mc_xsec) i.second->GetYaxis()->SetTitle(yaxis.c_str());
  } else {
    if (p.m_variable->Name() == "q2" || p.m_variable->Name() == "enu" ||
        p.m_variable->Name() == "ptmu") {
      yscale = 1.e39;
      stryscale = "-39";
    }
    data_xsec_w_tot_error->Scale(yscale);
    data_xsec_w_stat_error->Scale(yscale);
    data_xsec->Scale(yscale);
    for (auto i : mc_xsec) i.second->Scale(yscale);
    // Fix q2 first bin
    if (do_q2_rebin and p.m_variable->Name() == "q2") {
      // already divided by 0.25 - 0.006
      // but we really want to divide by 0.25
      // TODO sketchy AF
      double scale = (0.025 - 0.006) / 0.025;

      data_xsec_w_tot_error->SetBinContent(
          2, data_xsec_w_tot_error->GetBinContent(2) * scale);
      data_xsec_w_tot_error->SetBinError(
          2, data_xsec_w_tot_error->GetBinError(2) * scale);

      data_xsec_w_stat_error->SetBinContent(
          2, data_xsec_w_stat_error->GetBinContent(2) * scale);
      data_xsec_w_stat_error->SetBinError(
          2, data_xsec_w_stat_error->GetBinError(2) * scale);

      data_xsec->SetBinContent(2, data_xsec->GetBinContent(2) * scale);
      data_xsec->SetBinError(2, data_xsec->GetBinError(2) * scale);

      for (auto i : mc_xsec) {
        i.second->SetBinContent(2, i.second->GetBinContent(2) * scale);
        i.second->SetBinError(2, i.second->GetBinError(2) * scale);
        i.second->SetLineColor(8);
      }
    }
    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis =
        Form("#sigma (10^{%s} cm^{2}/nucleon)", stryscale.c_str());
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    for (auto i : mc_xsec) i.second->GetYaxis()->SetTitle(yaxis.c_str());
  }
  // Create a TObjArray from the input vector of model xsecs
  TObjArray* mc_xsec_arr = new TObjArray;
  mc_xsec_arr->Add(mc_xsec["MnvTune 4.7.1"]);
  for (const auto& i : mc_xsec) {
    if (i.first == "MnvTune 4.7.1") continue;
    mc_xsec_arr->Add(i.second);
  }
  gStyle->SetEndErrorSize(4);
  // Plot it
  const bool use_hist_titles = true;
  p.m_mnv_plotter.DrawDataMCVariations(data_xsec, mc_xsec_arr, pot_scale, "TR",
                                       use_hist_titles);
  data_xsec_w_stat_error->SetMarkerStyle(20);
  data_xsec_w_stat_error->SetMarkerSize(1);
  data_xsec_w_stat_error->SetLineColor(1);
  data_xsec_w_stat_error->SetLineWidth(3);
  data_xsec_w_stat_error->Draw("X0 E1 SAME");

  std::vector<PlotUtils::MnvH1D*> ln_mc_xsec_arr;
  for (const auto& i : mc_xsec) {
    //    PlotUtils::MnvH1D* mc_xsec_potnorm = i.second->Clone(Form("POTNorm%s",
    //    i.first.c_str())); mc_xsec_potnorm
    ln_mc_xsec_arr.push_back(LogHist(*i.second, p.m_variable->Name()));
  }

  // Add chi2 label
  {
    const bool use_data_error_mtx = true;
    const bool use_only_shape_errors = false;
    const bool use_model_stat =
        false;  // model statistical errors -- these are small if you use a
                // priori effden, but atm I'm not. So keep this off.
    // p.m_mnv_plotter.AddChi2Label(data_xsec, mc_xsec, pot_scale, "TR", 0.03,
    // -0.175, use_data_error_mtx, use_only_shape_errors); // this auto turns on
    // model stat errors, I've manually turned it off within the function.

    // std::cout << "use_overflow is " << p.m_mnv_plotter.chi2_use_overflow_err
    // << "\n";
    // Check chi2
    int ndf = -1;
    int nof = -1;
    // double pot_scale = 1.;
    std::cout << p.m_variable->Name() << "\n";
    std::cout
        << "Model & Conventional $\\chi^2$ & Log-normal $\\chi^2$ \\\\ \n";
    PlotUtils::MnvH1D* data_log = LogHist(*data_xsec, p.m_variable->Name());
    int auxnbins = data_xsec->GetXaxis()->GetNbins();
    std::cout << "Numero de bines = " << auxnbins << "\n";
    std::cout << "Fist bin = " << data_log->GetBinLowEdge(1) << "\n";
    std::cout << "last bin = " << data_log->GetBinLowEdge(auxnbins + 1) << "\n";
    int idx = 0;
    for (auto i : mc_xsec) {
      double chi2 = p.m_mnv_plotter.Chi2DataMC(
          data_xsec, i.second, ndf, pot_scale, use_data_error_mtx,
          use_only_shape_errors, use_model_stat);
      double log_chi2 = p.m_mnv_plotter.Chi2DataMC(
          data_log, ln_mc_xsec_arr[idx], nof, pot_scale, use_data_error_mtx,
          use_only_shape_errors, use_model_stat);
      std::cout << i.first << " & " << chi2 << " & " /*<< ndf << " & "
                 << chi2/ndf << " & " << nof << " & "*/
                << log_chi2 <<
          /*" & " << exp(log_chi2)/nof <<*/ "\\\\ \n";
      idx++;
    }
    // std::cout << "   chi2 = "         << chi2     << "\n";
    // std::cout << "   ndf = "          << ndf      << "\n";
    // std::cout << "   chi2/ndf = "     << chi2/ndf << "\n";

    // add label manually
    //    if (p.m_variable->Name() == "tpi") ndf = 6;
    //    if (p.m_variable->Name() == "enu") ndf = 6;
    //    if (p.m_variable->Name() == "pzmu") ndf = 9;
    //    if (p.m_variable->Name() == "pmu") ndf = 8;
    //    if (p.m_variable->Name() == "wexp") ndf = 4;
    //    char* words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf,
    //                       chi2 / (Double_t)ndf);
    //    int align = 33;
    //    p.m_mnv_plotter.AddPlotLabel(words, 0.8, 0.745, 0.03, 1, 62, align);
  }

  // POT info
  // -1 --> don't do mc POT
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, -1, 0.3, 0.88);
  else {
    // p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, -1, 0.3, 0.88);
    p.m_mnv_plotter.WriteNorm("POT-Normalized", 0.3, 0.88, 0.03);
    p.m_mnv_plotter.WriteNorm(Form("Data POT: %.2E", p.m_data_pot), 0.3,
                              0.88 - 0.03, 0.03);
  }

  // p.m_mnv_plotter.WritePreliminary(0.32, 0.812);

  // Change max number of y-axis digits
  // std::cout << "Old max digits = " << TGaxis::GetMaxDigits() << "\n";
  // TGaxis::SetMaxDigits(4);
  // std::cout << "New max digits = " << TGaxis::GetMaxDigits() << "\n";

  // Plot Title
  // minerva plots don't use titles
  // p.m_mnv_plotter.title_size = 0.05;
  // p.SetTitle("Cross Section " + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/CrossSection_%s_%s_%s%s%s_Nuisance", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "pdf");

  // This section is to plot the ratios, where the denominator is the
  // MnvGENIE vX.Y.Z cross section.
  if (plotRatios) {
    TCanvas cR("c2", "c2");
    double x1, y1, x2, y2;
    double legend_text_size = .025;

    p.m_mnv_plotter.DecodeLegendPosition(x1, y1, x2, y2, "TR",
                                         ratio_mc_xsec.size() + 1);

    std::cout << "y1 = " << y1 << "\n";
    TLegend* leg = new TLegend(0.15, 0.70, x2, y2);
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetFillColor(-1);
    leg->SetFillStyle(0);
    leg->SetTextSize(legend_text_size);
    leg->SetTextFont(62);

    double plotMin = 0.5;
    double plotMax = -1;
    if (p.m_variable->Name() == "thetamu_deg") plotMax = 1.8;
    if (p.m_variable->Name() == "q2") plotMax = 2.25;
    if (p.m_variable->Name() == "pzmu") plotMax = 1.6;
    if (p.m_variable->Name() == "ptmu") plotMax = 2.25;
    if (p.m_variable->Name() == "pmu") plotMax = 1.5;
    if (p.m_variable->Name() == "mixthetapi_deg") plotMax = 2;

    p.m_mnv_plotter.axis_title_offset_y = 1.2;
    p.m_mnv_plotter.axis_title_size_y = 0.06;
    p.m_mnv_plotter.DrawDataMCRatio(num_data, (TH1*)mnv_mc_xsec, 1., true,
                                    plotMin, plotMax);
    if (p.m_variable->Name() == "q2") cR.SetLogx();

    p.m_mnv_plotter.ApplyNextLineStyle(num_data, true, true);
    ratio_mc_xsec["MnvTune 4.7.1"]->Divide(ratio_mc_xsec["MnvTune 4.7.1"],
                                           mnv_mc_xsec);
    p.m_mnv_plotter.ApplyNextLineStyle(ratio_mc_xsec["MnvTune 4.7.1"], false,
                                       true);
    ratio_mc_xsec["MnvTune 4.7.1"]->Draw("SAME HIST");
    leg->AddEntry(ratio_mc_xsec["MnvTune 4.7.1"], Form("%s", "MnvTune 4.7.1"),
                  "l");
    //    for (int dum_i = 0; dum_i<3; dum_i++)
    //       p.m_mnv_plotter.ApplyNextLineStyle(ratio_mc_xsec["MnvTune 4.7.1"],
    //       false, true);
    for (auto i : ratio_mc_xsec) {
      if (i.first == "MnvTune 4.7.1") continue;
      i.second->Divide(i.second, mnv_mc_xsec);
      //      i.second->SetLineWidth(5);
      p.m_mnv_plotter.ApplyNextLineStyle(i.second, false, true);
      i.second->Draw("SAME HIST");
      leg->AddEntry(i.second, Form("%s", i.first.c_str()), "l");
    }

    num_data->Divide(num_data, mnv_mc_xsec);
    num_data->SetMarkerStyle(20);
    num_data->SetMarkerSize(1);
    num_data->SetLineWidth(3);
    num_data->SetLineStyle(1);
    num_data->SetLineColor(1);
    num_data->Draw("X0 E1 SAME");
    num_data_stat_error->Divide(num_data_stat_error, mnv_mc_xsec);
    num_data_stat_error->SetMarkerStyle(20);
    num_data_stat_error->SetMarkerSize(1);
    num_data_stat_error->SetLineWidth(3);
    num_data_stat_error->SetLineStyle(1);
    num_data_stat_error->SetLineColor(1);
    num_data_stat_error->Draw("X0 E1 SAME");

    leg->AddEntry(num_data, "Data", "ple");
    leg->Draw();
    std::string outfilename =
        Form("%s/Ratio_CrossSection_%s_%s_%s%s%s_Nuisance", outdir.c_str(),
             p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
             GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
             bwn_str.c_str());
    p.m_mnv_plotter.MultiPrint(&cR, outfilename, "pdf");
  }
  delete data_xsec;
  delete data_xsec_w_tot_error;
  delete data_xsec_w_stat_error;
}

void plot_one_model(Plotter p, MnvH1D* data, TH1D* mc, std::string outdir = ".",
                    double ymax = -1, bool do_log_scale = false,
                    bool do_bin_width_norm = true) {
  std::cout << "Plotting CrossSection " << p.m_variable->Name() << "\n";

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  PlotUtils::MnvH1D* data_xsec = nullptr;
  TH1D* mc_xsec = nullptr;

  if (p.m_variable->Name() == "q2") {
    // if (false) {
    data_xsec = RebinQ2Plot(*data);
    mc_xsec = RebinQ2Plot(*mc);
  } else if (p.m_variable->Name() == "enu" || p.m_variable->Name() == "pmu" ||
             p.m_variable->Name() == "ptmu" || p.m_variable->Name() == "pzmu") {
    data_xsec = RebinningtoGeV(*data, p.m_variable->Name());
    mc_xsec = RebinningtoGeV(*mc, p.m_variable->Name());
  } else {
    data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
    mc_xsec = (PlotUtils::MnvH1D*)mc->Clone("mc");
  }
  /*if (p.m_variable->Name() == "q2") {
    data_xsec = RebinQ2Plot(*data);
    mc_xsec = RebinQ2Plot(*mc);
  } else {
    data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
    mc_xsec = dynamic_cast<TH1D*>(mc->Clone("mc"));
  }*/

  data_xsec->SetTitle("Data");

  TCanvas canvas("c1", "c1");

  // Get Hists
  TH1D* data_xsec_w_tot_error = new TH1D(data_xsec->GetCVHistoWithError());
  TH1D* data_xsec_w_stat_error = new TH1D(data_xsec->GetCVHistoWithStatError());

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
  if (p.m_variable->Name() == "q2") {
    canvas.SetLogx();
  }

  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // Y-label
  p.m_mnv_plotter.axis_title_offset_y = 0.;

  // X label
  p.SetXLabel(data_xsec_w_tot_error);
  p.SetXLabel(data_xsec_w_stat_error);
  p.SetXLabel(mc_xsec);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = 1.;
  }

  // Bin Width Normalization, Y-axis label, and 10^-42 shift
  if (do_bin_width_norm) {
    // data_xsec_w_tot_error ->Scale(1.e38, "width");
    // data_xsec_w_stat_error->Scale(1.e38, "width");
    // mc_xsec_w_stat_error  ->Scale(1.e38, "width");

    data_xsec_w_tot_error->Scale(1.e42, "width");
    data_xsec_w_stat_error->Scale(1.e42, "width");
    mc_xsec->Scale(1.e42, "width");
    if (p.m_variable->Name() == "q2") {
      data_xsec_w_tot_error->SetBinContent(
          2, data_xsec_w_tot_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      data_xsec_w_stat_error->SetBinContent(
          2,
          data_xsec_w_stat_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      mc_xsec->SetBinContent(
          2, mc_xsec->GetBinContent(2) * (0.025 - 0.006) / 0.025);
    }

    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel +
                        " (10^{-42} cm^{2}/" + p.m_variable->m_units +
                        "/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    mc_xsec->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = true;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(data_xsec_w_tot_error, mc_xsec,
                                          pot_scale, "TR", use_hist_titles);

  /*
  // Add chi2 label
  {
    const bool use_data_error_mtx = true;
    const bool use_only_shape_errors = false;
    const bool use_model_stat =
        false;  // model statistical errors -- these are small if you use a
                // priori effden, but atm I'm not. So keep this off.
    // p.m_mnv_plotter.AddChi2Label(data_xsec, mc_xsec, pot_scale, "TR", 0.03,
    // -0.175, use_data_error_mtx, use_only_shape_errors); // this auto turns on
    // model stat errors, I've manually turned it off within the function.

    // std::cout << "use_overflow is " << p.m_mnv_plotter.chi2_use_overflow_err
    // << "\n";
    // Check chi2
    int ndf = -1;
    // double pot_scale = 1.;
    double chi2 = p.m_mnv_plotter.Chi2DataMC(
        data_xsec, mc_xsec, ndf, pot_scale, use_data_error_mtx,
        use_only_shape_errors, use_model_stat);

    // std::cout << p.m_variable->Name() << "\n";
    // std::cout << "   chi2 = "         << chi2     << "\n";
    // std::cout << "   ndf = "          << ndf      << "\n";
    // std::cout << "   chi2/ndf = "     << chi2/ndf << "\n";

    // add label manually
    if (p.m_variable->Name() == "tpi") ndf = 6;
    if (p.m_variable->Name() == "enu") ndf = 6;
    if (p.m_variable->Name() == "pzmu") ndf = 9;
    if (p.m_variable->Name() == "pmu") ndf = 8;
    if (p.m_variable->Name() == "wexp") ndf = 4;
    char* words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf,
                       chi2 / (Double_t)ndf);
    int align = 33;
    p.m_mnv_plotter.AddPlotLabel(words, 0.8, 0.745, 0.03, 1, 62, align);
  }
  */

  // POT info
  // -1 --> don't do mc POT
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, -1, 0.3, 0.88);
  else {
    // p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, -1, 0.3, 0.88);
    p.m_mnv_plotter.WriteNorm("POT-Normalized", 0.3, 0.88, 0.03);
    p.m_mnv_plotter.WriteNorm(Form("Data POT: %.2E", p.m_data_pot), 0.3,
                              0.88 - 0.03, 0.03);
  }

  //// Label re: error bars
  // p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.6, 0.772, 0.03);

  // Change max number of y-axis digits
  // std::cout << "Old max digits = " << TGaxis::GetMaxDigits() << "\n";
  // TGaxis::SetMaxDigits(4);
  // std::cout << "New max digits = " << TGaxis::GetMaxDigits() << "\n";

  // Plot Title
  // minerva plots don't use titles
  // p.m_mnv_plotter.title_size = 0.05;
  // p.SetTitle("Cross Section " + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/CrossSection_%s_%s_%s%s%s_Nuisance", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");

  delete data_xsec;
  delete mc_xsec;
  delete data_xsec_w_tot_error;
  delete data_xsec_w_stat_error;
}

void set_POT(TFile& fin, CCPi::MacroUtil& util) {
  util.m_mc_pot = -1;
  util.m_data_pot = -1;

  // Data
  // TODO make this safer
  PlotUtils::MnvH1D* h_data_pot = (PlotUtils::MnvH1D*)fin.Get("data_pot");
  float data_pot = h_data_pot->GetBinContent(1);
  util.m_data_pot = data_pot;

  // MC
  PlotUtils::MnvH1D* h_mc_pot = (PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float mc_pot = h_mc_pot->GetBinContent(1);
  util.m_mc_pot = mc_pot;

  // Ratio
  util.m_pot_scale = util.m_data_pot / util.m_mc_pot;
}

//==============================================================================
// Main
//==============================================================================
void plotNuisance(int signal_definition_int = 1, int plot_errors = 0) {
  // Data xsec file input
  TFile fin1(
      "DataXSecInputs_20241204_ALL_mix_CCPiWeightSys_p4.root",
      //   "DataXSecInputs_20241217_ALL_mix_thetapi_CCPiWeighSys_p4.root",
      "READ");
  // TFile
  // fin1("DataXSecInputs_20240610_ALL_thetapisig_NewEstimatorptmucut_Sys_p4.root",
  // "READ");
  std::cout << "Reading data input from " << fin1.GetName() << "\n";

  std::map<std::string, std::string> models;
  models["MnvTune 4.7.1"] =
      "DataXSecInputs_20241217_ALL_mix_thetapi_CCPiWeighSys_p4.root";

  models["GENIE v3 hA"] =
      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
      "nuisance_ME_FHC_tracker_GENIE_v3_0_6_G18_02a_02_11a_CH_50M.root";
  models["GENIE v3 hN"] =
      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
      "nuisance_ME_FHC_tracker_GENIE_v3_0_6_G18_02b_02_11a_CH_50M.root";
  models["GENIE v2.12.6"] =
      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
      "nuisance_ME_FHC_tracker_GENIE_v2_12_6.root";

  models["NEUT LFG"] =
      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
      "nuisance_ME_FHC_tracker_NEUT_v5_4_1_LFG_ma105.root";
  models["GiBUU"] =
      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
      "CCPion_GiBUU_T0_tracker.root";
  models["NuWro LFG"] =
      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
      "nuisance_ME_FHC_tracker_NuWro_1902_LFG.root";

  //  models["NEUT SF"] = // NEUT SF
  //      "/exp/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/nuisance/"
  //      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
  //      "nuisance_ME_FHC_tracker_NEUT_v5_4_1_SF_ma103.root";
  //  models["NuWro SF"] = // NuWro SF
  //      "/exp/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/nuisance/"
  //      "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/Nuisance/"
  //      "nuisance_ME_FHC_tracker_NuWro_1902_SF.root";

  const std::string plist = "ME1L";
  std::string data_file_list = GetPlaylistFile(plist, false, false);
  std::string mc_file_list = GetPlaylistFile(plist, true, false);

  // Macro Utility
  const std::string macro("Plot Nuisance");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);

  // Get POT from file, not from any chain
  set_POT(fin1, util);
  util.PrintMacroConfiguration(macro);

  // Initialize variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);

  // Load histogram contents into variables from input file.
  for (auto var : variables) {
    if (var->Name() == sidebands::kFitVarString) {
      var->m_hists.m_stacked_wsideband = StackedHistogram<WSidebandType>(
          var->m_hists.m_label, var->m_hists.m_xlabel,
          var->m_hists.m_bins_array, kNWSidebandTypes,
          sidebands::kWSideband_ColorScheme);
    }
    var->LoadMCHistsFromFile(fin1, util.m_error_bands);
    var->LoadDataHistsFromFile(fin1);
  }

  // remove unwanted variables
  ContainerEraser::erase_if(variables, [](Variable* v) {
    return v->Name() == "tpi_mbr" || v->Name() == "wexp_fit" ||
           v->Name() == "ehad";
  });

  const bool do_frac_unc = true;
  const bool include_stat = true;
  const bool do_cov_area_norm = false;
  const double ymax = -1.;
  const bool do_log_scale = false;
  bool plotRatios = true;
  for (auto reco_var : variables) {
    if (reco_var->m_is_true) continue;

    bool do_bin_width_norm = true;
    // Set up plotting properties
    Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                      do_cov_area_norm, include_stat, util.m_signal_definition);

    // Get nuisance prediction from input file
    std::string var_name;
    if (reco_var->Name() == "mixtpi")
      var_name = "tpi";
    else if (reco_var->Name() == "mixthetapi_deg")
      var_name = "thpi";
    else if (reco_var->Name() == "thetamu_deg")
      var_name = "thmu";
    else
      var_name = reco_var->Name();
    if (reco_var->Name() == "enu") do_bin_width_norm = false;
    PlotUtils::MnvH1D* Mnv_mc_cross_section = (PlotUtils::MnvH1D*)fin1.Get(
        Form("mc_cross_section_%s", reco_var->Name().c_str()));
    // Plot data vs a single model
    if (false) {
      auto model = models.begin();
      std::cout << model->first << " fin " << model->second << "\n";
      TFile fin(model->second.c_str(), "READ");
      TH1D* mc_cross_section = (TH1D*)fin.Get(var_name.c_str());
      assert(mc_cross_section);
      std::cout << mc_cross_section->Integral() << "  "
                << mc_cross_section->GetEntries() << "\n";
      mc_cross_section->SetTitle(model->first.c_str());
      plot_one_model(plot_info, reco_var->m_hists.m_cross_section,
                     mc_cross_section);
      delete mc_cross_section;
      fin.Close();
    }

    // Plot data vs all models on a single plot
    if (true) {
      // First get the model xsec plots from their files
      std::map<std::string, PlotUtils::MnvH1D*> mc_cross_sections;
      for (const auto& model : models) {
        std::string description = model.first;
        TFile fin(model.second.c_str(), "READ");
        std::cout << "Varname = " << var_name << "\n";
        TH1D* mc_cross_section = nullptr;
        if (description == "MnvTune 4.7.1")
          mc_cross_section = (TH1D*)Mnv_mc_cross_section->Clone("mnvTune");
        else
          mc_cross_section = (TH1D*)fin.Get(var_name.c_str());
        assert(mc_cross_section);
        if (description == "GiBUU") mc_cross_section->Scale(1e-38);
        mc_cross_section->SetTitle(description.c_str());
        PlotUtils::MnvH1D* mc_cross_section_mnvh1d =
            new PlotUtils::MnvH1D(*mc_cross_section);
        mc_cross_sections[model.first] = mc_cross_section_mnvh1d;
        fin.Close();
      }
      // Plot
      plot_all_models(plot_info, reco_var->m_hists.m_cross_section,
                      mc_cross_sections, Mnv_mc_cross_section, ".", -1, false,
                      do_bin_width_norm, plotRatios);
    }
  }
  fin1.Close();

  //============================================================================
}

#endif  // plotNuisance_C
