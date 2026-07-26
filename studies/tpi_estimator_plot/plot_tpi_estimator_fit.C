#include <iostream>

#include "../includes/Plotter.h"
#include "MnvPlotter.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"

//  KEY: TGraph	TpiEst;1	Graph
//  KEY: TGraphErrors	TpiEst_Errs;1	Graph
//  KEY: TGraphErrors	TpiEst_Errs_Fit;1	Graph
//  KEY: TGraph	FitCurve;1	Graph

/*
int PlotTogether(TH1* h1, std::string label1, TH1* h2, std::string label2,
                 std::string tag, double ymax = -1, bool do_log_scale = false,
                 bool do_fit = false, std::string ylabel = "") {
  std::cout << "PlotTogether" << std::endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  // TH1::SetDefaultSumw2();

  TCanvas cF("c4", "c4");
  cF.Update();

  TLegend* leg = new TLegend(0.75, 0.85, 0.95, 0.95, NULL, "brNDC");

  if (h1->GetMaximum() > h2->GetMaximum()) {
    std::cout << "h1 > h2  " << h1->GetMaximum() << "  " << h2->GetMaximum()
              << "\n";

    h1->SetLineWidth(3);
    h1->SetLineColor(kBlack);
    h1->Draw("HIST");

    cF.Update();

    if (ymax > 0) h1->GetYaxis()->SetRangeUser(0, ymax);
    h2->GetYaxis()->SetTitle(ylabel.c_str());
    h2->SetLineColor(kRed);
    h2->SetLineWidth(3);
    h2->Draw("HISTSAME");
  } else {
    std::cout << "h1 < h2  " << h1->GetMaximum() << "  " << h2->GetMaximum()
              << "\n";
    h2->SetLineWidth(3);
    h2->SetLineColor(kRed);
    h2->Draw("HIST");

    cF.Update();

    if (ymax > 0) h2->GetYaxis()->SetRangeUser(0, ymax);
    h1->GetYaxis()->SetTitle(ylabel.c_str());
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(3);
    h1->Draw("HISTSAME");
  }

  if (do_log_scale) cF.SetLogy();

  cF.Update();

  //// draw a line
  // TLine *line = new TLine(0,0,0,cF.GetUymax());
  // line->SetLineColor(kBlue);
  // line->Draw();

  // legend
  TLegendEntry* entry = leg->AddEntry("NULL", "", "h");
  leg->SetTextSize(0.035);
  {
    entry = leg->AddEntry("entry", label1.c_str(), "l");
    entry->SetLineWidth(5);
    entry->SetLineColor(kBlack);
  }
  {
    entry = leg->AddEntry("entry", label2.c_str(), "l");
    entry->SetLineWidth(5);
    entry->SetLineColor(kRed);
  }
  leg->Draw();

  cF.Print(Form("%s.png", tag.c_str()));
  // cF.Print(Form("%s.eps", tag.c_str()));

  delete leg;

  return 0;
}
*/

void PlotDataMCWithError(PlotUtils::MnvH1D* mc, PlotUtils::MnvH1D* data,
                         Plotter p, std::string tag) {
  std::cout << "Plotting\n";
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  const bool use_hist_titles = true;
  // PlotUtils::MnvH1D* tmp_bg = nullptr;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(data, mc, pot_scale, "TR");
  //, use_hist_titles,
  // p.m_variable->m_hists.m_bg.hist,
  // tmp_bg, NULL, p.m_do_cov_area_norm, p.m_include_stat);
  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "png");
}

int plot_tpi_estimator_fit() {
  // Infiles
  TFile fin("studies/tpi_estimator_plot/tpi_from_range_fit.root", "READ");
  std::cout << "Reading input from " << fin.GetName() << std::endl;

  TGraphErrors* data = (TGraphErrors*)fin.Get("TpiEst_Errs");
  TGraphErrors* fit = (TGraphErrors*)fin.Get("TpiEst_Errs_Fit");

  //===========================================================================
  // attempt to just plot the tgraph
  //===========================================================================
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas cF("c4", "c4");
  cF.Update();

  data->Draw("AP");
  fit->Draw("PSAME");

  data->SetTitle("");
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.4);
  data->GetXaxis()->SetLabelFont(63);  // Bold for x-axis labels
  data->GetYaxis()->SetLabelFont(63);  // Bold for y-axis labels
  data->GetXaxis()->SetTitleFont(63);
  data->GetYaxis()->SetTitleFont(63);

  cF.Update();
  cF.Draw();

  cF.Print("tpi_estimator_fit.png");

  /*
    //===========================================================================
    // attempt to convert tgraph into th1's
    //===========================================================================
    // convert data to MnvH1D
    int n_bins = data->GetN();
    PlotUtils::MnvH1D* data_hist = new PlotUtils::MnvH1D(
        "hist", "tpi_vs_range_hist", n_bins,
        data->GetX()[0] - 0.5 * (data->GetErrorX(0)),
        data->GetX()[n_bins - 1] + 0.5 * (data->GetErrorX(n_bins - 1)));

    for (int i = 0; i < n_bins; ++i) {
      double x = data->GetX()[i];
      double y = data->GetY()[i];
      double ex = data->GetErrorX(i);
      double ey = data->GetErrorY(i);

      int bin = data_hist->FindBin(x);   // Find corresponding bin for x
      data_hist->SetBinContent(bin, y);  // Set bin content
      data_hist->SetBinError(bin, ey);   // Set bin error
    }

    // convert fit to MnvH1D
    PlotUtils::MnvH1D* fit_hist = new PlotUtils::MnvH1D(
        "fit_hist", "tpi_vs_range_fit_hist", n_bins,
        fit->GetX()[0] - 0.5 * (fit->GetErrorX(0)),
        fit->GetX()[n_bins - 1] + 0.5 * (fit->GetErrorX(n_bins - 1)));

    for (int i = 0; i < n_bins; ++i) {
      double x = fit->GetX()[i];
      double y = fit->GetY()[i];
      double ex = fit->GetErrorX(i);
      double ey = fit->GetErrorY(i);

      int bin = fit_hist->FindBin(x);   // Find corresponding bin for x
      fit_hist->SetBinContent(bin, y);  // Set bin content
      fit_hist->SetBinError(bin, ey);   // Set bin error
    }

    PlotUtils::MnvPlotter mnv_plotter(kCCNuPionIncStyle);

    TCanvas canvas("c1", "c1");

    mnv_plotter.DrawDataMCWithErrorBand(data_hist, fit_hist, 1, "TR");

    mnv_plotter.MultiPrint(&canvas, "tpi_estimator_fit", "png");
  */

  // PlotUtils::MnvH1D* AaronXsecGENIE =
  // (PlotUtils::MnvH1D*)fin1.Get("q2_xsec"); Plotter plot_info(reco_var,
  // util.m_mc_pot, util.m_data_pot, do_frac_unc,
  //                   do_cov_area_norm, include_stat,
  //                   util.m_signal_definition);
  // Plot_CrossSection(plot_info, reco_var->m_hists.m_cross_section,
  //                   m_mc_cross_section, ".", -1, false, do_bin_width_norm);

  return 0;
}
