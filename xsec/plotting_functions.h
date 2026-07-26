#ifndef plotting_functions_h
#define plotting_functions_h

#include <cassert>  // !! must compile in debug mode root script.C++g
#include <iostream>
#include <sstream>

#include "../includes/myPlotStyle.h"
#include "PlotUtils/HistogramUtils.h"
#include "PlotUtils/MnvColors.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/TargetUtils.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif
#include "Constants.h"  // enum SignalDefinition
#include "PlotUtils/MnvVertErrorBand.h"
#include "Plotter.h"
#include "SignalDefinition.h"  // GetSignalFileTag
#include "Systematics.h"       // namespace systematics
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TPad.h"
#include "TPaveStats.h"
//#include "TStyle.h"
#include "TArrayD.h"
#include "TText.h"
#include "Variable.h"

class Variable;

PlotUtils::MnvH1D* LogHist(const PlotUtils::MnvH1D& old_hist, std::string var) {
  PlotUtils::MnvH1D* log_hist = nullptr;
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)old_hist.Clone("hist");
  log_hist = (PlotUtils::MnvH1D*)old_hist.Clone("log_hist");
  log_hist->Reset();
  log_hist->ClearAllErrorBands();
  int nbins = old_hist.GetNbinsX();
  for (int i = 1; i <= nbins; i++) {
    double val = old_hist.GetBinContent(i);
    double err = old_hist.GetBinError(i);
    if (val == 0.)
      continue;

    else {
      log_hist->SetBinContent(i, log(val));
      if (err != 0.)
        log_hist->SetBinError(i, err / val);
      else
        log_hist->SetBinError(i, 0.);
    }
  }

  // ASSERT CV
  for (int i = 1; i < nbins + 1; i++) {
    double val = old_hist.GetBinContent(i);
    double err = old_hist.GetBinError(i);
    if (val == 0)
      continue;

    else {
      val = log(val);
      if (err != 0)
        err = err / old_hist.GetBinContent(i);
      else
        err = 0.;
    }
    assert(log_hist->GetBinContent(i) == val);
    assert(log_hist->GetBinError(i) == abs(err));
  }
  // Universes
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();
    log_hist->AddVertErrorBand(error_name, n_univs);

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      TH1* univ_i_hist_new =
          log_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));

      for (int i = 0; i <= nbins; i++) {
        double unival = univ_i_hist_old->GetBinContent(i);
        if (unival == 0)
          continue;

        else {
          univ_i_hist_new->SetBinContent(
              i, log(univ_i_hist_old->GetBinContent(i)));
          if (univ_i_hist_old->GetBinError(i) != 0.)
            univ_i_hist_new->SetBinError(
                i, univ_i_hist_old->GetBinError(i) / unival);
          else
            univ_i_hist_new->SetBinError(i, 0.);
        }
      }
      delete univ_i_hist_old;
    }
  }
  // ASSERT UNIVERSES
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();
    // std::cout << error_name << "\n";

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      // std::cout << "  " << univ_i << "\n";

      TH1* univ_i_hist_new =
          log_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));

      for (int i = 0; i <= nbins; i++) {
        double val = univ_i_hist_old->GetBinContent(i);
        double err = univ_i_hist_old->GetBinError(i);
        if (val == 0.)
          continue;

        else {
          val = log(val);
          if (err != 0.)
            err = err / univ_i_hist_old->GetBinContent(i);
          else
            err = 0.;
        }
        /*std::cout << "    " << univ_i_hist_new->GetBinContent(i) <<
        " = " << val <<  " | "
                 << univ_i_hist_new->GetBinError(i) << " = " <<
                 err << "\n";*/
        assert(univ_i_hist_new->GetBinContent(i) == val);
        assert(univ_i_hist_new->GetBinError(i) == abs(err));
      }
      delete univ_i_hist_old;
    }
  }

  for (auto error_name : old_hist.GetLatErrorBandNames()) {
    int n_univs = old_hist.GetLatErrorBand(error_name)->GetNHists();
    log_hist->AddLatErrorBand(error_name, n_univs);
    // std::cout << "Lateral error band = " << error_name << "\n";
    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      TH1* lat_univ_i_hist_new =
          log_hist->GetLatErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetLatErrorBand(error_name)->GetHist(univ_i));

      for (int i = 0; i <= nbins; i++) {
        double unival = univ_i_hist_old->GetBinContent(i);
        if (unival == 0)
          continue;

        else {
          lat_univ_i_hist_new->SetBinContent(
              i, log(univ_i_hist_old->GetBinContent(i)));
          if (univ_i_hist_old->GetBinError(i) != 0.)
            lat_univ_i_hist_new->SetBinError(
                i, univ_i_hist_old->GetBinError(i) / unival);
          else
            lat_univ_i_hist_new->SetBinError(i, 0.);
        }
      }
      delete univ_i_hist_old;
    }
  }

  //  cout << "Found it " <<
  //  old_hist.HasErrorMatrix(Form("unfolding_cov_matrix_%s",
  //			 var.c_str())) << endl;
  // unfoldingCov
  // loop over error matrices ...
  std::vector<std::string> errorMatrixNames =
      old_hist.GetSysErrorMatricesNames();
  //  for( std::vector<std::string>::const_iterator i =
  //  errorMatrixNames.begin(); i != errorMatrixNames.end(); ++i ){
  for (auto error_name : old_hist.GetSysErrorMatricesNames()) {
    if (old_hist.HasErrorMatrix(error_name)) {
      TMatrixD covmx = old_hist.GetSysErrorMatrix(error_name);
      int nx = covmx.GetNrows();
      for (int ix = 0; ix < nx; ix++) {
        double vali = old_hist.GetBinContent(ix);
        for (int jx = 0; jx < nx; jx++) {
          double valj = old_hist.GetBinContent(jx);
          if (vali != 0. && valj != 0.)
            covmx[ix][jx] *= 1 / (vali * valj);
          else
            covmx[ix][jx] *= 0.;
        }
      }
      log_hist->FillSysErrorMatrix(error_name, covmx);
    }
  }

  return log_hist;
}
// Make the q2 plot with aaron appearance.
// add a 0th bin at the beginning from (0, epsilon) so it doesn't look weird in
// log scale
PlotUtils::MnvH1D* RebinQ2Plot(const PlotUtils::MnvH1D& old_hist) {
  // make a new bin array and convert mev^2 TO gev^2
  // old bins
  const TArrayD old_bins_array = *(old_hist.GetXaxis()->GetXbins());
  const int n_old_bins = old_bins_array.GetSize();

  //  std::cout << "size of old hist bin array " << n_old_bins << "\n";

  TArrayD new_bins_array = old_bins_array;

  // increase number of bins by 1.
  // This "Set" adds a new 0 at the beginning of the array.
  new_bins_array.Set(n_old_bins + 1);
  const int n_new_bins = new_bins_array.GetSize();
  // std::cout << "size of new hist bin array " << n_new_bins << "\n";

  // Scale to GeV^2
  for (int i = n_new_bins - 1; i >= 2; i--) {
    new_bins_array[i] = new_bins_array[i - 1] * 1.e-6;
  }

  // Setting the first bin edge to 0.006 makes the q2 plot look good in GeV^2
  // and log scale. Also, it's what Aaron uses.
  new_bins_array[1] = 6.e-3;  // <-- THIS DETERMINES HOW BIG FIRST BIN APPEARS
  new_bins_array[0] = 0.;

  // manually make the new MH1D
  std::string new_name = Form("%s_%s", old_hist.GetName(), "_rebin");
  PlotUtils::MnvH1D* new_hist = new PlotUtils::MnvH1D(
      new_name.c_str(), old_hist.GetTitle(), new_bins_array.GetSize() - 1,
      new_bins_array.GetArray());

  new_hist->SetLineColor(old_hist.GetLineColor());
  new_hist->SetLineStyle(old_hist.GetLineStyle());
  new_hist->SetLineWidth(old_hist.GetLineWidth());

  new_hist->SetMarkerColor(old_hist.GetMarkerColor());
  new_hist->SetMarkerStyle(old_hist.GetMarkerStyle());
  new_hist->SetMarkerSize(old_hist.GetMarkerSize());

  new_hist->SetTitle(old_hist.GetTitle());
  new_hist->GetXaxis()->SetTitle(old_hist.GetXaxis()->GetTitle());
  new_hist->GetYaxis()->SetTitle(old_hist.GetYaxis()->GetTitle());

  // finally, move contents, bin-by-bin, universe-by-universe from old to new
  // WARNING: THIS IS A PLOTTING HACK. THIS HIST'S 0TH AND 1ST BINS ARE NO
  // TECHNICALLY CORRECY
  // CV
  for (int i = 0; i < n_old_bins; i++) {
    int new_bin_idx = i + 1;
    new_hist->SetBinContent(new_bin_idx, old_hist.GetBinContent(i));
    new_hist->SetBinError(new_bin_idx, old_hist.GetBinError(i));
  }
  new_hist->SetBinContent(0, 0.);
  new_hist->SetBinError(0, 0.);

  // ASSERT CV
  for (int i = 0; i < n_old_bins; i++) {
    int new_bin_idx = i + 1;
    assert(new_hist->GetBinContent(new_bin_idx) == old_hist.GetBinContent(i));
    assert(new_hist->GetBinError(new_bin_idx) == old_hist.GetBinError(i));
  }
  // Universes
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();
    new_hist->AddVertErrorBand(error_name, n_univs);

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      TH1* univ_i_hist_new =
          new_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));

      for (int i = 0; i < n_old_bins; i++) {
        int new_bin_idx = i + 1;
        univ_i_hist_new->SetBinContent(new_bin_idx,
                                       univ_i_hist_old->GetBinContent(i));
        univ_i_hist_new->SetBinError(new_bin_idx,
                                     univ_i_hist_old->GetBinError(i));
      }

      univ_i_hist_new->SetBinContent(0, 0.);
      univ_i_hist_new->SetBinError(0, 0.);
      delete univ_i_hist_old;
    }
  }

  // ASSERT UNIVERSES
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();
    // std::cout << error_name << "\n";

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      // std::cout << "  " << univ_i << "\n";

      TH1* univ_i_hist_new =
          new_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));

      for (int i = 0; i < n_old_bins; i++) {
        int new_bin_idx = i + 1;
        // std::cout << "    " << univ_i_hist_new->GetBinContent(new_bin_idx) <<
        // " = " << univ_i_hist_old->GetBinContent(i) <<  " | "
        //          << univ_i_hist_new->GetBinError(new_bin_idx) << " = " <<
        //          univ_i_hist_old->GetBinError(i) << "\n";
        assert(univ_i_hist_new->GetBinContent(new_bin_idx) ==
               univ_i_hist_old->GetBinContent(i));
        assert(univ_i_hist_new->GetBinError(new_bin_idx) ==
               univ_i_hist_old->GetBinError(i));
      }
      delete univ_i_hist_old;
    }
  }
  for (auto error_name : old_hist.GetUncorrErrorNames()) {
    int n_univs = 1;
    new_hist->AddUncorrError(error_name);
    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      TH1* univ_i_hist_new = new_hist->GetUncorrError(error_name);
      TH1D* univ_i_hist_old = new TH1D(*old_hist.GetUncorrError(error_name));
      for (int i = 0; i < n_old_bins; i++) {
        int new_bin_idx = i + 1;
        univ_i_hist_new->SetBinContent(new_bin_idx,
                                       univ_i_hist_old->GetBinContent(i));
        univ_i_hist_new->SetBinError(new_bin_idx,
                                     univ_i_hist_old->GetBinError(i));
      }

      univ_i_hist_new->SetBinContent(0, 0.);
      univ_i_hist_new->SetBinError(0, 0.);
      delete univ_i_hist_old;
    }
  }
  for (auto error_name : old_hist.GetUncorrErrorNames()) {
    // for (auto error_name : error_names) {
    int n_univs = 1;

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      //     std::cout << "  " << univ_i << "\n";

      TH1* univ_i_hist_new = new_hist->GetUncorrError(error_name);
      TH1D* univ_i_hist_old = new TH1D(*old_hist.GetUncorrError(error_name));
      for (int i = 0; i < n_old_bins; i++) {
        int new_bin_idx = i + 1;
        // std::cout << "    " << univ_i_hist_new->GetBinContent(new_bin_idx) <<
        // " = " << univ_i_hist_old->GetBinContent(i) <<  " | "
        //          << univ_i_hist_new->GetBinError(new_bin_idx) << " = " <<
        //          univ_i_hist_old->GetBinError(i) << "\n";
        assert(univ_i_hist_new->GetBinContent(new_bin_idx) ==
               univ_i_hist_old->GetBinContent(i));
        assert(univ_i_hist_new->GetBinError(new_bin_idx) ==
               univ_i_hist_old->GetBinError(i));
      }
      delete univ_i_hist_old;
    }
  }

  return new_hist;
}

PlotUtils::MnvH1D* RebinningtoGeV(const PlotUtils::MnvH1D& old_hist,
                                  std::string var) {
  int nbins = old_hist.GetNbinsX();
  int x = 1;
  int start = 0;
  if (var == "pmu") {
    nbins = nbins - 1;
    x = 2;
    start = 1;
  }
  TArrayD old_bins_array = *(old_hist.GetXaxis()->GetXbins());
  TArrayD new_bins_array = old_bins_array;
  std::string name = old_hist.GetName();
  for (int i = 0; i < nbins + 1; i++) {
    new_bins_array[i] = new_bins_array[i + start] / 1000;
  }
  std::string new_name = Form("%s_%s", old_hist.GetName(), "_rebin");
  PlotUtils::MnvH1D* new_hist = new PlotUtils::MnvH1D(
      new_name.c_str(), old_hist.GetTitle(), new_bins_array.GetSize() - x,
      new_bins_array.GetArray());

  new_hist->SetLineColor(old_hist.GetLineColor());
  new_hist->SetLineStyle(old_hist.GetLineStyle());
  new_hist->SetLineWidth(old_hist.GetLineWidth());

  new_hist->SetFillColor(old_hist.GetFillColor());
  new_hist->SetFillStyle(old_hist.GetFillStyle());

  new_hist->SetMarkerColor(old_hist.GetMarkerColor());
  new_hist->SetMarkerStyle(old_hist.GetMarkerStyle());
  new_hist->SetMarkerSize(old_hist.GetMarkerSize());

  new_hist->SetTitle(old_hist.GetTitle());
  new_hist->GetXaxis()->SetTitle("GeV^{2}");
  new_hist->GetYaxis()->SetTitle(old_hist.GetYaxis()->GetTitle());

  // finally, move contents, bin-by-bin, universe-by-universe from old to new
  // WARNING: THIS IS A PLOTTING HACK. THIS HIST'S 0TH AND 1ST BINS ARE NO
  // TECHNICALLY CORRECY
  // CV
  for (int i = 1; i <= nbins + 1; i++) {
    new_hist->SetBinContent(i, old_hist.GetBinContent(i + start));
    new_hist->SetBinError(i, old_hist.GetBinError(i + start));
  }
  new_hist->SetBinContent(0, 0.);
  new_hist->SetBinError(0, 0.);

  // ASSERT CV
  for (int i = 1; i < nbins + 1; i++) {
    assert(new_hist->GetBinContent(i) == old_hist.GetBinContent(i + start));
    assert(new_hist->GetBinError(i) == old_hist.GetBinError(i + start));
  }
  // Universes
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();
    new_hist->AddVertErrorBand(error_name, n_univs);
    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      TH1* univ_i_hist_new =
          new_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));
      for (int i = 1; i < nbins + 1; i++) {
        univ_i_hist_new->SetBinContent(
            i, univ_i_hist_old->GetBinContent(i + start));
        univ_i_hist_new->SetBinError(i,
                                     univ_i_hist_old->GetBinError(i + start));
      }

      univ_i_hist_new->SetBinContent(0, 0.);
      univ_i_hist_new->SetBinError(0, 0.);
      delete univ_i_hist_old;
    }
  }
  // ASSERT UNIVERSES
  for (auto error_name : old_hist.GetVertErrorBandNames()) {
    // for (auto error_name : error_names) {
    int n_univs = old_hist.GetVertErrorBand(error_name)->GetNHists();

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      //      std::cout << "  " << univ_i << "\n";
      TH1* univ_i_hist_new =
          new_hist->GetVertErrorBand(error_name)->GetHist(univ_i);
      TH1D* univ_i_hist_old =
          new TH1D(*old_hist.GetVertErrorBand(error_name)->GetHist(univ_i));

      for (int i = 1; i < nbins + 1; i++) {
        int new_bin_idx = i;
        //         std::cout << "    " <<
        //         univ_i_hist_new->GetBinContent(new_bin_idx) << " = " <<
        //         univ_i_hist_old->GetBinContent(i) <<  " | "
        //                  << univ_i_hist_new->GetBinError(new_bin_idx) << " =
        //                  " << univ_i_hist_old->GetBinError(i) << "\n";
        assert(univ_i_hist_new->GetBinContent(new_bin_idx) ==
               univ_i_hist_old->GetBinContent(i + start));
        assert(univ_i_hist_new->GetBinError(new_bin_idx) ==
               univ_i_hist_old->GetBinError(i + start));
      }
      delete univ_i_hist_old;
    }
  }

  for (auto error_name : old_hist.GetUncorrErrorNames()) {
    int n_univs = 1;
    new_hist->AddUncorrError(error_name);
    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      TH1* univ_i_hist_new = new_hist->GetUncorrError(error_name);
      TH1D* univ_i_hist_old = new TH1D(*old_hist.GetUncorrError(error_name));
      for (int i = 1; i < nbins + 1; i++) {
        univ_i_hist_new->SetBinContent(
            i, univ_i_hist_old->GetBinContent(i + start));
        univ_i_hist_new->SetBinError(i,
                                     univ_i_hist_old->GetBinError(i + start));
      }

      univ_i_hist_new->SetBinContent(0, 0.);
      univ_i_hist_new->SetBinError(0, 0.);
      delete univ_i_hist_old;
    }
  }
  for (auto error_name : old_hist.GetUncorrErrorNames()) {
    // for (auto error_name : error_names) {
    int n_univs = 1;

    for (int univ_i = 0; univ_i < n_univs; ++univ_i) {
      //       std::cout << "  " << univ_i << "\n";

      TH1* univ_i_hist_new = new_hist->GetUncorrError(error_name);
      TH1D* univ_i_hist_old = new TH1D(*old_hist.GetUncorrError(error_name));
      for (int i = 1; i < nbins + 1; i++) {
        int new_bin_idx = i;
        // std::cout << "    " << univ_i_hist_new->GetBinContent(new_bin_idx) <<
        // " = " << univ_i_hist_old->GetBinContent(i) <<  " | "
        //          << univ_i_hist_new->GetBinError(new_bin_idx) << " = " <<
        //          univ_i_hist_old->GetBinError(i) << "\n";
        assert(univ_i_hist_new->GetBinContent(new_bin_idx) ==
               univ_i_hist_old->GetBinContent(i + start));
        assert(univ_i_hist_new->GetBinError(new_bin_idx) ==
               univ_i_hist_old->GetBinError(i + start));
      }
      delete univ_i_hist_old;
    }
  }
  return new_hist;
}

PlotUtils::MnvH1D* UndoBWN(PlotUtils::MnvH1D* h) {
  PlotUtils::MnvH1D* binw = h->Clone("BinWidth");
  for (int i = 1; i <= (int)h->GetNbinsX(); i++)
    binw->SetBinContent(i, h->GetBinWidth(i));

  h->MultiplySingle(h, binw);
  return h;
}

//==============================================================================
// Some Systematics General Functions
//==============================================================================

void SetErrorGroups(MnvPlotter& mnv_plotter, bool is_subgroups) {
  mnv_plotter.error_summary_group_map.clear();
  if (!is_subgroups) {
    mnv_plotter.error_summary_group_map["LowQ2Pi"].push_back("LowQ2Pi");
    mnv_plotter.error_summary_group_map["Muon"].push_back("Muon_Energy_MINOS");
    mnv_plotter.error_summary_group_map["Muon"].push_back(
        "Muon_Energy_MINERvA");
    mnv_plotter.error_summary_group_map["Muon"].push_back(
        "Muon_Energy_Resolution");
    mnv_plotter.error_summary_group_map["Muon"].push_back(
        "MINOS_Reconstruction_Efficiency");
    mnv_plotter.error_summary_group_map["Muon"].push_back(
        "MuonAngleXResolution");
    mnv_plotter.error_summary_group_map["Muon"].push_back(
        "MuonAngleYResolution");
    mnv_plotter.error_summary_group_map["Muon"].push_back("MuonResolution");
    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "MichelEfficiency");
    //  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_D2_MaRES");
    //  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_EP_MvRES");
    //  mnv_plotter.error_summary_group_map["GENIE"].push_back("GENIE_MaCCQE");
    mnv_plotter.error_summary_group_map["Flux"].push_back("Flux");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_CH");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_C");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_Fe");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_H2O");
    mnv_plotter.error_summary_group_map["Others"].push_back("Target_Mass_Pb");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_em");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_meson");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_other");
    mnv_plotter.error_summary_group_map["Others"].push_back("response_proton");
    mnv_plotter.error_summary_group_map["Others"].push_back("GEANT_Proton");
    mnv_plotter.error_summary_group_map["Others"].push_back("GEANT_Pion");
    mnv_plotter.error_summary_group_map["Others"].push_back("GEANT_Neutron");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "GENIE_D2_NormCCRES");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "DiffractiveModelUnc");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "CoherentPiUnc_C");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "CoherentPiUnc_CH");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "CoherentPiUnc_Fe");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "CoherentPiUnc_H2O");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "CoherentPiUnc_Pb");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "GENIE_Rvn1pi");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "GENIE_Rvp1pi");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "GENIE_Rvn2pi");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "GENIE_Rvp2pi");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "RPA_LowQ2");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "RPA_HighQ2");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "Low_Recoil_2p2h_Tune");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "CCPi+ Tune");
    for (auto g : systematics::kGenieSystematics_InteractionModel)
      mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(g);
    // for(auto g : systematics::kGenieSystematics_FSI)
    //  mnv_plotter.error_summary_group_map["Genie_FSI"].push_back(g);

    for (auto g : systematics::kGenieSystematics_FSI_nucleons)
      mnv_plotter.error_summary_group_map["GENIE FSI"].push_back(g);

    for (auto g : systematics::kGenieSystematics_FSI_pions)
      mnv_plotter.error_summary_group_map["GENIE FSI"].push_back(g);

    // mnv_plotter.error_summary_group_map["Detector"].push_back("EmuRangeCurve");
    // mnv_plotter.error_summary_group_map["Detector"].push_back("PartResp");
    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "TrackAngle");
    mnv_plotter.error_summary_group_map["Muon"].push_back("BeamAngle");
    mnv_plotter.error_summary_group_map["Muon"].push_back("BeamAngleX");
    mnv_plotter.error_summary_group_map["Muon"].push_back("BeamAngleY");
    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "Birks");
    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "BetheBloch");
    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "Mass");
    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "NodeCutEff");

    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "TpiFromMichelRangeFit");
    mnv_plotter.error_summary_group_map["Pion Reconstruction"].push_back(
        "UntrackedPi");
  }
  if (is_subgroups) {
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "NonRESPi");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "MnvTunes");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "CCQE");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "RESPi");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "Coherent-Diffractive");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "DIS-Hadronization");
    mnv_plotter.error_summary_group_map["Cross Section Models"].push_back(
        "Elastic");

    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvn1pi");
    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvp1pi");
    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvn2pi");
    mnv_plotter.error_summary_group_map["NonRESPi"].push_back("GENIE_Rvp2pi");
    mnv_plotter.error_summary_group_map["MnvTunes"].push_back(
        "Low_Recoil_2p2h_Tune");
    mnv_plotter.error_summary_group_map["MnvTunes"].push_back("RPA_LowQ2");
    mnv_plotter.error_summary_group_map["MnvTunes"].push_back("RPA_HighQ2");
    mnv_plotter.error_summary_group_map["MnvTunes"].push_back("CCPi+ Tune");
    mnv_plotter.error_summary_group_map["CCQE"].push_back("GENIE_MaCCQE");
    mnv_plotter.error_summary_group_map["CCQE"].push_back(
        "GENIE_VecFFCCQEshape");
    mnv_plotter.error_summary_group_map["CCQE"].push_back(
        "GENIE_CCQEPauliSupViaKF");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_D2_MaRES");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_EP_MvRES");
    mnv_plotter.error_summary_group_map["RESPi"].push_back("GENIE_NormNCRES");
    mnv_plotter.error_summary_group_map["RESPi"].push_back(
        "GENIE_Theta_Delta2Npi");
    mnv_plotter.error_summary_group_map["RESPi"].push_back(
        "GENIE_D2_NormCCRES");
    mnv_plotter.error_summary_group_map["Coherent-Diffractive"].push_back(
        "DiffractiveModelUnc");
    mnv_plotter.error_summary_group_map["Coherent-Diffractive"].push_back(
        "CoherentPiUnc_CH");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back(
        "GENIE_AhtBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back(
        "GENIE_BhtBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back(
        "GENIE_CV1uBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back(
        "GENIE_CV2uBY");
    mnv_plotter.error_summary_group_map["DIS-Hadronization"].push_back(
        "GENIE_NormDISCC");
    mnv_plotter.error_summary_group_map["Elastic"].push_back("GENIE_MaNCEL");
    mnv_plotter.error_summary_group_map["Elastic"].push_back("GENIE_EtaNCEL");
  }
  //-- define colors of the standard errors
  mnv_plotter.error_color_map.clear();

  /*
    //Systematic color scheme
    mnv_plotter.error_color_map["Flux"]                   = kYellow-3;
    mnv_plotter.error_color_map["Genie_FSI"]              = kGreen+2;
    mnv_plotter.error_color_map["Genie_InteractionModel"] = kPink+2;
    mnv_plotter.error_color_map["Detector"]               = kCyan+2;
    mnv_plotter.error_color_map["RPA"]                    = kOrange+2;
    mnv_plotter.error_color_map["NonResPi"]               = kRed+2;
    mnv_plotter.error_color_map["Low_Recoil_2p2h_Tune"]   = kViolet+2;
  */
}

void AddingGroups(MnvPlotter& mnv_plotter, PlotUtils::MnvH1D& h,
                  std::string subgroup) {
  SetErrorGroups(mnv_plotter, true);
  std::vector<std::string> group =
      mnv_plotter.error_summary_group_map[subgroup];
  TH1D* ErrGroup = (TH1D*)h.Clone(Form("%s", subgroup.c_str()));
  ErrGroup->Reset();
  for (int i = 0; i < (int)group.size(); ++i) {
    TH1D* hErr = dynamic_cast<TH1D*>(
        h.GetVertErrorBand(group[i])->GetErrorBand(false, false).Clone(uniq()));
    MnvHist::AddInQuadrature(ErrGroup, hErr);
  }
  h.AddUncorrError(subgroup, ErrGroup, true);
}
void Plot_ErrorGroup(Plotter p, PlotUtils::MnvH1D* h,
                     std::string error_group_name, std::string tag,
                     double ignore_threshold = 0., double ymax = -1.) {
  TCanvas canvas("c1", "c1");

  p.m_mnv_plotter.good_colors = MnvColors::GetColors(MnvColors::k36Palette);
  p.m_mnv_plotter.legend_n_columns = 2;
  p.m_mnv_plotter.legend_offset_x = -0.03;
  if (error_group_name == "GENIE FSI") p.m_mnv_plotter.legend_n_columns = 2;
  // Clone hist
  PlotUtils::MnvH1D* hist = nullptr;

  if (p.m_variable->Name() == "q2") {
    hist = RebinQ2Plot(*h);
  } else if (p.m_variable->Name() == "enu" || p.m_variable->Name() == "pmu" ||
             p.m_variable->Name() == "ptmu" || p.m_variable->Name() == "pzmu") {
    hist = RebinningtoGeV(*h, p.m_variable->Name());
  } else {
    hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  }
  // X label
  p.SetXLabel(hist);

  // For groups MnvPlotter sets the ymax to headroom*axis_maximum_group
  // Or for summary, ymax = axis_maximum
  if (ymax != -1.) {
    p.m_mnv_plotter.axis_maximum = ymax;
    p.m_mnv_plotter.axis_maximum_group = ymax;
    p.m_mnv_plotter.headroom = 1.;
  }

  // const char* legend_position = error_group_name == "" ? "N" : "TR";
  const char* legend_position = "TL";

  if (error_group_name == "LEGENDONLY") {
    p.m_mnv_plotter.axis_maximum = 1000;
    p.m_mnv_plotter.axis_maximum_group = 1000;
    p.m_mnv_plotter.headroom = 1.;

    p.m_mnv_plotter.DrawErrorSummary(hist, legend_position, p.m_include_stat,
                                     true, ignore_threshold,
                                     p.m_do_cov_area_norm, "", p.m_do_frac_unc);
  } else {
    // XXX WARNING: problems when do_cov_area_norm = true !!
    p.m_mnv_plotter.DrawErrorSummary(
        hist, legend_position, p.m_include_stat, true, ignore_threshold,
        p.m_do_cov_area_norm, error_group_name, p.m_do_frac_unc, "", true);
  }

  // Plot Title (For systematics, this has to be after drawing)
  p.SetTitle(tag);

  std::string outfile_name =
      Form("ErrorSummary_%s_%s_%s_%s_%s_%s", tag.c_str(),
           p.m_variable->Name().c_str(), p.m_do_frac_unc_str.c_str(),
           p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(),
           error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "pdf");
}

// Deprecated?
void Plot_ErrorSummary(Plotter p, PlotUtils::MnvH1D* hist, std::string tag) {
  SetErrorGroups(p.m_mnv_plotter, false);

  //  Plot_ErrorGroup(p, hist, "", tag.c_str(), 0.0, 0.7);
  //  Plot_ErrorGroup(p, hist, "Flux", tag.c_str(), 0.0, 0.3);
  //  Plot_ErrorGroup(p, hist, "Detector", tag.c_str(), 0.02, 0.3);
  //  Plot_ErrorGroup(p, hist, "Genie_FSI", tag.c_str(), 0.04, 0.15);
  //  Plot_ErrorGroup(p, hist, "Genie_InteractionModel", tag.c_str(), 0.04,
  //  0.4); Plot_ErrorGroup(p, hist, "NonResPi", tag.c_str(), 0.0, 0.1);
  //  //  Plot_ErrorGroup(p, hist, "2p2h", tag.c_str(), 0.0, 0.1);
  //  //  Plot_ErrorGroup(p, hist, "RPA", tag.c_str(), 0.0, 0.1);
  //  //  Plot_ErrorGroup(p, hist, "Michel", tag.c_str(), 0.0, 0.3);
  //  // Plot_ErrorGroup(p, hist, "GENIE", tag.c_str(), 0.0, 0.3);
  //  //  Plot_ErrorGroup(p, hist, "Target", tag.c_str(), 0.0, 0.3);
  //  Plot_ErrorGroup(p, hist, "Others", tag.c_str(), 0.0, 0.3);
  //  Plot_ErrorGroup(p, hist, "Cross Section Models", tag.c_str(), 0.0, 0.3);
  //  //  Plot_ErrorGroup(p, hist, "PhysicsModel", tag.c_str(), 0.0, 0.3);

  Plot_ErrorGroup(p, hist, "", tag.c_str(), 0.0, 0.7);
  //  Plot_ErrorGroup(p, hist, "Detector", "CrossSection", 0.0, 0.1);
  Plot_ErrorGroup(p, hist, "Pion Reconstruction", tag.c_str(), 0.0, 0.3);
  Plot_ErrorGroup(p, hist, "Flux", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "GENIE FSI", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "Muon", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "Others", tag.c_str(), 0.02, 0.3);

  AddingGroups(p.m_mnv_plotter, *hist, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *hist, "RESPi");
  AddingGroups(p.m_mnv_plotter, *hist, "CCQE");
  AddingGroups(p.m_mnv_plotter, *hist, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *hist, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *hist, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *hist, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  Plot_ErrorGroup(p, hist, "Cross Section Models", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "RESPi", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "NonRESPi", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "CCQE", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "Coherent-Diffractive", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "DIS-Hadronization", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "MnvTunes", tag.c_str(), 0.02, 0.3);
  Plot_ErrorGroup(p, hist, "Elastic", tag.c_str(), 0.02, 0.3);
}

//==============================================================================
// Event Selection Plots
//==============================================================================
void PlotVar_Selection(Plotter p, double ymax = -1., bool do_log_scale = false,
                       bool do_bg = true, bool do_tuned_bg = false,
                       bool do_bin_width_norm = true) {
  std::cout << "Plotting Selection " << p.m_variable->Name() << std::endl;
  TCanvas canvas("c1", "c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_selection_data);
  assert(p.m_variable->m_hists.m_selection_mc.hist);

  PlotUtils::MnvH1D* mc = nullptr;
  PlotUtils::MnvH1D* data = nullptr;
  PlotUtils::MnvH1D* tmp_bg = nullptr;

  if (p.m_variable->Name() == "q2") {
    // if (false) {
    mc = RebinQ2Plot(*p.m_variable->m_hists.m_selection_mc.hist);
    data = RebinQ2Plot(*p.m_variable->m_hists.m_selection_data);
    if (do_bg) {
      if (do_tuned_bg)
        tmp_bg = RebinQ2Plot(*p.m_variable->m_hists.m_tuned_bg);
      else
        tmp_bg = RebinQ2Plot(*p.m_variable->m_hists.m_bg.hist);
    }

  } else if (p.m_variable->Name() == "enu" || p.m_variable->Name() == "pmu" ||
             p.m_variable->Name() == "ptmu" || p.m_variable->Name() == "pzmu") {
    data = RebinningtoGeV(*p.m_variable->m_hists.m_selection_data,
                          p.m_variable->Name());
    mc = RebinningtoGeV(*p.m_variable->m_hists.m_selection_mc.hist,
                        p.m_variable->Name());
    if (do_bg) {
      if (do_tuned_bg) {
        tmp_bg = RebinningtoGeV(*p.m_variable->m_hists.m_tuned_bg,
                                p.m_variable->Name());
      } else
        tmp_bg = RebinningtoGeV(*p.m_variable->m_hists.m_bg.hist,
                                p.m_variable->Name());
    }
  } else {
    data = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_selection_data->Clone(
        "data");
    mc = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_selection_mc.hist->Clone(
        "mc");
    if (do_bg) {
      if (do_tuned_bg)
        tmp_bg = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_tuned_bg->Clone(
            "bg_tmp");
      else
        tmp_bg = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg.hist->Clone(
            "bg_tmp");
    }
  }

  // Log Scale
  if (do_log_scale) {
    canvas.SetLogy();
    p.m_mnv_plotter.axis_minimum = 1;
  }
  if (p.m_variable->Name() == "q2") {
    canvas.SetLogx();
  }

  // Y-axis limit
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;

  // X label
  // p.SetXLabel(p.m_variable->m_hists.m_selection_mc.hist);
  p.SetXLabel(mc);
  // p.SetXLabel(data);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm)
    pot_scale = 1.;
  else
    pot_scale = p.m_data_pot / p.m_mc_pot;

  // Bin Width Normalization
  if (do_bin_width_norm) {
    if (tmp_bg) tmp_bg->Scale(1., "width");
    if (data) data->Scale(1., "width");
    mc->Scale(1., "width");

    if (p.m_variable->Name() == "q2") {
      data->SetBinContent(2, data->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      tmp_bg->SetBinContent(2,
                            tmp_bg->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      mc->SetBinContent(2, mc->GetBinContent(2) * (0.025 - 0.006) / 0.025);
    }

    // Y label
    std::string yaxis = "N Events / " + p.m_variable->m_units;
    mc->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = false;
  // note: this function applies a POT scale to the bg hist
  p.m_mnv_plotter.DrawDataMCWithErrorBand(
      data, mc, pot_scale, "TR", use_hist_titles, tmp_bg, NULL,
      p.m_do_cov_area_norm, p.m_include_stat);

  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);
  else
    p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);

  // Plot Title
  p.m_mnv_plotter.title_size = 0.05;
  p.SetTitle("Selection");

  std::string logy_str = do_log_scale ? "_logscale" : "";

  std::string bg_str = do_tuned_bg ? "_tunedBG" : "_untunedBG";

  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("Selection_%s_%s_%s%s%s%s", p.m_variable->Name().c_str(),
           p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bg_str.c_str(), bwn_str.c_str());
  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotVar_ErrorSummary(Plotter p) {
  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_selection_mc.hist);

  SetErrorGroups(p.m_mnv_plotter, false);

  PlotUtils::MnvH1D* sel =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_selection_mc.hist->Clone(
          uniq());
  //  Plot_ErrorGroup(p, sel, "", "Sel", 0.0, 0.35);
  //  Plot_ErrorGroup(p, sel, "LEGENDONLY", "Sel", 0.0, 0.2);
  //  //  Plot_ErrorGroup(p, sel, "2p2h", "Sel", 0.0, 0.01);
  //  Plot_ErrorGroup(p, sel, "Detector", "Sel", 0.0, 0.15);
  //  Plot_ErrorGroup(p, sel, "Flux", "Sel", 0.0, 0.15);
  //  Plot_ErrorGroup(p, sel, "Genie_FSI_nucleons", "Sel", 0.0, 0.06);
  //  Plot_ErrorGroup(p, sel, "Genie_FSI_pions", "Sel", 0.0, 0.2);
  //  Plot_ErrorGroup(p, sel, "Genie_InteractionModel", "Sel", 0.0, 0.2);
  //  Plot_ErrorGroup(p, sel, "Muon", "Sel", 0.0, 0.14);
  //  Plot_ErrorGroup(p, sel, "NonResPi", "Sel", 0.0, 0.08);
  //  //  Plot_ErrorGroup(p, sel, "RPA", "Sel", 0.0, 0.012);
  //  //  Plot_ErrorGroup(p, sel, "Michel", "Sel", 0.0, 0.15);
  //  //  Plot_ErrorGroup(p, sel, "GENIE", "Sel", 0.0, 0.30);
  //  //  Plot_ErrorGroup(p, sel, "Target", "Sel", 0.0, 0.15);
  //  Plot_ErrorGroup(p, sel, "Others", "Sel", 0.0, 0.05);
  //  Plot_ErrorGroup(p, sel, "Cross Section Models", "Sel", 0.0, 0.15);
  //  //  Plot_ErrorGroup(p, sel, "PhysicsModel", "Sel", 0.0, 0.15);
  //  Plot_ErrorGroup(p, sel, "LEGENDONLY", "Sel", 0.0, 0.1);
  Plot_ErrorGroup(p, sel, "", "Sel", 0.0, 0.4);
  //  Plot_ErrorGroup(p, sel, "Detector", "Sel", 0.0, 0.1);
  Plot_ErrorGroup(p, sel, "Pion Reconstruction", "Sel", 0.0, 0.1);
  Plot_ErrorGroup(p, sel, "Flux", "Sel", 0.0, 0.15);
  Plot_ErrorGroup(p, sel, "GENIE FSI", "Sel", 0.0, 0.2);
  Plot_ErrorGroup(p, sel, "Muon", "Sel", 0.0, 0.2);
  Plot_ErrorGroup(p, sel, "Others", "Sel", 0.0, 0.03);

  AddingGroups(p.m_mnv_plotter, *sel, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *sel, "RESPi");
  AddingGroups(p.m_mnv_plotter, *sel, "CCQE");
  AddingGroups(p.m_mnv_plotter, *sel, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *sel, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *sel, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *sel, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  Plot_ErrorGroup(p, sel, "Cross Section Models", "Sel", 0.0, 0.12);
  Plot_ErrorGroup(p, sel, "RESPi", "Sel", 0.0, 0.07);
  Plot_ErrorGroup(p, sel, "NonRESPi", "Sel", 0.0, 0.1);
  Plot_ErrorGroup(p, sel, "CCQE", "Sel", 0.0, 0.05);
  Plot_ErrorGroup(p, sel, "Coherent-Diffractive", "Sel", 0.0, 0.1);
  Plot_ErrorGroup(p, sel, "DIS-Hadronization", "Sel", 0.0, 0.02);
  Plot_ErrorGroup(p, sel, "MnvTunes", "Sel", 0.0, 0.01);
  Plot_ErrorGroup(p, sel, "Elastic", "Sel", 0.0, 0.01);
}

//==============================================================================
// Background-Subtracted
//==============================================================================
void Plot_BGSub(Plotter p, std::string outdir = ".", double ymax = -1,
                bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << "Plotting BG-subtracted Data " << p.m_variable->Name()
            << std::endl;

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_bg_subbed_data);
  assert(p.m_variable->m_hists.m_effnum.hist);

  TCanvas canvas("c1", "c1");

  // Get Hists
  PlotUtils::MnvH1D* tmp_bg_subbed_data = nullptr;
  PlotUtils::MnvH1D* tmp_effnum = nullptr;
  if (p.m_variable->Name() == "q2") {
    tmp_bg_subbed_data = RebinQ2Plot(*p.m_variable->m_hists.m_bg_subbed_data);
    tmp_effnum = RebinQ2Plot(*p.m_variable->m_hists.m_effnum.hist);
  } else if (p.m_variable->Name() == "enu" || p.m_variable->Name() == "pmu" ||
             p.m_variable->Name() == "ptmu" || p.m_variable->Name() == "pzmu") {
    tmp_bg_subbed_data = RebinningtoGeV(*p.m_variable->m_hists.m_bg_subbed_data,
                                        p.m_variable->Name());
    tmp_effnum = RebinningtoGeV(*p.m_variable->m_hists.m_effnum.hist,
                                p.m_variable->Name());
  } else {
    tmp_bg_subbed_data =
        (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg_subbed_data->Clone(
            "unfolded");
    tmp_effnum = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_effnum.hist->Clone(
        "effnum_true");
  }

  TH1D* bg_sub_data_w_tot_error =
      new TH1D(tmp_bg_subbed_data->GetCVHistoWithError());
  TH1D* bg_sub_data_w_stat_error =
      new TH1D(tmp_bg_subbed_data->GetCVHistoWithStatError());
  TH1D* effnum_w_stat_error = new TH1D(tmp_effnum->GetCVHistoWithStatError());

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

  // X label
  p.SetXLabel(bg_sub_data_w_tot_error);
  p.SetXLabel(bg_sub_data_w_stat_error);
  p.SetXLabel(effnum_w_stat_error);

  // Overall Normalization
  double pot_scale = -99.;
  // PlotUtils::MnvH1D* tmp_bg = nullptr;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  if (do_bin_width_norm) {
    bg_sub_data_w_tot_error->Scale(1., "width");
    bg_sub_data_w_stat_error->Scale(1., "width");
    effnum_w_stat_error->Scale(1., "width");

    if (p.m_variable->Name() == "q2") {
      bg_sub_data_w_tot_error->SetBinContent(
          2,
          bg_sub_data_w_tot_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      bg_sub_data_w_stat_error->SetBinContent(
          2,
          bg_sub_data_w_stat_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      effnum_w_stat_error->SetBinContent(
          2, effnum_w_stat_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
    }
    // Y label
    std::string yaxis = "N Events / " + p.m_variable->m_units;
    effnum_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = false;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_tot_error,
                                          effnum_w_stat_error, pot_scale, "TR",
                                          use_hist_titles);

  // p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_stat_error,
  // effnum_w_stat_error, pot_scale, "TR");
  // p.m_mnv_plotter.DrawDataMC(bg_sub_data_w_tot_error, effnum_w_stat_error,
  // pot_scale, "TR");
  // p.m_mnv_plotter.DrawDataMCWithErrorBand(bg_sub_data_w_error, effnum,
  // pot_scale, "TR",
  //                                        use_hist_titles, NULL, NULL,
  //                                        p.m_do_cov_area_norm,
  //                                        p.m_include_stat);

  // POT info
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);
  else
    p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);

  // Label re: error bars
  p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.7, 0.75, 0.03);

  // Plot Title
  p.m_mnv_plotter.title_size = 0.04;
  p.SetTitle(
      "Background Subtracted");  // + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/BGSubData_%s_%s_%s%s%s", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());
  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotBGSub_ErrorGroup(Plotter p, std::string error_group_name,
                          double ignore_threshold = 0., double ymax = -1.) {
  TCanvas canvas("c1", "c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_bg_subbed_data);

  // Get the histo
  PlotUtils::MnvH1D* bg_sub_data =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg_subbed_data->Clone(
          "bg_sub_data");

  // X label
  p.SetXLabel(bg_sub_data);

  // Y-axis limit
  // For groups MnvPlotter sets the ymax to headroom*axis_maximum_group
  // Or for summary, ymax = axis_maximum
  if (ymax != -1.) {
    p.m_mnv_plotter.axis_maximum = ymax;
    p.m_mnv_plotter.axis_maximum_group = ymax;
    p.m_mnv_plotter.headroom = 1.;
  }

  // Draw
  // XXX WARNING: potential problems when do_cov_area_norm = true
  p.m_mnv_plotter.DrawErrorSummary(bg_sub_data, "TR", p.m_include_stat, true,
                                   ignore_threshold, p.m_do_cov_area_norm,
                                   error_group_name, p.m_do_frac_unc);

  // Plot Title (For systematics, this has to be after drawing)
  p.SetTitle("Background-Subtracted Data");

  // Save to file
  std::string outfile_name = Form(
      "ErrorSummary_BGSubData_%s_%s_%s_%s_%s", p.m_variable->Name().c_str(),
      p.m_do_frac_unc_str.c_str(), p.m_do_cov_area_norm_str.c_str(),
      GetSignalFileTag(p.m_signal_definition).c_str(),
      error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotBGSub_ErrorSummary(Plotter p) {
  SetErrorGroups(p.m_mnv_plotter, false);

  PlotUtils::MnvH1D* bg_sub_data =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg_subbed_data->Clone(
          "bg_sub_data");

  // name, ignore threshold, ymax
  // PlotBGSub_ErrorGroup(p, "",                       0.0,  -1.);// // plot all
  // groups together PlotBGSub_ErrorGroup(p, "Flux",                   0.0,
  // 0.1);// PlotBGSub_ErrorGroup(p, "Detector",               0.0,   0.3);//
  // PlotBGSub_ErrorGroup(p, "Genie_FSI",              0.01,  0.1);//
  // PlotBGSub_ErrorGroup(p, "Genie_InteractionModel", 0.02,   0.2);//
  // PlotBGSub_ErrorGroup(p, "NonResPi",               0.0,   0.1);//
  // PlotBGSub_ErrorGroup(p, "2p2h",                   0.0,   0.1);//
  // PlotBGSub_ErrorGroup(p, "RPA",                    0.0,   0.1);//
  //  Plot_ErrorGroup(p, bg_sub_data, "LEGENDONLY", "BGSub", 0.0);
  //  Plot_ErrorGroup(p, bg_sub_data, "", "BGSub", 0.0,
  //                  -1.);  // // plot all groups together
  //  Plot_ErrorGroup(p, bg_sub_data, "Flux", "BGSub", 0.0, 0.1); //
  //  Plot_ErrorGroup(p, bg_sub_data, "Detector", "BGSub", 0.0, 0.3); //
  //  Plot_ErrorGroup(p, bg_sub_data, "Genie_FSI_nucleons", "BGSub", 0.0, 0.1);
  //  // Plot_ErrorGroup(p, bg_sub_data, "Genie_FSI_pions", "BGSub", 0.0, 0.1);
  //  // Plot_ErrorGroup(p, bg_sub_data, "Genie_InteractionModel", "BGSub", 0.0,
  //                  0.2);                                            //
  //  Plot_ErrorGroup(p, bg_sub_data, "NonResPi", "BGSub", 0.0, 0.1);  //
  //  //  Plot_ErrorGroup(p, bg_sub_data, "2p2h", "BGSub", 0.0, 0.02);     //
  //  //  Plot_ErrorGroup(p, bg_sub_data, "RPA", "BGSub", 0.0, 0.1);       //
  //  //  Plot_ErrorGroup(p, bg_sub_data, "Michel", "BGSub", 0.0, 0.3);
  //  //  Plot_ErrorGroup(p, bg_sub_data, "GENIE", "BGSub", 0.0, 0.3);
  //  //  Plot_ErrorGroup(p, bg_sub_data, "Target", "BGSub", 0.0, 0.3);
  //  Plot_ErrorGroup(p, bg_sub_data, "Others", "BGSub", 0.0, 0.3);
  //  Plot_ErrorGroup(p, bg_sub_data, "Cross Section Models", "BGSub", 0.0,
  //  0.05);
  //  //  Plot_ErrorGroup(p, bg_sub_data, "PhysicsModel", "BGSub", 0.0, 0.02);
  //  Plot_ErrorGroup(p, bg_sub_data, "Muon", "BGSub", 0.0, 0.3);

  //  Plot_ErrorGroup(p, bg_sub_data, "LEGENDONLY", "BGSub", 0.0, 0.1);
  Plot_ErrorGroup(p, bg_sub_data, "", "BGSub", 0.0, 0.3);
  //  Plot_ErrorGroup(p, bg_sub_data, "Detector", "BGSub", 0.0, 0.1);
  Plot_ErrorGroup(p, bg_sub_data, "Pion Reconstruction", "BGSub", 0.0, 0.1);
  Plot_ErrorGroup(p, bg_sub_data, "Flux", "BGSub", 0.0, 0.2);
  Plot_ErrorGroup(p, bg_sub_data, "GENIE FSI", "BGSub", 0.0, 0.1);
  Plot_ErrorGroup(p, bg_sub_data, "Muon", "BGSub", 0.0, 0.3);
  Plot_ErrorGroup(p, bg_sub_data, "Others", "BGSub", 0.0, 0.20);

  AddingGroups(p.m_mnv_plotter, *bg_sub_data, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *bg_sub_data, "RESPi");
  AddingGroups(p.m_mnv_plotter, *bg_sub_data, "CCQE");
  AddingGroups(p.m_mnv_plotter, *bg_sub_data, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *bg_sub_data, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *bg_sub_data, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *bg_sub_data, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  Plot_ErrorGroup(p, bg_sub_data, "Cross Section Models", "BGSub", 0.0, 0.25);
  Plot_ErrorGroup(p, bg_sub_data, "RESPi", "BGSub", 0.0, 0.25);
  Plot_ErrorGroup(p, bg_sub_data, "NonRESPi", "BGSub", 0.0, 0.1);
  Plot_ErrorGroup(p, bg_sub_data, "CCQE", "BGSub", 0.0, 0.2);
  Plot_ErrorGroup(p, bg_sub_data, "Coherent-Diffractive", "BGSub", 0.0, 0.05);
  Plot_ErrorGroup(p, bg_sub_data, "DIS-Hadronization", "BGSub", 0.0);
  Plot_ErrorGroup(p, bg_sub_data, "MnvTunes", "BGSub", 0.0, 0.05);
  Plot_ErrorGroup(p, bg_sub_data, "Elastic", "BGSub", 0.0);
}

//==============================================================================
// Unfolded
//==============================================================================
void Plot_Unfolded(Plotter p, MnvH1D* data, MnvH1D* mc,
                   std::string outdir = ".", double ymax = -1,
                   bool do_log_scale = false, bool do_bin_width_norm = true) {
  std::cout << "Plotting Unfolded " << p.m_variable->Name() << std::endl;

  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  TCanvas canvas("c1", "c1");

  // Get Hists
  PlotUtils::MnvH1D* unfolded = nullptr;
  PlotUtils::MnvH1D* effnum_true = nullptr;
  if (p.m_variable->Name() == "q2") {
    unfolded = RebinQ2Plot(*data);
    effnum_true = RebinQ2Plot(*mc);
  } else if (p.m_variable->Name() == "enu" || p.m_variable->Name() == "pmu" ||
             p.m_variable->Name() == "ptmu" || p.m_variable->Name() == "pzmu") {
    unfolded = RebinningtoGeV(*data, p.m_variable->Name());
    effnum_true = RebinningtoGeV(*mc, p.m_variable->Name());
  } else {
    unfolded = (PlotUtils::MnvH1D*)data->Clone("unfolded");
    effnum_true = (PlotUtils::MnvH1D*)mc->Clone("effnum_true");
  }

  TH1D* unfolded_w_tot_error = new TH1D(unfolded->GetCVHistoWithError());
  TH1D* unfolded_w_stat_error = new TH1D(unfolded->GetCVHistoWithStatError());
  TH1D* effnum_true_w_stat_error =
      new TH1D(effnum_true->GetCVHistoWithStatError());

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

  // X label
  p.SetXLabel(unfolded_w_tot_error);
  p.SetXLabel(unfolded_w_stat_error);
  p.SetXLabel(effnum_true_w_stat_error);

  // Overall Normalization
  double pot_scale = -99.;
  if (p.m_do_cov_area_norm) {
    pot_scale = 1.;
  } else {
    pot_scale = p.m_data_pot / p.m_mc_pot;
  }

  // Bin Width Normalization
  if (do_bin_width_norm) {
    unfolded_w_tot_error->Scale(1., "width");
    unfolded_w_stat_error->Scale(1., "width");
    effnum_true_w_stat_error->Scale(1., "width");

    if (p.m_variable->Name() == "q2") {
      unfolded_w_tot_error->SetBinContent(
          2, unfolded_w_tot_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      unfolded_w_stat_error->SetBinContent(
          2, unfolded_w_stat_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      effnum_true_w_stat_error->SetBinContent(
          2,
          effnum_true_w_stat_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
    }
    // Y label
    std::string yaxis = "N Events / " + p.m_variable->m_units;
    effnum_true_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }

  // Draw
  const bool use_hist_titles = false;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(unfolded_w_tot_error,
                                          effnum_true_w_stat_error, pot_scale,
                                          "TR", use_hist_titles);

  // POT info
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);
  else
    p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, p.m_mc_pot, 0.3, 0.88);

  // Label re: error bars
  p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.7, 0.75, 0.03);

  // Plot Title
  p.m_mnv_plotter.title_size = 0.05;
  p.SetTitle("Unfolded ");  // + GetSignalName(p.m_signal_definition));

  // Print .png
  std::string logy_str = do_log_scale ? "_logscale" : "";
  std::string bwn_str = do_bin_width_norm ? "_BWN" : "";

  std::string outfile_name =
      Form("%s/Unfolded_%s_%s_%s%s%s", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotUnfolded_ErrorSummary(Plotter p) {
  SetErrorGroups(p.m_mnv_plotter, false);
  PlotUtils::MnvH1D* unf =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_unfolded->Clone(uniq());
  //  Plot_ErrorGroup(p, unf, "LEGENDONLY", "Unfolded", 0.0);
  //  Plot_ErrorGroup(p, unf, "", "Unfolded", 0.0);
  //  //  Plot_ErrorGroup(p, unf, "2p2h", "Unfolded", 0.0, 0.02);
  //  Plot_ErrorGroup(p, unf, "Detector", "Unfolded", 0.0, 0.1);
  //  Plot_ErrorGroup(p, unf, "Flux", "Unfolded", 0.0, 0.06);
  //  Plot_ErrorGroup(p, unf, "Genie_FSI_nucleons", "Unfolded", 0.0, 0.1);
  //  Plot_ErrorGroup(p, unf, "Genie_FSI_pions", "Unfolded", 0.0, 0.1);
  //  Plot_ErrorGroup(p, unf, "Genie_InteractionModel", "Unfolded", 0.0, 0.15);
  //  Plot_ErrorGroup(p, unf, "NonResPi", "Unfolded", 0.0, 0.1);
  //  //  Plot_ErrorGroup(p, unf, "RPA", "Unfolded", 0.0, 0.02);
  //  //  Plot_ErrorGroup(p, unf, "Michel", "Unfolded", 0.0, 0.1);
  //  //  Plot_ErrorGroup(p, unf, "GENIE", "Unfolded", 0.0, 0.26);
  //  //  Plot_ErrorGroup(p, unf, "Target", "Unfolded", 0.0, 0.1);
  //  Plot_ErrorGroup(p, unf, "Others", "Unfolded", 0.0, 0.24);
  //  Plot_ErrorGroup(p, unf, "Cross Section Models", "Unfolded", 0.0, 0.02);
  //  //  Plot_ErrorGroup(p, unf, "PhysicsModel", "Unfolded", 0.0, 0.02);
  //  Plot_ErrorGroup(p, unf, "Muon", "Unfolded", 0.0, 0.2);

  //  Plot_ErrorGroup(p, unf, "LEGENDONLY", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "", "Unfolded", 0.0, 0.4);
  //  Plot_ErrorGroup(p, unf, "Detector", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "Pion Reconstruction", "Unfolded", 0.0, 0.08);
  Plot_ErrorGroup(p, unf, "Flux", "Unfolded", 0.0, 0.15);
  Plot_ErrorGroup(p, unf, "GENIE FSI", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "Muon", "Unfolded", 0.0, 0.25);
  Plot_ErrorGroup(p, unf, "Others", "Unfolded", 0.0, 0.15);

  AddingGroups(p.m_mnv_plotter, *unf, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *unf, "RESPi");
  AddingGroups(p.m_mnv_plotter, *unf, "CCQE");
  AddingGroups(p.m_mnv_plotter, *unf, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *unf, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *unf, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *unf, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  Plot_ErrorGroup(p, unf, "Cross Section Models", "Unfolded", 0.0, 0.2);
  Plot_ErrorGroup(p, unf, "RESPi", "Unfolded", 0.0, 0.15);
  Plot_ErrorGroup(p, unf, "NonRESPi", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "CCQE", "Unfolded", 0.0, 0.2);
  Plot_ErrorGroup(p, unf, "Coherent-Diffractive", "Unfolded", 0.0, 0.05);
  Plot_ErrorGroup(p, unf, "DIS-Hadronization", "Unfolded", 0.0, 0.1);
  Plot_ErrorGroup(p, unf, "MnvTunes", "Unfolded", 0.0, 0.03);
  Plot_ErrorGroup(p, unf, "Elastic", "Unfolded", 0.0);
}

//==============================================================================
// Cross Section
//==============================================================================
void Plot_CrossSection(Plotter p, MnvH1D* data, MnvH1D* mc,
                       std::string outdir = ".", double ymax = -1,
                       bool do_log_scale = false,
                       bool do_bin_width_norm = true) {
  std::cout << "Plotting CrossSection " << p.m_variable->Name() << std::endl;
  std::string name = p.m_variable->Name();
  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  double uncfactor;
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

  PlotUtils::MnvH1D* data_xsec = nullptr;
  PlotUtils::MnvH1D* mc_xsec = nullptr;

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

  TCanvas canvas("c1", "c1");

  // Get Hists
  TH1D* data_xsec_w_tot_error = new TH1D(data_xsec->GetCVHistoWithError());
  TH1D* data_xsec_w_stat_error = new TH1D(data_xsec->GetCVHistoWithStatError());
  TH1D* mc_xsec_w_stat_error = new TH1D(mc_xsec->GetCVHistoWithStatError());
  /*
    if (p.m_variable->Name() == "pmu"){
      data_xsec_w_tot_error->GetXaxis()->SetRange(0., 20.);
      data_xsec_w_stat_error->GetXaxis()->SetRange(0., 20.);
      mc_xsec_w_stat_error->GetXaxis()->SetRange(0., 20.);
    }*/
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
  p.m_mnv_plotter.axis_title_offset_y = 1.5;

  // X label
  p.SetXLabel(data_xsec_w_tot_error);
  p.SetXLabel(data_xsec_w_stat_error);
  p.SetXLabel(mc_xsec_w_stat_error);

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
    mc_xsec_w_stat_error->Scale(1.e42, "width");

    // already divided by 0.25 - 0.006
    // but we really want to divide by 0.25
    // TODO sketchy AF
    if (p.m_variable->Name() == "q2") {
      data_xsec_w_tot_error->SetBinContent(
          2, data_xsec_w_tot_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      data_xsec_w_stat_error->SetBinContent(
          2,
          data_xsec_w_stat_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
      mc_xsec_w_stat_error->SetBinContent(
          2, mc_xsec_w_stat_error->GetBinContent(2) * (0.025 - 0.006) / 0.025);
    }
    // if (this is q2) scale the first bin differently

    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel +
                        " (10^{-42} cm^{2}/" + p.m_variable->m_units +
                        "/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    mc_xsec_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }

  else {
    // data_xsec_w_tot_error ->Scale(1.e38, "width");
    // data_xsec_w_stat_error->Scale(1.e38, "width");
    // mc_xsec_w_stat_error  ->Scale(1.e38, "width");

    data_xsec_w_tot_error->Scale(1.e42);
    data_xsec_w_stat_error->Scale(1.e42);
    mc_xsec_w_stat_error->Scale(1.e42);

    // Y label
    // std::string yaxis = "d#sigma/d" + p.m_variable->m_hists.m_xlabel + "
    // (10^{-38} cm^{2}/" + p.m_variable->m_units + "/nucleon)";
    std::string yaxis = "#sigma (10^{-42} cm^{2}/nucleon)";
    p.m_mnv_plotter.axis_title_size_y = 0.04;
    mc_xsec_w_stat_error->GetYaxis()->SetTitle(yaxis.c_str());
  }
  /*
  // Print xsec and error for each bin (AFTER BWN)
  int low_edge = -99;
  int up_edge = -99;
  double val = -99.;
  double err = -99.;
  double frac_err = -99.;
  for (int i = 0; i <= data_xsec_w_tot_error->GetNbinsX(); ++i ){
    low_edge = data_xsec_w_tot_error->GetXaxis()->GetBinLowEdge(i);
    up_edge  = data_xsec_w_tot_error->GetXaxis()->GetBinUpEdge(i);
    val      = data_xsec_w_tot_error->GetBinContent(i);
    err      = data_xsec_w_tot_error->GetBinError(i);
    frac_err = err/val;

    std::cout << i << "  " << low_edge << "  " << up_edge << "  " << val << "  "
  << err << "  " << frac_err << "\n";
  }
  */

  // Draw
  const bool use_hist_titles = false;
  p.m_mnv_plotter.DrawDataMCWithErrorBand(data_xsec_w_tot_error,
                                          mc_xsec_w_stat_error, pot_scale, "TR",
                                          use_hist_titles);

  data_xsec_w_stat_error->SetMarkerStyle(20);
  data_xsec_w_stat_error->SetMarkerSize(1.0);
  data_xsec_w_stat_error->SetMarkerColor(1);
  data_xsec_w_stat_error->SetLineWidth(1);
  data_xsec_w_stat_error->SetLineStyle(1);
  data_xsec_w_stat_error->SetLineColor(1);
  data_xsec_w_stat_error->DrawCopy("SAME E1 X0");
  // Add chi2 label
  {
    const bool use_data_error_mtx = true;
    const bool use_only_shape_errors = false;
    const bool use_model_stat =
        false;  // model statistical errors -- these are small if you use a
                // priori effden, but atm I'm not. So keep this off.
    //    p.m_mnv_plotter.AddChi2Label(data_xsec, mc_xsec, pot_scale, "TR",
    //    0.03, -0.175, use_data_error_mtx, use_only_shape_errors); // this auto
    //    turns on
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
    /*    if (p.m_variable->Name() == "tpi") ndf = 6;
        if (p.m_variable->Name() == "enu") ndf = 6;
        if (p.m_variable->Name() == "pzmu") ndf = 9;
        if (p.m_variable->Name() == "pmu") ndf = 8;
        if (p.m_variable->Name() == "wexp") ndf = 4;*/
    char* words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf,
                       chi2 / (Double_t)ndf);
    int align = 33;
    p.m_mnv_plotter.AddPlotLabel(words, 0.8, 0.745, 0.03, 1, 62, align);
  }

  // POT info
  // -1 --> don't do mc POT
  if (p.m_do_cov_area_norm)
    p.m_mnv_plotter.AddAreaNormBox(p.m_data_pot, -1, 0.3, 0.88);
  else {
    //  p.m_mnv_plotter.AddPOTNormBox(p.m_data_pot, 0, 0.3, 0.88);
    p.m_mnv_plotter.WriteNorm("POT-Normalized", 0.3, 0.88, 0.03);
    p.m_mnv_plotter.WriteNorm(Form("Data POT: %.2E", p.m_data_pot), 0.3,
                              0.88 - 0.03, 0.03);
  }
  //  p.m_mnv_plotter.WritePreliminary(0.33, 0.812);
  // Label re: error bars
  p.m_mnv_plotter.AddPlotLabel("MC stat-only errors", 0.7, 0.772, 0.03);

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
      Form("%s/CrossSection_%s_%s_%s%s%s", outdir.c_str(),
           p.m_variable->Name().c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(), logy_str.c_str(),
           bwn_str.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}

void PlotCrossSection_ErrorSummary(Plotter p) {
  SetErrorGroups(p.m_mnv_plotter, false);
  PlotUtils::MnvH1D* xsec =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_cross_section->Clone(uniq());
  std::string name = p.m_variable->Name();

  double uncfactor;
  if (name == "mixtpi") {
    // uncfactor = 5.7;// Without CCPi weight
    uncfactor = 6.5;  // Wight CCPi Weight
    xsec->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "mixthetapi_deg") {
    // uncfactor = 10.2;// Without CCPi weight
    uncfactor = 6.2;  // Wight CCPi Weight
    xsec->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "q2") {
    // uncfactor = 6.9;// Without CCPi weight
    uncfactor = 7.5;  // Wight CCPi Weight
    xsec->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "ptmu") {
    // uncfactor = 7.9;// Without CCPi weight
    uncfactor = 7.3;  // Wight CCPi Weight
    xsec->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }
  if (name == "thetamu_deg") {
    uncfactor = 15.2;  // Wight CCPi Weight
    xsec->ModifyStatisticalUnc(uncfactor,
                               Form("unfolding_cov_matrix_%s", name.c_str()));
  }

  // for (auto b : xsec->GetErrorBandNames()) std::cout << b << "\n";

  double detector_threshold = 0.0, detector_ymax = .15;
  double total_max = 0.3;
  double FSI_threshold = 0.0, FSI_ymax = 0.1;
  double Int_threshold = 0.015, Int_ymax = 0.15;
  if (name == "enu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.0085;
    FSI_ymax = 0.045;
    Int_ymax = 0.08;
  } else if (name == "pmu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.011;
    FSI_ymax = 0.045;
    Int_ymax = 0.08;
  } else if (name == "ptmu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.015;
    FSI_ymax = 0.045;
    Int_threshold = 0.025;
    Int_ymax = 0.2;
    total_max = 0.5;
  } else if (name == "pzmu") {
    detector_ymax = 0.2;
    FSI_threshold = 0.015;
    FSI_ymax = 0.045;
    Int_ymax = 0.09;
  } else if (name == "q2") {
    detector_ymax = 0.25;
    FSI_threshold = 0.015;
    FSI_ymax = 0.065;
    Int_threshold = 0.025;
    Int_ymax = 0.3;
    total_max = 0.5;
  } else if (name == "thetamu_deg") {
    detector_ymax = 0.25;
    FSI_threshold = 0.015;
    FSI_ymax = 0.065;
    Int_threshold = 0.025;
    Int_ymax = 0.15;
  } else if (name == "thetapi_deg") {
    FSI_threshold = 0.015;
    FSI_ymax = 0.065;
    Int_threshold = 0.02;
    Int_ymax = 0.12;
  } else if (name == "tpi") {
    detector_ymax = 0.3;
    FSI_threshold = 0.01;
    FSI_ymax = 0.085;
    Int_ymax = 0.12;
  } else if (name == "wexp") {
    detector_ymax = 0.6;
    FSI_threshold = 0.008;
    FSI_ymax = 0.06;
    Int_ymax = 0.1;
  }

  const double ymax = 0.5;
  p.m_include_stat = true;
  // name, ignore threshold, ymax
  //  Plot_ErrorGroup(p, xsec, "LEGENDONLY", "CrossSection", 0.0, 0.1);
  Plot_ErrorGroup(p, xsec, "", "CrossSection", 0.0, total_max);
  //  Plot_ErrorGroup(p, xsec, "Detector", "CrossSection", 0.0, 0.1);
  Plot_ErrorGroup(p, xsec, "Pion Reconstruction", "CrossSection", 0.0, 0.8);
  Plot_ErrorGroup(p, xsec, "Flux", "CrossSection", 0.0, 0.2);
  Plot_ErrorGroup(p, xsec, "GENIE FSI", "CrossSection", 0.0, 0.1);
  Plot_ErrorGroup(p, xsec, "Muon", "CrossSection", 0.0, 0.3);
  Plot_ErrorGroup(p, xsec, "Others", "CrossSection", 0.0, 0.15);

  AddingGroups(p.m_mnv_plotter, *xsec, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *xsec, "RESPi");
  AddingGroups(p.m_mnv_plotter, *xsec, "CCQE");
  AddingGroups(p.m_mnv_plotter, *xsec, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *xsec, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *xsec, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *xsec, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  Plot_ErrorGroup(p, xsec, "Cross Section Models", "CrossSection", 0.0, 0.25);
  Plot_ErrorGroup(p, xsec, "RESPi", "CrossSection", 0.0, 0.15);
  Plot_ErrorGroup(p, xsec, "NonRESPi", "CrossSection", 0.0, 0.1);
  Plot_ErrorGroup(p, xsec, "CCQE", "CrossSection", 0.0, 0.2);
  Plot_ErrorGroup(p, xsec, "Coherent-Diffractive", "CrossSection", 0.0, 0.05);
  Plot_ErrorGroup(p, xsec, "DIS-Hadronization", "CrossSection", 0.0, 0.05);
  Plot_ErrorGroup(p, xsec, "MnvTunes", "CrossSection", 0.0, 0.3);
  Plot_ErrorGroup(p, xsec, "Elastic", "CrossSection", 0.0);
}

void PlotMatrix(TMatrixD mtx, std::string name, std::string tag) {
  myPlotStyle();
  setCorrelationPalette();
  TCanvas* c = new TCanvas(uniq(), "scaled matrix");
  mtx.Draw("colz");
  c->Print(Form("%s_%s.png", name.c_str(), tag.c_str()));
  // c->Print(Form("CovMatrix_%s_%s.png", name.c_str(), tag.c_str()));
}

void PrintChi2Info(Plotter p, MnvH1D* data, MnvH1D* mc) {
  // Make sure we remembered to load the source histos from the input file.
  assert(data);
  assert(mc);

  PlotUtils::MnvH1D* data_xsec = (PlotUtils::MnvH1D*)data->Clone("data");
  PlotUtils::MnvH1D* mc_xsec = (PlotUtils::MnvH1D*)mc->Clone("mc");

  TCanvas canvas("c1", "c1");

  myPlotStyle();
  setCorrelationPalette();

  const bool use_data_error_mtx = true;
  const bool use_only_shape_errors = false;
  bool use_model_stat = false;  // model statistical errors, should be small

  int ndf = -1;
  double pot_scale = 1.;
  double chi2 = p.m_mnv_plotter.Chi2DataMC(
      data_xsec, mc_xsec, ndf, pot_scale, use_data_error_mtx,
      use_only_shape_errors, use_model_stat);
  use_model_stat = true;
  double chi2_mcstat = p.m_mnv_plotter.Chi2DataMC(
      data_xsec, mc_xsec, ndf, pot_scale, use_data_error_mtx,
      use_only_shape_errors, use_model_stat);
  // std::cout << p.m_variable->Name() << "\n";
  // std::cout << "   chi2 = "         << chi2     << "\n";
  // std::cout << "   ndf = "          << ndf      << "\n";
  //"   chi2/ndf = "     <<
  std::cout << p.m_variable->Name() << "  " << chi2 / ndf << "  "
            << chi2_mcstat / ndf << "\n";
}

//==============================================================================
// W Sideband Fit
//==============================================================================
void PlotWSidebandFit_ErrorGroup(Plotter p, std::string error_group_name,
                                 PlotUtils::MnvH1D* h, std::string tag) {
  TCanvas canvas("c1", "c1");
  p.m_mnv_plotter.good_colors = MnvColors::GetColors(MnvColors::k36Palette);
  // Make sure we remembered to load the source histos from the input file.
  assert(h);

  // Get the hist
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");

  // if (error_group_name != "")
  p.m_mnv_plotter.axis_maximum = 0.6;

  // X label
  // hist->GetXaxis()->SetTitle();

  double ignore_threshold;
  if (error_group_name == "")
    ignore_threshold = 0.;
  else
    ignore_threshold = 0.0;

  // XXX WARNING: potential problems when do_cov_area_norm = true
  p.m_mnv_plotter.DrawErrorSummary(hist, "TR", p.m_include_stat, true,
                                   ignore_threshold, p.m_do_cov_area_norm,
                                   error_group_name, p.m_do_frac_unc);

  // Plot Title (For systematics, this has to be after drawing)
  p.SetTitle("W Sideband Fit");

  std::string outfile_name =
      Form("ErrorSummary_%s_%s_%s_%s_%s", tag.c_str(),
           p.m_do_frac_unc_str.c_str(), p.m_do_cov_area_norm_str.c_str(),
           GetSignalFileTag(p.m_signal_definition).c_str(),
           error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");

  p.m_mnv_plotter.axis_maximum = 0.6;
}

void PlotWSidebandFit_ErrorSummary(Plotter p, PlotUtils::MnvH1D* hist,
                                   std::string tag) {
  SetErrorGroups(p.m_mnv_plotter, false);

  //  PlotWSidebandFit_ErrorGroup(p, "LEGENDONLY", hist, tag);  // plot all
  //  groups together
  //  PlotWSidebandFit_ErrorGroup(p, "", hist, tag);  // plot all groups
  //  together PlotWSidebandFit_ErrorGroup(p, "Flux", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Detector", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Genie_FSI_nucleons", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Genie_FSI_pions", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Genie_InteractionModel", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "NonResPi", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "2p2h", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "RPA", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Michel", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "GENIE", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Target", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Others", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Cross Section Models", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "PhysicsModel", hist, tag);
  //  PlotWSidebandFit_ErrorGroup(p, "Muon", hist, tag);

  PlotWSidebandFit_ErrorGroup(p, "", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Pion Reconstruction", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Flux", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "GENIE FSI", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Muon", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Others", hist, tag);

  AddingGroups(p.m_mnv_plotter, *hist, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *hist, "RESPi");
  AddingGroups(p.m_mnv_plotter, *hist, "CCQE");
  AddingGroups(p.m_mnv_plotter, *hist, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *hist, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *hist, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *hist, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  PlotWSidebandFit_ErrorGroup(p, "Cross Section Models", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "RESPi", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "NonRESPi", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "CCQE", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Coherent-Diffractive", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "DIS-Hadronization", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "MnvTunes", hist, tag);
  PlotWSidebandFit_ErrorGroup(p, "Elastic", hist, tag);
}

void PlotWSidebandStacked(
    const Variable* variable, const PlotUtils::MnvH1D* h_data,
    PlotUtils::MnvH1D* loW, PlotUtils::MnvH1D* midW, PlotUtils::MnvH1D* hiW,
    const TObjArray& array_mc, float data_pot, float mc_pot,
    SignalDefinition signal_definition, std::string tag = "", double ymax = -1,
    bool do_bin_width_norm = true, bool do_postfit = false) {
  // Never don't clone when plotting
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");
  //  TObjArray* array = new TObjArray();
  array.RemoveRange(0, 3);
  array.Compress();
  std::cout << "Plotting " << variable->Name() << std::endl;

  PlotUtils::MnvH1D* data = nullptr;
  PlotUtils::MnvH1D* h_low_W = nullptr;
  PlotUtils::MnvH1D* h_mid_W = nullptr;
  PlotUtils::MnvH1D* h_hi_W = nullptr;

  if (variable->Name() == "q2") {
    // if (false) {
    data = RebinQ2Plot(*h_data);
    for (int i = 0; i <= 3; ++i) {
      array[i] = RebinQ2Plot(*dynamic_cast<PlotUtils::MnvH1D*>(array[i]));
    }
  } else if (variable->Name() == "enu" || variable->Name() == "pmu" ||
             variable->Name() == "ptmu" || variable->Name() == "pzmu" ||
             variable->Name() == "wexp" || variable->Name() == "wexp_fit") {
    data = RebinningtoGeV(*h_data, variable->Name());
    for (int i = 0; i <= 3; ++i) {
      array[i] = RebinningtoGeV(*dynamic_cast<PlotUtils::MnvH1D*>(array[i]),
                                variable->Name());
    }
  } else {
    data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  }

  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  if (ymax > 0) mnvPlotter.axis_maximum = ymax;

  double pot_scale = data_pot / mc_pot;
  std::string fit_str = do_postfit ? "Post" : "Pre";
  std::string label =
      Form("Breakdown_%sWSideband_%s_%s_PN_%s", fit_str.c_str(),
           variable->m_label.c_str(),
           GetSignalFileTag(signal_definition).c_str(), tag.c_str());
  TCanvas cE("c1", "c1");

  std::string y_label = "Events";
  // bin width norm
  if (do_bin_width_norm) {
    data->Scale(1., "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1., "width");
    //      array->Add(dynamic_cast<PlotUtils::MnvH1D*>(h));
    y_label = "Events / " + variable->m_units;
  }
  if (do_postfit) {
    dynamic_cast<PlotUtils::MnvH1D*>(array[1])->Scale(loW->GetBinContent(1));
    dynamic_cast<PlotUtils::MnvH1D*>(array[2])->Scale(midW->GetBinContent(1));
    dynamic_cast<PlotUtils::MnvH1D*>(array[3])->Scale(hiW->GetBinContent(1));
  }

  mnvPlotter.DrawDataStackedMC(
      data, &array, pot_scale, "TR", "Data", -1, -1, 1001,
      Form("%s (%s)", variable->m_hists.m_xlabel.c_str(),
           variable->m_units.c_str()),
      y_label.c_str());

  double arrow_height = data->GetBinContent(data->GetMaximumBin()) * 1.1;
  double arrow_location = signal_definition.m_w_max / 1000 + 0.1;
  mnvPlotter.AddCutArrow(arrow_location, 0.0, arrow_height, 0.2, "R");
  //  mnvPlotter.WritePreliminary("TL");
  // mnvPlotter.AddPOTNormBox(data_pot, 1.0, 300, 0.85);
  mnvPlotter.WriteNorm("POT-Normalized", 0.3, 0.88, 0.03);
  mnvPlotter.WriteNorm(Form("Data POT: %.2E", data_pot), 0.3, 0.88 - 0.03,
                       0.03);
  mnvPlotter.AddHistoTitle(tag.c_str());
  mnvPlotter.MultiPrint(&cE, label, "pdf");
}

void PlotFittedW(const Variable* variable, const CVUniverse& universe,
                 const PlotUtils::MnvH1D* loW_fit,
                 const PlotUtils::MnvH1D* midW_fit,
                 const PlotUtils::MnvH1D* hiW_fit, float data_pot, float mc_pot,
                 SignalDefinition signal_definition, bool do_prefit = false,
                 std::string tag = "", double ymax = -1,
                 bool do_bin_width_norm = true) {
  // Setup
  std::cout << "Plotting " << variable->Name() << std::endl;
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  double pot_scale = data_pot / mc_pot;
  std::string fit_str = do_prefit ? "Pre" : "Post";
  std::string label =
      Form("%sWFit_%s_%s_PN_%s", fit_str.c_str(), variable->m_label.c_str(),
           GetSignalFileTag(signal_definition).c_str(), tag.c_str());
  TCanvas cE("c1", "c1");

  // Never don't clone when plotting
  PlotUtils::MnvH1D* h_data =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_data->Clone("data");
  PlotUtils::MnvH1D* h_sig =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_sig
          .univHist(&universe)
          ->Clone("sig");
  PlotUtils::MnvH1D* h_loW =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_loW
          .univHist(&universe)
          ->Clone("loW");
  PlotUtils::MnvH1D* h_midW =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_midW
          .univHist(&universe)
          ->Clone("midW");
  PlotUtils::MnvH1D* h_hiW =
      (PlotUtils::MnvH1D*)variable->m_hists.m_wsidebandfit_hiW
          .univHist(&universe)
          ->Clone("hiW");

  // Apply fit
  if (do_prefit) {
    ;
  } else {
    h_loW->Scale(loW_fit->GetBinContent(1));
    h_midW->Scale(midW_fit->GetBinContent(1));
    h_hiW->Scale(hiW_fit->GetBinContent(1));
    // std::cout << "low fit = " << loW_fit->GetBinContent(1) << "\n";
    // std::cout << "Middle fit = " << loW_fit->GetBinContent(1) << "\n";
    // std::cout << "High fit = " << loW_fit->GetBinContent(1) << "\n";
  }

  std::string y_label = "Events";
  // bin width norm
  if (do_bin_width_norm) {
    h_data->Scale(1., "width");
    h_sig->Scale(1., "width");
    h_loW->Scale(1., "width");
    h_midW->Scale(1., "width");
    h_hiW->Scale(1., "width");
    y_label = "Events / MeV";
  }

  if (ymax < 0) ymax = h_data->GetBinContent(h_data->GetMaximumBin()) * 1.6;
  if (ymax > 0) mnvPlotter.axis_maximum = ymax;

  // Prepare stack
  std::string legend_name =
      GetTruthClassification_LegendLabel(kWSideband_Signal);
  h_sig->SetTitle(legend_name.c_str());

  legend_name = GetTruthClassification_LegendLabel(kWSideband_Low);
  h_loW->SetTitle(legend_name.c_str());

  legend_name = GetTruthClassification_LegendLabel(kWSideband_Mid);
  h_midW->SetTitle(legend_name.c_str());

  legend_name = GetTruthClassification_LegendLabel(kWSideband_High);
  h_hiW->SetTitle(legend_name.c_str());

  SetHistColorScheme(h_sig, int(kWSideband_Signal),
                     sidebands::kWSideband_ColorScheme);
  SetHistColorScheme(h_loW, int(kWSideband_Low),
                     sidebands::kWSideband_ColorScheme);
  SetHistColorScheme(h_midW, int(kWSideband_Mid),
                     sidebands::kWSideband_ColorScheme);
  SetHistColorScheme(h_hiW, int(kWSideband_High),
                     sidebands::kWSideband_ColorScheme);

  TObjArray* array = new TObjArray();
  array->Add(h_sig);
  array->Add(h_loW);
  array->Add(h_midW);
  array->Add(h_hiW);

  // Draw
  mnvPlotter.DrawDataStackedMC(h_data, array, pot_scale, "TR", "Data", -1, -1,
                               1001,
                               Form("%s %s", variable->m_hists.m_xlabel.c_str(),
                                    variable->m_units.c_str()),
                               y_label.c_str());

  //  mnvPlotter.WritePreliminary("TL");
  mnvPlotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);

  std::ostringstream oss;
  oss << fit_str << "fit " << tag;
  std::string title = oss.str();

  mnvPlotter.AddHistoTitle(title.c_str());
  mnvPlotter.MultiPrint(&cE, label, "pdf");

  delete h_data;
  delete array;
}

//==============================================================================
// Backgrounds
//==============================================================================
/*
void PlotBG_ErrorGroup(Plotter p, std::string error_group_name,
                       bool do_tuned = false, double ignore_threshold = 0.,
                       double ymax = -1.) {
  TCanvas canvas ("c1","c1");

  // Make sure we remembered to load the source histos from the input file.
  assert(p.m_variable->m_hists.m_tuned_bg);
  assert(p.m_variable->m_hists.m_bg.hist);

  PlotUtils::MnvH1D* tuned_bg =
(PlotUtils::MnvH1D*)p.m_variable->m_hists.m_tuned_bg->Clone("tuned_bg");
  PlotUtils::MnvH1D* bg       =
(PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg.hist->Clone("bg");

  // X label
  p.SetXLabel(tuned_bg);
  p.SetXLabel(bg);

  // For groups MnvPlotter sets the ymax to headroom*axis_maximum_group
  // Or for summary, ymax = axis_maximum
  if (ymax != -1.) {
    p.m_mnv_plotter.axis_maximum = ymax;
    p.m_mnv_plotter.axis_maximum_group = ymax;
    p.m_mnv_plotter.headroom = 1.;
  }

  // Draw
  // XXX WARNING: potential problems when do_cov_area_norm = true
  if (do_tuned) {
    p.m_mnv_plotter.DrawErrorSummary(tuned_bg, "TR", p.m_include_stat, true,
                                     ignore_threshold, p.m_do_cov_area_norm,
                                     error_group_name, p.m_do_frac_unc);
  }
  else {
    p.m_mnv_plotter.DrawErrorSummary(bg, "TR", p.m_include_stat, true,
                                     ignore_threshold, p.m_do_cov_area_norm,
                                     error_group_name, p.m_do_frac_unc);
  }

  // Plot Title (For systematics, this has to be after drawing)
  if (do_tuned)
    p.SetTitle("Tuned Background");
  else
    p.SetTitle("Untuned Background");


  std::string tuned_str = do_tuned ? "BGTuned" : "BGUntuned";

  std::string outfile_name = Form("ErrorSummary_%s_%s_%s_%s_%s_%s",
                                  tuned_str.c_str(),
                                  p.m_variable->Name().c_str(),
                                  p.m_do_frac_unc_str.c_str(),
                                  p.m_do_cov_area_norm_str.c_str(),
                                  GetSignalFileTag(p.m_signal_definition).c_str(),
                                  error_group_name.c_str());

  p.m_mnv_plotter.MultiPrint(&canvas, outfile_name, "png");
}
*/

void PlotBG_ErrorSummary(Plotter p, bool do_tuned = false) {
  SetErrorGroups(p.m_mnv_plotter, false);
  PlotUtils::MnvH1D* bg;

  if (do_tuned)
    bg =
        (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_tuned_bg->Clone("tuned_bg");
  else
    bg = (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_bg.hist->Clone("bg");

  double detector_threshold = 0., detector_ymax = 0.2;
  double FSI_threshold = 0., FSI_ymax = 0.2;
  double Int_threshold = 0., Int_ymax = 0.15;
  std::string name = p.m_variable->Name();
  if (name == "enu") {
    Int_threshold = 0.0;
    Int_ymax = 0.12;
  } else if (name == "pmu") {
    Int_threshold = 0.0;
    Int_ymax = 0.12;
  } else if (name == "ptmu") {
    Int_threshold = 0.0;
    Int_ymax = 0.3;
  } else if (name == "pzmu") {
    Int_threshold = 0.0;
    Int_ymax = 0.1;
  } else if (name == "q2") {
    Int_threshold = 0.0;
    Int_ymax = 0.2;
  } else if (name == "thetamu_deg") {
    Int_threshold = 0.0;
    Int_ymax = 0.15;
  } else if (name == "thetapi_deg") {
    Int_threshold = 0.0;
    Int_ymax = 0.12;
  } else if (name == "tpi") {
    Int_threshold = 0.0;
    Int_ymax = 0.12;
  } else if (name == "wexp") {
    Int_threshold = 0.0;
    Int_ymax = 0.2;
  }

  std::string tuned_str = do_tuned ? "BGTuned" : "BGUntuned";

  // name, ignore threshold, ymax
  //  Plot_ErrorGroup(p, bg, "LEGENDONLY", tuned_str, 0.0);
  //  Plot_ErrorGroup(p, bg, "", tuned_str, 0.0);
  //  //  Plot_ErrorGroup(p, bg, "2p2h", tuned_str, 0.0, 0.05);
  //  Plot_ErrorGroup(p, bg, "Detector", tuned_str, detector_threshold,
  //                  detector_ymax);
  //  Plot_ErrorGroup(p, bg, "Flux", tuned_str, 0.0, 0.15);
  //  Plot_ErrorGroup(p, bg, "Genie_FSI_pions", tuned_str, FSI_threshold,
  //  FSI_ymax); Plot_ErrorGroup(p, bg, "Genie_FSI_nucleons", tuned_str,
  //  FSI_threshold,
  //                  FSI_ymax);
  //  Plot_ErrorGroup(p, bg, "Genie_InteractionModel", tuned_str, Int_threshold,
  //                  Int_ymax);
  //  Plot_ErrorGroup(p, bg, "NonResPi", tuned_str, 0.0, 0.2);
  //  //  Plot_ErrorGroup(p, bg, "RPA", tuned_str, 0.0, 0.05);
  //  //  Plot_ErrorGroup(p, bg, "Michel", tuned_str, 0.0, 0.05);
  //  //  Plot_ErrorGroup(p, bg, "GENIE", tuned_str, 0.0, 0.25);
  //  //  Plot_ErrorGroup(p, bg, "Target", tuned_str, 0.0, 0.15);
  //  Plot_ErrorGroup(p, bg, "Others", tuned_str, 0.0, 0.30);
  //  Plot_ErrorGroup(p, bg, "Cross Section Models", tuned_str, 0.0, 0.05);
  //  //  Plot_ErrorGroup(p, bg, "PhysicsModel", tuned_str, 0.0, 0.05);
  //  Plot_ErrorGroup(p, bg, "Muon", tuned_str, 0.0, 0.05);

  //  Plot_ErrorGroup(p, bg, "LEGENDONLY", tuned_str, 0.0, 0.1);
  Plot_ErrorGroup(p, bg, "", tuned_str, 0.0, 0.3);
  //  Plot_ErrorGroup(p, bg, "Detector", tuned_str, 0.0, 0.1);
  Plot_ErrorGroup(p, bg, "Pion Reconstruction", tuned_str, 0.0, 0.06);
  Plot_ErrorGroup(p, bg, "Flux", tuned_str, 0.0, 0.2);
  Plot_ErrorGroup(p, bg, "GENIE FSI", tuned_str, 0.0, 0.15);
  Plot_ErrorGroup(p, bg, "Muon", tuned_str, 0.0, 0.1);
  Plot_ErrorGroup(p, bg, "Others", tuned_str, 0.0, 0.05);

  AddingGroups(p.m_mnv_plotter, *bg, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *bg, "RESPi");
  AddingGroups(p.m_mnv_plotter, *bg, "CCQE");
  AddingGroups(p.m_mnv_plotter, *bg, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *bg, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *bg, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *bg, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  Plot_ErrorGroup(p, bg, "Cross Section Models", tuned_str, 0.0, 0.25);
  Plot_ErrorGroup(p, bg, "RESPi", tuned_str, 0.0, 0.25);
  Plot_ErrorGroup(p, bg, "NonRESPi", tuned_str, 0.0, 0.15);
  Plot_ErrorGroup(p, bg, "CCQE", tuned_str, 0.0, 0.1);
  Plot_ErrorGroup(p, bg, "Coherent-Diffractive", tuned_str, 0.0, 0.08);
  Plot_ErrorGroup(p, bg, "DIS-Hadronization", tuned_str, 0.0, 0.05);
  Plot_ErrorGroup(p, bg, "MnvTunes", tuned_str, 0.0, 0.02);
  Plot_ErrorGroup(p, bg, "Elastic", tuned_str, 0.0);
}

/*
void PlotBG_ErrorSummary(Plotter p, bool do_tuned = false) {
  SetErrorGroups(p.m_mnv_plotter);

  //name, ignore threshold, ymax
  PlotBG_ErrorGroup(p, "",                       do_tuned, 0.0,  0.6); // plot
all groups together PlotBG_ErrorGroup(p, "Flux",                   do_tuned,
0.0,  0.2); PlotBG_ErrorGroup(p, "Detector",               do_tuned, 0.0,  0.6);
  PlotBG_ErrorGroup(p, "Genie_FSI",              do_tuned, 0.02, 0.6);
  PlotBG_ErrorGroup(p, "Genie_InteractionModel", do_tuned, 0.02, 0.5);
  PlotBG_ErrorGroup(p, "NonResPi",               do_tuned, 0.0,  0.1);
  PlotBG_ErrorGroup(p, "2p2h",                   do_tuned, 0.0,  0.1);
  PlotBG_ErrorGroup(p, "RPA",                    do_tuned, 0.0,  0.1);
  //PlotBG_ErrorGroup(p, "LowQ2Pi");
}
*/

//==============================================================================
// Hack-y functions
//==============================================================================
void PlotTotalError(PlotUtils::MnvH1D* hist, std::string method_str,
                    Plotter p) {
  TCanvas cF("c4", "c4");
  TH1D* hTotalErr = (TH1D*)hist
                        ->GetTotalError(p.m_include_stat, p.m_do_frac_unc,
                                        p.m_do_cov_area_norm)
                        .Clone(Form("h_total_err_errSum_%d", __LINE__));
  hTotalErr->SetTitle(Form("Total Uncertainty (%s)", method_str.c_str()));
  if (p.m_do_frac_unc) hTotalErr->GetYaxis()->SetRangeUser(0, 1);
  // hTotalErr->Scale(0.408);
  hTotalErr->Draw();
  cF.Print(Form("TotalUncertainty_%s_%s_%s.png", p.m_do_frac_unc_str.c_str(),
                p.m_do_cov_area_norm_str.c_str(), method_str.c_str()));
}

void PlotStatError(PlotUtils::MnvH1D* hist, Plotter p, std::string tag,
                   double ymax = -1., std::string ylabel = "") {
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  p.SetXLabel(hist);
  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;
  // Y-axis label
  if (ylabel != "") hist->GetYaxis()->SetTitle(ylabel.c_str());

  TH1D* stat_error = (TH1D*)hist->GetStatError(p.m_do_frac_unc).Clone(uniq());

  stat_error->Draw();

  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "png");
}

void PlotVertBand(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* hist, Plotter p) {
  TCanvas cF("c4", "c4");
  TH1* h1 = (TH1*)hist->GetVertErrorBand(band.c_str())
                ->GetErrorBand(p.m_do_frac_unc, p.m_do_cov_area_norm)
                .Clone(Form("Pmu_%s_%s", band.c_str(), method_str.c_str()));
  h1->SetTitle(Form("%s Uncertainty (%s)", band.c_str(), method_str.c_str()));
  p.SetTitle();
  p.SetXLabel(hist);
  // h1->Scale(0.408);
  h1->Draw("h");
  cF.Print(Form("%s_band_%s.png", band.c_str(), method_str.c_str()));
}

void PlotVertBandAllUniverses(std::string band, std::string method_str,
                              PlotUtils::MnvH1D* hist, Plotter p) {
  TCanvas cF("c4", "c4");
  p.SetTitle();
  p.SetXLabel(hist);
  hist->GetVertErrorBand(band.c_str())->DrawAll("", true);
  cF.Print(
      Form("%s_band_all_universes_%s.png", band.c_str(), method_str.c_str()));
}

void PlotVertUniverse(std::string band, unsigned int universe,
                      std::string method_str, PlotUtils::MnvH1D* hist) {
  TCanvas cF("c1", "c1");
  TH1D* h1 = hist->GetVertErrorBand(band.c_str())->GetHist(universe);

  h1->SetLineColor(kBlack);
  h1->SetLineStyle(1);
  h1->Draw("hist");
  cF.Print(Form("Enu_%s_band_universe%i_%s.png", band.c_str(), universe + 1,
                method_str.c_str()));
}

void PlotDataMC(PlotUtils::MnvH1D* mc, PlotUtils::MnvH1D* data, Plotter p,
                std::string tag) {
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  p.m_mnv_plotter.DrawDataMC(data, mc, pot_scale, "TR");
  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "png");
}

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

void PlotTH1_1(TH1* h1, std::string tag, double ymax = -1,
               bool do_log_scale = false, bool do_fit = false) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  // TH1::SetDefaultSumw2();

  TCanvas cF("c4", "c4");

  h1->SetTitle(tag.c_str());
  h1->Draw("HIST");

  if (do_log_scale) cF.SetLogy();

  cF.Update();

  if (ymax > 0) h1->GetYaxis()->SetRangeUser(0, ymax);

  if (do_fit) {
    // gaussian fit
    h1->Fit("gaus", "0");
    h1->Draw("");
    cF.Update();
    TF1* fit1 = (TF1*)h1->GetFunction("gaus");
    fit1->SetLineColor(kRed);
    fit1->Draw("SAME");
  }

  /*
  // draw a line
  TLine *line = new TLine(0,0,0,cF.GetUymax());
  line->SetLineColor(kBlue);
  line->Draw();
  */

  cF.Print(Form("%s.png", tag.c_str()));
  // cF.Print(Form("%s.eps", tag.c_str()));
}

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

void PlotMC(PlotUtils::MnvH1D* doomyhist, Plotter p, std::string tag,
            double ymax = -1., std::string ylabel = "") {
  TCanvas canvas("c1", "c1");
  double pot_scale = p.m_data_pot / p.m_mc_pot;
  PlotUtils::MnvH1D* hist = nullptr;
  if (p.m_variable->Name() == "q2_true") {
    hist = RebinQ2Plot(*doomyhist);
    canvas.SetLogx();
  } else if (p.m_variable->Name() == "enu_true" ||
             p.m_variable->Name() == "pmu_true" ||
             p.m_variable->Name() == "ptmu_true" ||
             p.m_variable->Name() == "pzmu_true") {
    hist = RebinningtoGeV(*doomyhist, p.m_variable->Name());
  } else {
    hist = (PlotUtils::MnvH1D*)doomyhist->Clone("hist");
  }

  p.SetXLabel(hist);
  // Y-axis range
  if (ymax > 0) p.m_mnv_plotter.axis_maximum = ymax;
  // Y-axis label
  if (ylabel != "") hist->GetYaxis()->SetTitle(ylabel.c_str());
  // PlotUtils::MnvH1D* tmp_bg = nullptr;
  p.m_mnv_plotter.DrawMCWithErrorBand(
      hist);  // I think that this call only shows stat errors.
  p.m_mnv_plotter.MultiPrint(&canvas, tag.c_str(), "pdf");
}

void PlotRatio(PlotUtils::MnvH1D* num, PlotUtils::MnvH1D* denom, std::string v,
               double norm, std::string l, bool fixRange, std::string ylabel,
               std::string xlabel) {
  // char* vchar = &v[0];
  std::string label(Form("Ratio_%s", v.c_str()));
  // char* labchar = &label[0];
  const bool drawSysLines = false;
  const bool drawOneLine = true;
  double Min = -1., Max = -1.;

  if (fixRange) {
    Min = 0.7;
    Max = 1.3;
  }
  const double plotMin = Min;
  const double plotMax = Max;
  const bool covAreaNormalize = false;
  double titleSize = 0.05;

  cout << "Plotting ratio " << label << endl;

  TCanvas* c2 = new TCanvas();
  if (v == "Q2") c2->SetLogx();
  denom->GetXaxis()->SetTitle(xlabel.c_str());
  //  num->GetYaxis()->SetTitle(xlabel.c_str());
  PlotUtils::MnvPlotter* ratio = new PlotUtils::MnvPlotter();
  ratio->PlotUtils::MnvPlotter::DrawDataMCRatio(
      num, denom, norm, drawSysLines, drawOneLine, plotMin, plotMax,
      ylabel.c_str(), covAreaNormalize);
  ratio->AddHistoTitle(Form("%s %s", label.c_str(), l.c_str()), titleSize);
  c2->Print(Form("%s_%s.pdf", label.c_str(), l.c_str()));
}

void PlotRatio1(PlotUtils::MnvH1D* num, PlotUtils::MnvH1D* denom,
                const char* label, bool fixRange = true) {
  cout << "Plotting ratio " << label << endl;
  TCanvas* c = new TCanvas;
  PlotUtils::MnvH1D* ratio = (PlotUtils::MnvH1D*)num->Clone(uniq());
  ratio->Divide(num, denom);

  // ratio->GetXaxis()->SetRangeUser(0, 2);
  ratio->GetYaxis()->SetRangeUser(0, 1);

  // double mindiff=fabs(ratio->GetMinimum()-1);
  // double maxdiff=fabs(ratio->GetMaximum()-1);
  // double ydiff=std::max(mindiff, maxdiff);
  ratio->SetMarkerStyle(kFullCircle);
  ratio->Draw("P");

  if (fixRange) {
    double ydiff = 0.2;
    ratio->SetMinimum(1.0 - ydiff);
    ratio->SetMaximum(1.0 + ydiff);
  }
  c->Print(Form("%s.pdf", label));
}

void PlotRatioVec(std::vector<PlotUtils::MnvH1D*> num, PlotUtils::MnvH1D* denom,
                  std::string v, double norm, std::string l, bool fixRange,
                  std::string ylabel, std::string xlabel) {
  // char* vchar = &v[0];
  std::string label(Form("Ratio_%s", v.c_str()));
  // char* labchar = &label[0];
  const bool drawSysLines = false;
  const bool drawOneLine = true;
  double Min = -1., Max = -1.;
  auto legend = new TLegend(0.75, 0.8, 1., 0.9);
  if (fixRange) {
    Min = 0.7;
    Max = 2.;
  }
  const double plotMin = Min;
  const double plotMax = Max;
  const bool covAreaNormalize = false;
  double titleSize = 0.05;
  cout << "Plotting ratio " << label << endl;
  std::vector<TH1*> ratios;
  TH1* dummyratio;
  for (int i = 0; i < (int)num.size(); ++i) {
    dummyratio = (TH1*)num[i]->Clone(uniq());
    ratios.push_back(dummyratio);
  }

  TCanvas* c2 = new TCanvas();
  if (v == "q2" || v == "Q2") c2->SetLogx();

  denom->GetXaxis()->SetTitle(xlabel.c_str());
  PlotUtils::MnvPlotter* ratio = new PlotUtils::MnvPlotter();
  ratio->PlotUtils::MnvPlotter::DrawDataMCRatio(
      num[0], denom, norm, drawSysLines, drawOneLine, plotMin, plotMax,
      ylabel.c_str(), covAreaNormalize);
  ratio->AddHistoTitle(Form("%s %s", label.c_str(), l.c_str()), titleSize);

  for (int i = 0; i < (int)ratios.size(); ++i) {
    ratios[i]->Divide((TH1*)denom);
    ratios[i]->SetMarkerStyle(20);
    ratios[i]->SetMarkerSize(1.0);
    ratios[i]->SetLineWidth(3);
    ratios[i]->SetLineColor(30 + (i * 10));

    ratios[i]->Draw("SAME");
    legend->AddEntry(ratios[i], ratios[i]->GetTitle(), "lep");
  }
  legend->Draw();
  cout << "Plotting ratio 2" << endl;
  c2->Print(Form("%s_%s.pdf", label.c_str(), l.c_str()));
}

//==============================================================================
// Migration & Efficiency
//==============================================================================
// Helpers
TH2D* GetHistWithUnderOverFlow(TH2D* h) {
  // Create new binning with under/overflow
  UInt_t nx = h->GetXaxis()->GetNbins() + 2;
  Double_t* xbins = new Double_t[nx + 1];
  for (UInt_t i = 0; i < nx; i++)
    xbins[i] = h->GetXaxis()->GetBinLowEdge(i + 0);
  xbins[nx] = xbins[nx - 1] + h->GetXaxis()->GetBinWidth(nx);

  UInt_t ny = h->GetYaxis()->GetNbins() + 2;
  Double_t* ybins = new Double_t[ny + 1];
  for (UInt_t i = 0; i < ny; i++)
    ybins[i] = h->GetYaxis()->GetBinLowEdge(i + 0);
  ybins[ny] = ybins[ny - 1] + h->GetYaxis()->GetBinWidth(ny);

  // Create new histogram with under/overflow
  TH2D* htmp = new TH2D(h->GetName(), h->GetTitle(), nx, xbins, ny, ybins);
  htmp->Sumw2();

  // Fill the new histogram including the overflows
  for (UInt_t i = 1; i <= nx; i++) {
    double xbincenter = htmp->GetXaxis()->GetBinCenter(i);
    for (UInt_t j = 1; j <= ny; j++) {
      double ybincenter = htmp->GetYaxis()->GetBinCenter(j);

      int binnumber = htmp->FindBin(xbincenter, ybincenter);
      int oldbinnumber = h->FindBin(xbincenter, ybincenter);
      // std::cout << binnumber << ":(" << xbincenter << "," << ybincenter << ")
      // ";
      double bincontent = h->GetBinContent(oldbinnumber);
      double binerror = h->GetBinError(oldbinnumber);
      htmp->SetBinContent(binnumber, bincontent);
      htmp->SetBinError(binnumber, binerror);
    }
    // std::cout << "\n";
  }

  // Fill underflow specially
  double xval_firstbin =
      h->GetXaxis()->GetBinLowEdge(1);  // i.e. high edge of underflow
  double yval_firstbin = h->GetYaxis()->GetBinLowEdge(1);

  // Pretty weird way of finding new underflow bin number IMO.
  // I found it on the internet. Might depend on binning or on SetCanExtend?
  double xval_underflow = xval_firstbin - 1.;
  double yval_underflow = yval_firstbin - 1.;

  // New hist's bin number, corresponding to old hist's underflow bin
  int underflowbin = htmp->FindBin(xval_underflow, yval_underflow);

  // Get underflow content and set
  double underflowbincontent = h->GetBinContent(0);
  double underflowbinerror = h->GetBinError(0);
  htmp->SetBinContent(underflowbin, underflowbincontent);
  htmp->SetBinError(underflowbin, underflowbinerror);

  // Restore the number of entries (for when using weights!=1)
  htmp->SetEntries(h->GetEffectiveEntries());
  return htmp;
}

TH2D* RowNormalize(TH2D* h) {
  int first_bin = 0;
  int last_bin_X = h->GetXaxis()->GetNbins() + 1;
  int last_bin_Y = h->GetYaxis()->GetNbins() + 1;

  TH2D* tmp = (TH2D*)h->Clone();
  tmp->Reset();

  for (int y = first_bin; y <= last_bin_Y; ++y) {
    Double_t norm = 0.;
    for (int x = first_bin; x <= last_bin_X; ++x)
      norm += h->GetBinContent(x, y);
    if (fabs(norm) > 1E-8) {
      for (int x = first_bin; x <= last_bin_X; ++x) {
        double percentage = 100 * h->GetBinContent(x, y) / norm;
        tmp->SetBinContent(x, y, percentage);
      }
    }
  }
  return tmp;
}

void PlotMigration_AbsoluteBins(PlotUtils::MnvH2D* hist, std::string name) {
  TCanvas c("c1", "c1");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  mnv_plotter.SetRedHeatPalette();
  bool draw_as_matrix = true;
  gStyle->SetHistMinimumZero(kFALSE);
  mnv_plotter.DrawNormalizedMigrationHistogram(hist, draw_as_matrix);
  c.Update();
  c.Print(Form("Migration_AbsBins_%s.png", name.c_str()));
}

void PlotMigration_VariableBins(PlotUtils::MnvH2D* hist, std::string name) {
  TGaxis::SetExponentOffset(-0.035, -0.048, "x");
  // TH2D* htmp = GetHistWithUnderOverFlow(hist);
  TH2D* htmp = (TH2D*)hist->Clone();
  TH2D* htmp2 = RowNormalize(htmp);
  htmp2->GetXaxis()->SetTitle(Form("%s %s", name.c_str(), "Reco"));
  htmp2->GetYaxis()->SetTitle(Form("%s %s", name.c_str(), "True"));
  TCanvas c("c1", "c1");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  mnv_plotter.SetRedHeatPalette();
  bool draw_as_matrix = false;
  gStyle->SetHistMinimumZero(kFALSE);
  mnv_plotter.DrawNormalizedMigrationHistogram(htmp2, draw_as_matrix);
  // gStyle->SetPaintTextFormat("2.0f");
  // htmp2->SetMarkerSize(2);
  // htmp2->Draw("colz text");
  c.Update();
  c.Print(Form("Migration_VarBins_%s.png", name.c_str()));
  // c.SetLogz();
  // c.Update();
  // c.Print("WMigrationMatrix_Wbins_logz.png");
  TGaxis::SetExponentOffset(0, 0, "x");
}

void PlotEfficiency_ErrorSummary(Plotter p) {
  SetErrorGroups(p.m_mnv_plotter, false);
  PlotUtils::MnvH1D* eff =
      (PlotUtils::MnvH1D*)p.m_variable->m_hists.m_efficiency->Clone(uniq());
  //  Plot_ErrorGroup(p, eff, "LEGENDONLY", "Eff", 0.0, 0.3);
  //  Plot_ErrorGroup(p, eff, "", "Eff", 0.0, 0.1);
  //  Plot_ErrorGroup(p, eff, "Flux", "Eff", 0.0, 0.01);
  //  Plot_ErrorGroup(p, eff, "Detector", "Eff", 0.0, 0.01);
  //  Plot_ErrorGroup(p, eff, "Genie_FSI_pions", "Eff", 0.0, 0.1);
  //  Plot_ErrorGroup(p, eff, "Genie_FSI_nucleons", "Eff", 0.0, 0.1);
  //  Plot_ErrorGroup(p, eff, "Genie_InteractionModel", "Eff", 0.0, 0.2);
  //  Plot_ErrorGroup(p, eff, "NonResPi", "Eff", 0.0, 0.05);
  //  //  Plot_ErrorGroup(p, eff, "2p2h", "Eff", 0.0, 0.001);
  //  //  Plot_ErrorGroup(p, eff, "RPA", "Eff", 0.0, 0.01);
  //  // Plot_ErrorGroup(p, eff, "Michel", "Eff", 0.0, 0.15);
  //  //  Plot_ErrorGroup(p, eff, "GENIE", "Eff", 0.0, 0.15);
  //  //  Plot_ErrorGroup(p, eff, "Target", "Eff", 0.0, 0.15);
  //  Plot_ErrorGroup(p, eff, "Others", "Eff", 0.0, 0.03);
  //  Plot_ErrorGroup(p, eff, "Cross Section Models", "Eff", 0.0, 0.04);
  //  //  Plot_ErrorGroup(p, eff, "PhysicsModel", "Eff", 0.0, 0.15);
  //  Plot_ErrorGroup(p, eff, "Muon", "Eff", 0.0, 0.05);

  double max = 0.15;
  if (p.m_variable->Name() == "ptmu_true") max = 0.15;
  //  Plot_ErrorGroup(p, eff, "LEGENDONLY", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "", "Eff", 0.0, max);
  //  Plot_ErrorGroup(p, eff, "Detector", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "Pion Reconstruction", "Eff", 0.0, 0.04);
  Plot_ErrorGroup(p, eff, "Flux", "Eff", 0.0, 0.01);
  Plot_ErrorGroup(p, eff, "GENIE FSI", "Eff", 0.0, 0.02);
  Plot_ErrorGroup(p, eff, "Muon", "Eff", 0.0, 0.08);
  Plot_ErrorGroup(p, eff, "Others", "Eff", 0.0, 0.03);

  AddingGroups(p.m_mnv_plotter, *eff, "NonRESPi");
  AddingGroups(p.m_mnv_plotter, *eff, "RESPi");
  AddingGroups(p.m_mnv_plotter, *eff, "CCQE");
  AddingGroups(p.m_mnv_plotter, *eff, "Coherent-Diffractive");
  AddingGroups(p.m_mnv_plotter, *eff, "DIS-Hadronization");
  AddingGroups(p.m_mnv_plotter, *eff, "MnvTunes");
  AddingGroups(p.m_mnv_plotter, *eff, "Elastic");

  SetErrorGroups(p.m_mnv_plotter, true);
  Plot_ErrorGroup(p, eff, "Cross Section Models", "Eff", 0.0, 0.05);
  Plot_ErrorGroup(p, eff, "RESPi", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "NonRESPi", "Eff", 0.0, 0.02);
  Plot_ErrorGroup(p, eff, "CCQE", "Eff", 0.0, 0.01);
  Plot_ErrorGroup(p, eff, "Coherent-Diffractive", "Eff", 0.0, 0.02);
  Plot_ErrorGroup(p, eff, "DIS-Hadronization", "Eff", 0.0, 0.0002);
  Plot_ErrorGroup(p, eff, "MnvTunes", "Eff", 0.0, 0.1);
  Plot_ErrorGroup(p, eff, "Elastic", "Eff", 0.0);
}

void PlotStackedXSec(Variable* variable, const PlotUtils::MnvH1D* h_data,
                     const TObjArray& array_mc, float data_pot, float mc_pot,
                     SignalDefinition signal_definition, std::string tag = "",
                     double ymax = -1, bool draw_arrow = false,
                     bool do_bin_width_norm = true) {
  // Never don't clone when plotting
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)h_data->Clone("data");
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  mnvPlotter.axis_minimum = 0.001;
  if (ymax > 0) mnvPlotter.axis_maximum = ymax;
  if (tag == "FSP") mnvPlotter.legend_offset_x = .15;
  if (tag == "Hadrons") mnvPlotter.legend_offset_x = .06;
  if (tag == "Wexp") mnvPlotter.legend_offset_x = .06;
  mnvPlotter.legend_text_size = 0.0405;

  double pot_scale = 1;
  std::string label = Form("Breakdown_CrossSection_%s_%s_%s",
                           GetSignalFileTag(signal_definition).c_str(),
                           variable->m_label.c_str(), tag.c_str());

  std::string y_label = "Events";

  // bin width norm
  if (do_bin_width_norm) {
    data->Scale(1.e42, "width");
    for (auto h : array)
      dynamic_cast<PlotUtils::MnvH1D*>(h)->Scale(1.e42, "width");
    y_label = "d#sigma/d" + variable->m_hists.m_xlabel + " (10^{-42} cm^{2}/" +
              variable->m_units + "/nucleon)";
  }

  TCanvas cE("c1", "c1");
  //  cE.SetLogy();
  mnvPlotter.DrawDataStackedMC(data, &array, pot_scale, "TR", "Data", -1, -1,
                               1001, variable->m_hists.m_xlabel.c_str(),
                               y_label.c_str());
  if (draw_arrow) {
    // double arrow_height = data->GetBinContent(data->GetMaximumBin()) *
    //                      data->GetNormBinWidth()/data->GetBinWidth(data->GetMaximumBin());
    double arrow_height = 250;
    double arrow_location = signal_definition.m_w_max;
    mnvPlotter.AddCutArrow(arrow_location, 0.0, arrow_height, 200., "L");
  }

  // mnvPlotter.WritePreliminary("TL");
  mnvPlotter.AddPOTNormBox(data_pot, mc_pot, 0.3, 0.85);
  mnvPlotter.AddHistoTitle("Cross Section");
  //  cE.SetLogy();
  mnvPlotter.MultiPrint(&cE, label, "png");
  // mnvPlotter.MultiPrint(&cE, label, "eps");
}
#endif  // plotting_functions_h
