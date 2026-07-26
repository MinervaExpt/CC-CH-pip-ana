#ifndef Plotter_H
#define Plotter_H

#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif
#include "Constants.h"         // enum SignalDefinition
#include "SignalDefinition.h"  // GetSignalName
#include "TPad.h"
#include "Variable.h"

class Plotter {
 public:
  // Constructor with Var
  Plotter(Variable* variable, float mc_pot, float data_pot, bool do_frac_unc,
          bool do_cov_area_norm, bool include_stat,
          SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable(variable),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Constructor without var
  Plotter(float mc_pot, float data_pot, bool do_frac_unc, bool do_cov_area_norm,
          bool include_stat, SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable(nullptr),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Members
  MnvPlotter m_mnv_plotter;
  Variable* m_variable;
  float m_mc_pot;
  float m_data_pot;
  bool m_do_frac_unc;
  bool m_do_cov_area_norm;
  bool m_include_stat;
  SignalDefinition m_signal_definition;
  std::string m_do_frac_unc_str;
  std::string m_do_cov_area_norm_str;

  // Add X label
  void SetXLabel(PlotUtils::MnvH1D* hist) {
    std::string label;
    label = m_variable->m_hists.m_xlabel + " (" + m_variable->m_units + ")";
    //     m_variable->Name() == "q2"
    //         ? m_variable->m_hists.m_xlabel + "(GeV^{2})"
    //         : m_variable->m_hists.m_xlabel + " (" + m_variable->m_units +
    //         ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add X label
  void SetXLabel(TH1* hist) {
    std::string label =
        m_variable->m_hists.m_xlabel + " (" + m_variable->m_units + ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add title to the MnvPlotter
  void SetTitle() {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    std::string title = GetSignalName(m_signal_definition);
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }

  void SetTitle(std::string title) {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }
};

#endif  // Plotter_H
