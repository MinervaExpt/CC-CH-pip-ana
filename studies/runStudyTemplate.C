//==============================================================================
// Template for a generic loop data/mc and plot script.
// Assume that NO systematics are being analyzed.
// Good for stacked histograms, branch validation, residuals, etc.
//==============================================================================
#include <iostream>
#include <vector>

//#include "ccpion_common.h"  // GetPlaylistFile
#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/HadronVariable.h"
#include "includes/MacroUtil.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "plotting_functions.h"

// Forward declare my variables because we're hiding the header.
class Variable;
class HadronVariable;

namespace run_study_template {
//==============================================================================
// Do some event processing (e.g. make cuts, get best pion) and fill hists
//==============================================================================
void FillVars(CCPiEvent& event, const std::vector<Variable*>& variables) {
  CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;

  if (universe->ShortName() != "cv") return;

  // Process Event

  // Check cuts, check is_sideband, and get pion candidates
  PassesCutsInfo cv_cuts_info = PassesCuts(event);

  // Set all of this info to the event and/or universe
  std::tie(event.m_passes_cuts, event.m_is_w_sideband,
           event.m_passes_all_cuts_except_w, event.m_reco_pion_candidate_idxs) =
      cv_cuts_info.GetAll();
  event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
  universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);

  // Need to re-call this because the node cut efficiency systematic
  // needs a pion candidate to calculate its weight.
  if (is_mc) event.m_weight = universe->GetWeight();

  if (event.m_passes_cuts) ccpi_event::FillStackedHists(event, variables);
}

//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;
  HVar* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg",
               CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);
  std::vector<Var*> variables = {thetapi_deg, pmu};
  return variables;
}
}  // namespace run_study_template

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFill(const CCPi::MacroUtil& util, CVUniverse* universe,
                 const EDataMCTruth& type, std::vector<Variable*>& variables) {
  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);

    // For mc, get weight, check signal, and sideband
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    // WRITE THE FILL FUNCTION
    run_study_template::FillVars(event, variables);
  }  // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runStudyTemplate(std::string plist = "ME1L") {
  bool is_mc = true;
  const bool use_xrootd = true;
  const bool do_test_playlist = true;
  std::string mc_file_list =
      CCPi::GetPlaylistFile(plist, is_mc, do_test_playlist, use_xrootd);
  is_mc = false;
  std::string data_file_list =
      CCPi::GetPlaylistFile(plist, is_mc, do_test_playlist, use_xrootd);
  const int signal_definition_int = SignalDefinition::OnePiTracked().m_id;
  const bool is_grid = false;
  const bool do_truth = false;
  const bool do_systematics = false;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.m_name = "runStudyTemplate";
  util.PrintMacroConfiguration();

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<Variable*> variables = run_study_template::GetVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  //=========================================
  // Loop and Fill
  //=========================================
  LoopAndFill(util, util.m_data_universe, kData, variables);
  LoopAndFill(util, util.m_error_bands.at("cv").at(0), kMC, variables);

  for (auto v : variables) {
    std::string tag = v->Name();
    double ymax = -1;
    bool do_bwn = true;
    std::cout << "Plotting" << std::endl;
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kS),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
               v->Name(), "SSB", ymax, do_bwn);
  }

  std::cout << "Success" << std::endl;
}
