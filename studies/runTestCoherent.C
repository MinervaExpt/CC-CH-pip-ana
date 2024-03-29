//==============================================================================
// Does the new |t| variable work?
//==============================================================================
#include <iostream>
#include <vector>

#include "ccpion_common.h"  // GetPlaylistFile
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

namespace run_test_coherent {
//==============================================================================
// Do some event processing (e.g. make cuts, get best pion) and fill hists
//==============================================================================
void FillVars(CCPiEvent& event, const std::vector<Variable*>& variables) {
  const CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;

  if (universe->ShortName() != "cv") return;

  event.m_passes_cuts = PassesCuts(event, event.m_is_w_sideband);
  event.m_highest_energy_pion_idx = GetHighestEnergyPionCandidateIndex(event);
  if (event.m_passes_cuts) ccpi_event::FillStackedHists(event, variables);
}

//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  typedef Variable Var;
  typedef HadronVariable HVar;
  HVar* t = new HVar("t", "#t", "mev", 30, 0, 4.e6, &CVUniverse::Gett);
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);
  std::vector<Var*> variables = {t, pmu};
  return variables;
}
}  // namespace run_test_coherent

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
    run_test_coherent::FillVars(event, variables);
  }  // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runTestCoherent(std::string plist = "ME1L") {
  //=========================================
  // Input tuples
  //=========================================
  bool is_mc = true;
  std::string mc_file_list, data_file_list;
  mc_file_list = GetPlaylistFile(plist, is_mc);
  is_mc = false;
  data_file_list = GetPlaylistFile(plist, is_mc);

  //=========================================
  // Init macro utility
  //=========================================
  const int signal_definition_int = 0;
  const std::string macro("runTestCoherent");
  const bool is_grid = false;
  const bool do_truth = false;
  const bool do_systematics = false;

  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<Variable*> variables = run_test_coherent::GetVariables();
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
    PlotCutVar(v, v->m_hists.m_selection_data, v->GetStackArray(kCCQE),
               util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
               v->Name(), "Channel", ymax, do_bwn);
  }

  std::cout << "Success" << std::endl;
}
