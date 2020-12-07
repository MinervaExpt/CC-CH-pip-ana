#ifndef runEffPurTable_C
#define runEffPurTable_C

#include "includes/CCPiMacroUtil.h"
#include "includes/CVUniverse.h"
#include "includes/common_stuff.h" // typedefs EventCount
#include "includes/CCPiEvent.h"
#include "includes/Cuts.h"
#include "../event_selection/EventSelectionTable.h"


//==============================================================================
// Loop and fill
//==============================================================================
void FillCounters(const CCPiMacroUtil& util, CVUniverse* universe,
                  const EDataMCTruth& type,
                  std::pair<EventCount*, EventCount*>& counters) {

  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  for(Long64_t i_event=0; i_event < n_entries; ++i_event){
    if (i_event%500000==0) std::cout << (i_event/1000) << "k " << std::endl;
    universe->SetEntry(i_event);
    CCPiEvent event(is_mc, is_truth, util.m_signal_definition, universe);
    ccpi_event::FillCounters(event, counters); // Does a lot of work
  } // events
  std::cout << "*** Done ***\n\n";
}


//==============================================================================
// Main
//==============================================================================
void runEffPurTable(int signal_definition_int = 0, const char* plist = "ALL") {
  // INIT MACRO UTILITY OBJECT
  const std::string macro("runEffPurTable");
  bool do_data = true, do_mc = true, do_truth = true;
  bool do_systematics = false, do_grid = false;
  CCPiMacroUtil util(signal_definition_int, plist, do_data, do_mc, do_truth,
                     do_systematics, do_grid);
  util.PrintMacroConfiguration(macro);

  // EFFICIENCY/PURITY COUNTERS
  // typdef EventCount map<ECut, double>
  EventCount n_remaining_sig, n_remaining_bg, n_remaining_data;
  std::pair<EventCount*, EventCount*> signal_bg_counters(&n_remaining_sig,
                                                         &n_remaining_bg);
  std::pair<EventCount*, EventCount*> data_count(&n_remaining_data, NULL);

  FillCounters(util, util.m_data_universe,  kData,  data_count);
  FillCounters(util, util.m_error_bands.at("cv").at(0), kMC,
               signal_bg_counters);
  FillCounters(util, util.m_error_bands_truth.at("cv").at(0), kTruth,
               signal_bg_counters);

  PrintEffPurTable(n_remaining_sig, n_remaining_bg, n_remaining_data,
                   util.m_data_pot, util.m_mc_pot);
}

#endif