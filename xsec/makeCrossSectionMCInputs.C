#ifndef makeXsecMCInputs_C
#define makeXsecMCInputs_C

#include <cassert>
#include <ctime>
#include <functional>

#include "ccpion_common.h"
#include "includes/Binning.h"
#include "includes/CCPiEvent.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"
#include "includes/Cuts.h"
#include "includes/HadronVariable.h"
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/Variable.h"
#include "includes/common_functions.h"  // GetVar, WritePOT

//==============================================================================
// Helper Functions
//==============================================================================
namespace make_xsec_mc_inputs {
typedef Variable Var;
typedef HadronVariable HVar;

std::vector<Variable*> GetOnePiVariables(bool include_truth_vars = true) {
  const int nadphibins = 16;
  const double adphimin = -CCNuPionIncConsts::PI;
  const double adphimax = CCNuPionIncConsts::PI;

  HVar* tpi = new HVar("tpi", "T_{#pi}", "MeV", CCPi::GetBinning("tpi"),
                       &CVUniverse::GetTpi);

  HVar* tpi_mbr = new HVar("tpi_mbr", "T_{#pi} (MBR)", tpi->m_units,
                           CCPi::GetBinning("tpi"), &CVUniverse::GetTpiMBR);

  HVar* tpi_res = new HVar("tpi_res", "T_{#pi} (Residual)", tpi->m_units, 100,
                           -10., 10., &CVUniverse::GetTpiFResidual);

  HVar* thetapi_deg =
      new HVar("thetapi_deg", "#theta_{#pi}", "deg",
               CCPi::GetBinning("thetapi_deg"), &CVUniverse::GetThetapiDeg);

  Var* pmu = new Var("pmu", "p_{#mu}", "GeV/c", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);

  Var* thetamu_deg =
      new Var("thetamu_deg", "#theta_{#mu}", "deg",
              CCPi::GetBinning("thetamu_deg"), &CVUniverse::GetThetamuDeg);

  Var* enu = new Var("enu", "E_{#nu}", "GeV", CCPi::GetBinning("enu"),
                     &CVUniverse::GetEnu);

  Var* q2 = new Var("q2", "Q^{2}", "GeV^{2}/c^{2}", CCPi::GetBinning("q2"),
                    &CVUniverse::GetQ2);

  Var* q2_Aaron = new Var("q2_Aaron", "Q^{2}", "MeV^{2}",
                          CCPi::GetBinning("q2"), &CVUniverse::GetQ2);

  Var* q2_NoAaron = new Var("q2_NoAaron", "Q^{2}", "MeV^{2}",
                            CCPi::GetBinning("q2"), &CVUniverse::GetQ2);

  Var* wexp = new Var("wexp", "W_{exp}", "GeV/c^{2}", CCPi::GetBinning("wexp"),
                      &CVUniverse::GetWexp);

  Var* wexp_fit =
      new Var(sidebands::kFitVarString, wexp->m_hists.m_xlabel, wexp->m_units,
              CCPi::GetBinning("wexp_fit"), &CVUniverse::GetWexp);

  Var* ptmu = new Var("ptmu", "p^{T}_{#mu}", "GeV/c", CCPi::GetBinning("ptmu"),
                      &CVUniverse::GetPTmu);

  Var* pzmu = new Var("pzmu", "p^{||}_{#mu}", "GeV/c", CCPi::GetBinning("pzmu"),
                      &CVUniverse::GetPZmu);

  Var* mehreen_tpi =
      new Var("mtpi", "Mehreen T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
              &CVUniverse::GetTpiTrackless);

  Var* mehreen_thetapi_deg = new Var("mthetapi_deg", "Mehreen #theta_{#pi}",
                                     "deg", CCPi::GetBinning("thetapi_deg"),
                                     &CVUniverse::GetThetapitracklessDeg);

  HVar* mixthetapi_deg = new HVar("mixthetapi_deg", "#theta_{#pi}", "deg",
                                  CCPi::GetBinning("thetapi_deg"),
                                  &CVUniverse::GetMixedThetapiDeg);

  HVar* mixtpi = new HVar("mixtpi", "T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
                          &CVUniverse::GetMixedTpi);

  HVar* mixtpi_Aaron =
      new HVar("mixtpi_Aaron", "T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
               &CVUniverse::GetMixedTpi);

  HVar* mixtpi_NoAaron =
      new HVar("mixtpi_NoAaron", "T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
               &CVUniverse::GetMixedTpi);

  HVar* bkdtrackedtpi =
      new HVar("bkdtrackedtpi", "T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
               &CVUniverse::GetMixedTpi);

  HVar* bkdtracklesstpi =
      new HVar("bkdtracklesstpi", "T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
               &CVUniverse::GetMixedTpi);

  HVar* bkdmixtpi =
      new HVar("bkdmixtpi", "T_{#pi}", "MeV", CCPi::GetBinning("mtpi"),
               &CVUniverse::GetMixedTpi);
  // True Variables
  bool is_true = true;
  HVar* tpi_true =
      new HVar("tpi_true", "T_{#pi} True", tpi->m_units,
               tpi->m_hists.m_bins_array, &CVUniverse::GetTpiTrue, is_true);

  HVar* thetapi_deg_true =
      new HVar("thetapi_deg_true", "#theta_{#pi} True", thetapi_deg->m_units,
               thetapi_deg->m_hists.m_bins_array,
               &CVUniverse::GetThetapiTrueDeg, is_true);

  Var* pmu_true =
      new Var("pmu_true", "p_{#mu} True", pmu->m_units,
              pmu->m_hists.m_bins_array, &CVUniverse::GetPmuTrue, is_true);

  Var* thetamu_deg_true =
      new Var("thetamu_deg_true", "#theta_{#mu} True", thetamu_deg->m_units,
              thetamu_deg->m_hists.m_bins_array, &CVUniverse::GetThetamuTrueDeg,
              is_true);

  Var* enu_true =
      new Var("enu_true", "E_{#nu} True", enu->m_units,
              enu->m_hists.m_bins_array, &CVUniverse::GetEnuTrue, is_true);

  Var* q2_true =
      new Var("q2_true", "Q^{2} True", q2->m_units, q2->m_hists.m_bins_array,
              &CVUniverse::GetQ2True, is_true);
  Var* q2_Aaron_true =
      new Var("q2_Aaron_true", "Q^{2} True", q2->m_units,
              q2->m_hists.m_bins_array, &CVUniverse::GetQ2True, is_true);
  Var* q2_NoAaron_true =
      new Var("q2_NoAaron_true", "Q^{2} True", q2->m_units,
              q2->m_hists.m_bins_array, &CVUniverse::GetQ2True, is_true);

  Var* wexp_true =
      new Var("wexp_true", "W_{exp} True", wexp->m_units,
              wexp->m_hists.m_bins_array, &CVUniverse::GetWexpTrue, is_true);

  Var* ptmu_true =
      new Var("ptmu_true", "p^{T}_{#mu} True", ptmu->m_units,
              ptmu->m_hists.m_bins_array, &CVUniverse::GetPTmuTrue, is_true);

  Var* pzmu_true =
      new Var("pzmu_true", "p^{||}_{#mu} True", pzmu->m_units,
              pzmu->m_hists.m_bins_array, &CVUniverse::GetPZmuTrue, is_true);

  Var* mehreen_tpi_true = new Var("mtpi_true", "Mehreen T_{#pi} True", "MeV",
                                  mehreen_tpi->m_hists.m_bins_array,
                                  &CVUniverse::GetTrueTpi, is_true);

  Var* mehreen_thetapi_deg_true =
      new Var("mthetapi_deg_true", "Mehreen #theta_{#pi} True", "deg",
              mehreen_thetapi_deg->m_hists.m_bins_array,
              &CVUniverse::GetThetapitracklessTrueDeg, is_true);

  HVar* mixthetapi_deg_true =
      new HVar("mixthetapi_deg_true", "#theta_{#pi} True",
               mixthetapi_deg->m_units, mixthetapi_deg->m_hists.m_bins_array,
               &CVUniverse::GetMixedThetapiTrueDeg, is_true);

  HVar* mixtpi_true = new HVar("mixtpi_true", "T_{#pi} True", mixtpi->m_units,
                               mixtpi->m_hists.m_bins_array,
                               &CVUniverse::GetMixedTpiTrue, is_true);

  HVar* mixtpi_Aaron_true = new HVar(
      "mixtpi_Aaron_true", "T_{#pi} True", mixtpi->m_units,
      mixtpi->m_hists.m_bins_array, &CVUniverse::GetMixedTpiTrue, is_true);

  HVar* mixtpi_NoAaron_true = new HVar(
      "mixtpi_NoAaron_true", "T_{#pi} True", mixtpi->m_units,
      mixtpi->m_hists.m_bins_array, &CVUniverse::GetMixedTpiTrue, is_true);
  HVar* bkdtrackedtpi_true = new HVar(
      "bkdtrackedtpi_true", "T_{#pi} True", mixtpi->m_units,
      mixtpi->m_hists.m_bins_array, &CVUniverse::GetMixedTpiTrue, is_true);

  HVar* bkdtracklesstpi_true = new HVar(
      "bkdtracklesstpi_true", "T_{#pi} True", mixtpi->m_units,
      mixtpi->m_hists.m_bins_array, &CVUniverse::GetMixedTpiTrue, is_true);

  HVar* bkdmixtpi_true = new HVar("bkdmixtpi_true", "T_{#pi} True",
                                  mixtpi->m_units, mixtpi->m_hists.m_bins_array,
                                  &CVUniverse::GetMixedTpiTrue, is_true);
  // Ehad variables
  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                      &CVUniverse::GetEhad);
  Var* ehad_true =
      new Var("ehad_true", "ehad True", "MeV", ehad->m_hists.m_bins_array,
              &CVUniverse::GetEhadTrue);
  ehad_true->m_is_true = true;

  std::vector<Var*> variables = {// tpi, tpi_mbr,
                                 /* thetapi_deg,*/ pmu, thetamu_deg, enu, q2,
                                 //      tpi_res,
                                 wexp, wexp_fit, ptmu, pzmu, ehad, mixtpi,
                                 mixthetapi_deg};
  if (include_truth_vars) {
    //    variables.push_back(tpi_true);
    //    variables.push_back(thetapi_deg_true);
    variables.push_back(pmu_true);
    variables.push_back(thetamu_deg_true);
    variables.push_back(enu_true);
    variables.push_back(q2_true);
    //    variables.push_back(q2_Aaron_true);
    //    variables.push_back(q2_NoAaron_true);
    variables.push_back(wexp_true);
    variables.push_back(ptmu_true);
    variables.push_back(pzmu_true);
    variables.push_back(ehad_true);
    //    variables.push_back(mehreen_tpi_true);
    variables.push_back(mixtpi_true);
    //    variables.push_back(mixtpi_Aaron_true);
    //    variables.push_back(mixtpi_NoAaron_true);
    //        variables.push_back(bkdtrackedtpi_true);
    //        variables.push_back(bkdtracklesstpi_true);
    //    variables.push_back(bkdmixtpi_true);
    //    variables.push_back(mehreen_thetapi_deg_true);
    variables.push_back(mixthetapi_deg_true);
  }

  return variables;
}

std::map<std::string, Variable*> GetOnePiVariables_Map(
    bool include_truth_vars = true) {
  std::map<std::string, Var*> var_map;
  std::vector<Var*> var_vec = GetOnePiVariables(include_truth_vars);
  for (auto v : var_vec) var_map[v->Name()] = v;
  return var_map;
}

}  // namespace make_xsec_mc_inputs

std::vector<Variable*> GetAnalysisVariables(
    const SignalDefinition& signal_definition,
    const bool include_truth_vars = false) {
  using GetVariablesFn = std::function<std::vector<Variable*>(bool)>;
  std::map<int, GetVariablesFn> get_variables{
      {SignalDefinition::OnePi().m_id, make_xsec_mc_inputs::GetOnePiVariables},
      {SignalDefinition::OnePiTracked().m_id,
       make_xsec_mc_inputs::GetOnePiVariables},
      {SignalDefinition::Nuke().m_id, make_xsec_mc_inputs::GetOnePiVariables},
      {SignalDefinition::OnePiThetaPi().m_id,
       make_xsec_mc_inputs::GetOnePiVariables},
      {SignalDefinition::OnePiUntracked().m_id,
       make_xsec_mc_inputs::GetOnePiVariables}};
  return get_variables.at(signal_definition.m_id)(include_truth_vars);
}

void SavingStacked(TFile& fout, TObjArray plotsArray, std::string var,
                   std::string type) {
  int size = plotsArray.GetEntries();
  for (int i = 0; i < size; ++i) {
    fout.cd();
    TObject* obj =
        plotsArray.At(i)->Clone(Form("%s_%s_%d", var.c_str(), type.c_str(), i));
    PlotUtils::MnvH1D* h = dynamic_cast<PlotUtils::MnvH1D*>(obj);
    h->Write();
    fout.Flush();
  }
}

void SyncAllHists(Variable& v) {
  v.m_hists.m_selection_mc.SyncCVHistos();
  v.m_hists.m_selection_mc_tracked.SyncCVHistos();
  v.m_hists.m_selection_mc_untracked.SyncCVHistos();
  v.m_hists.m_selection_mc_mixed.SyncCVHistos();
  v.m_hists.m_selection_mc_no_tpi_weight.SyncCVHistos();
  v.m_hists.m_selection_mc_tracked_no_tpi_weight.SyncCVHistos();
  v.m_hists.m_selection_mc_untracked_no_tpi_weight.SyncCVHistos();
  v.m_hists.m_selection_mc_mixed_no_tpi_weight.SyncCVHistos();
  v.m_hists.m_bg.SyncCVHistos();
  v.m_hists.m_bg_loW.SyncCVHistos();
  v.m_hists.m_bg_midW.SyncCVHistos();
  v.m_hists.m_bg_hiW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_sig.SyncCVHistos();
  v.m_hists.m_wsidebandfit_loW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_midW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_hiW.SyncCVHistos();
  v.m_hists.m_effnum.SyncCVHistos();
  v.m_hists.m_effden.SyncCVHistos();
}

// Given a macro, make an output filename with a timestamp
std::string GetOutFilename(const CCPi::MacroUtil& util, const int run = 0) {
  auto time = std::time(nullptr);
  char tchar[100];
  std::strftime(tchar, sizeof(tchar), "%F", std::gmtime(&time));  // YYYY-MM-dd
  const std::string timestamp = tchar;
  const std::string outfile_format = "%s_%d%d%d%d_%s_%d_%s.root";
  return std::string(Form(outfile_format.c_str(), util.m_name.c_str(),
                          util.m_signal_definition.m_id,
                          int(util.m_do_systematics), int(util.m_do_truth),
                          int(util.m_is_grid), util.m_plist_string.c_str(), run,
                          timestamp.c_str()));
}

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillMCXSecInputs(const UniverseMap& error_bands,
                             const Long64_t n_entries, const bool is_truth,
                             const SignalDefinition& signal_definition,
                             std::vector<Variable*>& variables) {
  const bool is_mc = true;
  const bool onlytracked = signal_definition.m_do_tracked_michel_reco &&
                           !signal_definition.m_do_untracked_michel_reco;
  const bool onlyuntracked = !signal_definition.m_do_tracked_michel_reco &&
                             signal_definition.m_do_untracked_michel_reco;
  int selcount = 0;
  if (onlyuntracked && onlytracked) {
    std::cout << "Invalid configuration\n";
    std::exit(1);
  }

  for (auto band : error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes) universe->SetTruth(is_truth);
  }

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % (n_entries / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    // if (i_event == 20000) break;
    //      if(i_event%1000==0) std::cout << i_event << " / " << n_entries <<
    //      "\r"
    //      << std::flush;
    //    Variables that hold info about whether the CVU passes cuts
    PassesCutsInfo cv_cuts_info;
    bool checked_cv = false;
    assert(!error_bands.at("cv").empty() &&
           "\"cv\" error band is empty!  Can't set Model weight.");
    auto& cvUniv = error_bands.at("cv").at(0);
    cvUniv->SetEntry(i_event);
    if (is_truth) {
      for (auto error_band : error_bands) {  // Loop for truth
        std::vector<CVUniverse*> universes = error_band.second;
        for (auto universe : universes) {
          universe->SetEntry(i_event);
          CCPiEvent event(is_mc, is_truth, signal_definition, universe);
          universe->SetIsSignal(event.m_is_signal);
          event.m_weight = universe->GetWeight();
          universe->SetPassesTrakedTracklessCuts(true, true, true, true, true,
                                                 true);
          /*if (event.m_is_signal && (universe->ShortName() == "cv" ||
          universe->ShortName() == "CCPi+ Tune")){
            int idx = universe->GetHighestEnergyTruePionIndex();
            std::cout << "Event = " << i_event
                      << " Universe " <<universe->ShortName()<<" Tpi weight = "
            << universe->GetChargedPionTuneWeight(universe->GetTpiTrue(idx),
          universe->GetQ2True()) << " IsSignal? = " << event.m_is_signal <<"\n";
            std::cout << "Weight = " << universe->GetWeight() << "\n";
            std::cout << "Q^2 = " << universe->GetQ2True() << "\n";
            std::cout << "tpi = " << universe->GetTpiTrue(idx) << "\n";
            std::cout << "untraked tpi = " << universe->GetTrueTpi() << "\n";
          }*/
          //          if (event.m_is_signal) std::cout << "Event = " << i_event
          //          << " Q2 = " <<
          //		   universe->GetQ2True()/1000000 << " Weight ="
          //                   << universe->GetWeight() << "\n";
          ccpi_event::FillTruthEvent(event, variables);
        }
      }
    } else {
      LowRecoilPion::Cluster d;
      LowRecoilPion::Cluster c(*cvUniv, 0);
      LowRecoilPion::Michel<CVUniverse> m(*cvUniv, 0);
      LowRecoilPion::MichelEvent<CVUniverse> trackless_michels;
      bool good_trackless_michels;
      if (onlytracked) {
        good_trackless_michels = false;
      } else {
        good_trackless_michels =
            LowRecoilPion::hasMichel<CVUniverse,
                                     LowRecoilPion::MichelEvent<CVUniverse>>::
                hasMichelCut(*cvUniv, trackless_michels);
        // good_trackless_michels = BestMichelDistance2DCut(*universe,
        // trackless_michels);
        good_trackless_michels =
            good_trackless_michels &&
            LowRecoilPion::BestMichelDistance2D<
                CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::
                BestMichelDistance2DCut(*cvUniv, trackless_michels);
        // good_trackless_michels = MichelRangeCut(*universe,
        // trackless_michels);
        good_trackless_michels =
            good_trackless_michels &&
            LowRecoilPion::GetClosestMichel<
                CVUniverse, LowRecoilPion::MichelEvent<CVUniverse>>::
                GetClosestMichelCut(*cvUniv, trackless_michels);
      }
      // Loop universes, make cuts, and fill
      for (auto error_band : error_bands) {
        std::vector<CVUniverse*> universes = error_band.second;
        for (auto universe : universes) {
          universe->SetEntry(i_event);
          // std::cout << universe->ShortName() << "\n";
          // if (universe->GetDouble("mc_incoming") == 12 &&
          //    universe->ShortName() == "cv")
          //  universe->PrintArachneLink();

          // CCPiEvent keeps track of lots of event properties
          CCPiEvent event(is_mc, is_truth, signal_definition, universe);
          universe->SetIsSignal(event.m_is_signal);
          event.m_weight = universe->GetWeight();

          //===============
          // CHECK CUTS
          //===============

          universe->SetVtxMichels(trackless_michels);
          bool pass = true;
          pass = pass && universe->GetNMichels() == 1;
          pass =
              pass && universe->GetTpiTrackless() > signal_definition.m_tpi_min;
          pass =
              pass && universe->GetTpiTrackless() < signal_definition.m_tpi_max;
          pass = pass && universe->GetPmu() > signal_definition.m_PmuMinCutVal;
          pass = pass && universe->GetPmu() < signal_definition.m_PmuMaxCutVal;
          pass = pass &&
                 universe->GetNIsoProngs() < signal_definition.m_IsoProngCutVal;
          pass =
              pass && universe->IsInHexagon(universe->GetVecElem("vtx", 0),
                                            universe->GetVecElem("vtx", 1),
                                            signal_definition.m_ApothemCutVal);
          pass = pass && universe->GetVecElem("vtx", 2) >
                             signal_definition.m_ZVtxMinCutVal;
          pass = pass && universe->GetVecElem("vtx", 2) <
                             signal_definition.m_ZVtxMaxCutVal;
          // pass = pass && universe->GetBool("isMinosMatchTrack");
          pass = pass && universe->GetInt("isMinosMatchTrack") == 1;
          pass = pass && universe->GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
          pass =
              pass && universe->GetThetamu() < signal_definition.m_thetamu_max;
          pass = pass && universe->GetPTmu() < signal_definition.m_ptmu_max;
          pass =
              pass && universe->GetTracklessWexp() > signal_definition.m_w_min;
          // implementing multipion cut
          int unique_michel_idx_untracked = -1;
          if (trackless_michels.m_idx != -1) {
            unique_michel_idx_untracked =
                trackless_michels.m_nmichels[trackless_michels.m_idx].tuple_idx;
          }
          LowRecoilPion::MichelEvent<CVUniverse> dummy_trackless_michel =
              event.m_universe->GetVtxMichels();
          //          std::cout << "trackless_michels.m_idx = " <<
          //          trackless_michels.m_idx
          //                    << " event.m_universe.m_vtx_michels.m_idx  = "
          //                    << dummy_trackless_michel.m_idx << "\n";
          std::vector<int> unique_michel_idx_tracked;
          endpoint::MichelMap tracked_michels = GetTrackedPionCandidates(event);
          for (auto candidate : tracked_michels) {
            unique_michel_idx_tracked.push_back(candidate.first);
          }
          //	for (int i = 0; i < (int)unique_michel_idx_tracked.size(); ++i)
          //          std::cout << "Untracked index = " <<
          //          unique_michel_idx_untracked << " Tracked index = " <<
          //          unique_michel_idx_tracked[i] << "\n";

          //===============
          // CHECK CUTS
          //===============
          // Universe only affects weights
          // Check Cuts -- computationally expensive
          //
          // This looks complicated for optimization reasons.
          // Namely, for all vertical-only universes (meaning only the event
          // weight differs from CV) no need to recheck cuts.
          PassesCutsInfo cuts_info;
          if (universe->IsVerticalOnly()) {
            if (!checked_cv) {
              cv_cuts_info = PassesCuts(event);
              checked_cv = true;
            }
            assert(checked_cv);
            cuts_info = cv_cuts_info;
          } else {
            cuts_info = PassesCuts(event);
          }

          // Save results of cuts to Event and universe
          std::tie(event.m_passes_cuts, event.m_is_w_sideband,
                   event.m_passes_all_cuts_except_w,
                   event.m_reco_pion_candidate_idxs) = cuts_info.GetAll();

          event.m_highest_energy_pion_idx =
              GetHighestEnergyPionCandidateIndex(event);

          universe->SetPionCandidates(event.m_reco_pion_candidate_idxs);

          // Re-call GetWeight because the node cut efficiency systematic
          // needs a pion candidate to calculate its weight.
          event.m_weight = universe->GetWeight();
          /*
          if ((good_trackless_michels && pass) || event.m_passes_cuts ||
              event.m_is_w_sideband || event.m_passes_all_cuts_except_w) {
            std::cout << "pass = " << good_trackless_michels << "\n";
            std::cout << "event.m_passes_cuts = " << event.m_passes_cuts
                      << "\n";
            std::cout << "event.m_is_w_sideband = " << event.m_is_w_sideband
                      << "\n";
            std::cout << "event.m_passes_all_cuts_except_w = "
                      << event.m_passes_all_cuts_except_w << "\n";
            std::cout << "Untracked index = " << unique_michel_idx_untracked
                      << "\n";
            for (int i = 0; i < (int)unique_michel_idx_tracked.size(); ++i)
              std::cout << "Untracked index = " << unique_michel_idx_untracked
                        << " Tracked index = " << unique_michel_idx_tracked[i]
                        << "\n";
          }
          */
          // These conditions are used to make the tracked or untracked dta
          // selection
          if (onlyuntracked) {
            event.m_passes_cuts = false;
            event.m_is_w_sideband = false;
            event.m_passes_all_cuts_except_w = false;
          }
          if (onlytracked) {
            good_trackless_michels = good_trackless_michels && false;
            pass = pass && false;
          }
          universe->SetVtxMichels(trackless_michels);
          event.m_passes_trackless_cuts_except_w = pass;
          event.m_passes_trackless_sideband = false;
          if (pass &&
              universe->GetTracklessWexp() > signal_definition.m_w_max) {
            if (universe->GetTracklessWexp() >= sidebands::kSidebandCutVal)
              event.m_passes_trackless_sideband = true;
            pass = false;
          }
          event.m_passes_trackless_cuts = good_trackless_michels && pass;
          event.m_passes_trackless_sideband =
              event.m_passes_trackless_sideband && good_trackless_michels;
          event.m_passes_trackless_cuts_except_w =
              event.m_passes_trackless_cuts_except_w && good_trackless_michels;
          universe->SetPassesTrakedTracklessCuts(
              event.m_passes_cuts, event.m_passes_trackless_cuts,
              event.m_is_w_sideband, event.m_passes_trackless_sideband,
              event.m_passes_all_cuts_except_w,
              event.m_passes_trackless_cuts_except_w);
          /*if (event.m_is_signal){
            std::cout << "Event = " << i_event << " q2 = " <<
          universe->GetQ2True()/1000000
                      << " Weight = " << universe->GetWeight() << "\n";
          }*/
          /*if (event.m_passes_cuts && (universe->ShortName() == "cv" ||
          universe->ShortName() == "CCPi+ Tune")){
            int idx = universe->GetHighestEnergyTruePionIndex();
            std::cout << "Event = " << i_event
                      << " Universe " << universe->ShortName() << " Tpi = " <<
                      universe->GetMixedTpi(event.m_highest_energy_pion_idx) <<
          "\n"; std::cout << "Event = " << i_event
                      << " Universe " <<universe->ShortName()<<" Tpi weight = "
           << universe->GetChargedPionTuneWeight(universe->GetTpiTrue(idx),
          universe->GetQ2True()) << " IsSignal? = " << event.m_is_signal <<"\n";
            std::cout << "Weight = " << universe->GetWeight() << "\n";
            std::cout << "Q^2 = " << universe->GetQ2True() << "\n";
            std::cout << "tpi = " << universe->GetTpiTrue(idx) << "\n";
            std::cout << "untraked tpi = " << universe->GetTrueTpi() << "\n";
          }*/
          /*
          if (event.m_passes_cuts){
            std::cout << i_event << "  " << universe->GetInt("mc_run")
            << "  " << universe->GetInt("mc_subrun") << "  " <<
            universe->GetInt("mc_nthEvtInFile") + 1 << "  "  <<
            universe->GetVecElem("slice_numbers", 0) << "  " <<
            universe->GetDouble("iso_prongs_count") << "  " <<
            universe->GetDouble("n_nonvtx_iso_blobs_all") << "  " <<
            event.m_is_signal << "\n"; selcount++;
          }
          */
          //        std::cout << "Event = " << i_event << "\n";

          /*

                    //if (event.m_passes_trackless_cuts ||
                    //event.m_passes_trackless_cuts_except_w ||
                    //event.m_passes_all_cuts_except_w) {
                    if (event.m_is_signal) {
                    //if (event.m_passes_trackless_cuts) {
                    if (universe->GetWexp() > 1400){ std::cout << "Event = "
                    << i_event << "\n"; std::cout << "event.m_passes_cuts = " <<
                      event.m_passes_cuts << "\n"; std::cout <<
                      "event.m_passes_all_cuts_except_w = " <<
                      event.m_passes_all_cuts_except_w << "\n"; std::cout <<
                      "event.m_passes_trackless_cuts = " <<
          event.m_passes_trackless_cuts
                    <<
                      "\n"; std::cout << "event.m_passes_trackless_cuts_except_w
          = " << event.m_passes_trackless_cuts_except_w << "\n"; std::cout <<
                    "GetWexp = " << universe->GetWexp() << "\n"; std::cout <<
                    "GetTracklessWexp = "
                      << universe->GetTracklessWexp() << "\n"; std::cout <<
                    "GetTrackedWexp = " << universe->GetTrackedWexp() << "\n";
                      //            std::cout << "GetQ2 = " << universe->GetQ2()
          << "\n";
                      //            std::cout <<
                                  std::cout << "Is signal = " <<
          event.m_is_signal <<
                    "\n";
                                }
                      //          std::cout <<
          "event.m_passes_trackless_sideband = " <<
          //            event.m_passes_trackless_sideband << "\n";
                      //          std::cout << "event.m_is_w_sideband = " <<
          //            event.m_is_w_sideband << "\n";
                      //          std::cout << "Tpi mixed = " <<
          //            universe->GetMixedTpi(event.m_highest_energy_pion_idx)
          << "\n";
                      //          std::cout << "Tpi tracked = " <<
          //            universe->GetTpi(event.m_highest_energy_pion_idx) <<
          "\n";
                      //          std::cout << "Tpi tracked true = " <<
          // universe->GetTpiTrue(universe->GetHighestEnergyTruePionIndex()) <<
          //            "\n";
                      //          std::cout << "Tpi trackless = " <<
          //            universe->GetTpiTrackless() << "\n";
                      //          std::cout << "Tpi trackless true = " <<
          //            universe->GetTrueTpi() << "\n";
                      //          std::cout << "trackless_michels.m_bestdist = "
          <<
          //            trackless_michels.m_bestdist << "\n";
                      //        std::cout << "universe->m_vtx_michels.m_bestdist
          = " <<
          //            universe->GetBestDistance() << "\n";
                      //        std::cout <<
          //            "universe->GetTpiUntracked(trackless_michels.m_bestdist)
          = " <<
          //            universe->GetTpiUntracked(trackless_michels.m_bestdist)
          << "\n"; std::cout << "universe->GetTpiTrackless() = " <<
                      universe->GetTpiTrackless() << "\n";
                                universe->PrintArachneLink();
                      //          std::cout << "thetapi reco = " <<
          //            universe->GetThetapitracklessDeg() << "\n";
                      //          std::cout << "thetapi truth = " <<
          //            universe->GetThetapitracklessTrueDeg() << "\n";
                    }
               */
          //===============
          // FILL RECO
          //===============
          /*if (event.m_passes_cuts || event.m_passes_trackless_cuts){
          //if (i_event == 7002 || i_event == 8789){
            std::cout << i_event << ", ";
            std::cout << "Event = " << i_event << "\n";
            std::cout << "Pass tracked cuts? " << event.m_passes_cuts << "\n";
            std::cout << "Pass the untracked cuts? " <<
          event.m_passes_trackless_cuts << "\n"; std::cout << "Pass is tracked
          sideband? " << event.m_is_w_sideband << "\n"; std::cout << "Pass is
          untracked sideband? " << event.m_passes_trackless_sideband << "\n";
            std::cout << "Pass cuts exept W tracked " <<
          event.m_passes_all_cuts_except_w << "\n"; std::cout << "Pass cuts
          exept W untracked " << event.m_passes_trackless_cuts_except_w << "\n";
            std::cout << "Tracked Wexp = " << universe->GetTrackedWexp() <<
          "\n"; std::cout << "Emu = " << universe->GetEmu() << " CalRecoilE = "
          << universe->GetCalRecoilEnergy() << " TrackRecE = " <<
          universe->GetTrackRecoilEnergy() << " thetamu =" <<
          universe->GetThetamu() << "\n";
          }*/
          ccpi_event::FillRecoEvent(event, variables);
        }  // universes
      }    // error bands
    }      // RECO
  }        // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void makeCrossSectionMCInputs(int signal_definition_int = 0,
                              std::string plist = "ME1A",
                              bool do_systematics = false,
                              bool do_truth = false,
                              const bool do_test_playlist = false,
                              bool is_grid = false, std::string input_file = "",
                              int run = 0) {
  // 1. Input Data
  std::string mc_file_list;
  const bool is_mc = true;
  const bool use_xrootd = true;
  assert(!(is_grid && input_file.empty()) &&
         "On the grid, infile must be specified.");
  mc_file_list = input_file.empty() ? GetPlaylistFile(plist, is_mc, use_xrootd)
                                    : input_file;

  // INIT MACRO UTILITY
  const std::string macro("MCXSecInputs");
  // std::string a_file =
  // "root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr/minerva/persistent/users/bmesserl/pions//20200713/merged/mc/ME1A/CCNuPionInc_mc_AnaTuple_run00110000_Playlist.root";
  // 2. Macro Utility/Manager -- keep track of macro-level things, creat input
  // data chain
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, plist, do_truth,
                       is_grid, do_systematics);
  util.m_name = "MCXSecInputs";
  util.PrintMacroConfiguration();

  // 3. Prepare Output
  std::string outfile_name = GetOutFilename(util, run);
  std::cout << "Saving output to " << outfile_name << "\n\n";
  TFile fout(outfile_name.c_str(), "RECREATE");

  // 4. Initialize Variables (and the histograms that they own)
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  // 5. Loop MC Reco -- process events and fill histograms owned by variables
  bool is_truth = false;
  LoopAndFillMCXSecInputs(util.m_error_bands, util.GetMCEntries(), is_truth,
                          util.m_signal_definition, variables);

  // 6. Loop Truth
  if (util.m_do_truth) {
    is_truth = true;
    LoopAndFillMCXSecInputs(util.m_error_bands_truth, util.GetTruthEntries(),
                            is_truth, util.m_signal_definition, variables);
  }

  // 7. Write to file
  std::cout << "Synching and Writing\n\n";
  WritePOT(fout, is_mc, util.m_mc_pot);
  fout.cd();
  for (auto v : variables) {
    SyncAllHists(*v);
    v->WriteMCHists(fout);
    /*    if (util.m_do_truth && v->m_is_true){
          SavingStacked(fout, v->GetStackArray(kOtherInt), v->Name(), "FSP");
          SavingStacked(fout, v->GetStackArray(kCCQE), v->Name(), "Int");
          SavingStacked(fout, v->GetStackArray(kPim), v->Name(), "Hadrons");
          SavingStacked(fout, v->GetStackArray(kOnePion), v->Name(), "Npi");
          SavingStacked(fout, v->GetStackArray(kOnePi0), v->Name(), "Npi0");
          SavingStacked(fout, v->GetStackArray(kOnePip), v->Name(), "Npip");
          SavingStacked(fout, v->GetStackArray(kWSideband_Low), v->Name(),
       "WSB"); SavingStacked(fout, v->GetStackArray(kB_Meson), v->Name(),
       "Msn"); SavingStacked(fout, v->GetStackArray(kB_HighW), v->Name(),
       "WBG");
        }*/
  }
}

#endif  // makeXsecMCInputs_C
