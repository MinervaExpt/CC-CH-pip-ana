#ifndef plotCrossSectionFromFile_C
#define plotCrossSectionFromFile_C

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMatrixD.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif  // MNVROOT6

#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/TargetUtils.h"
#include "includes/Constants.h"
#include "includes/MacroUtil.h"
#include "includes/Plotter.h"
#include "includes/Variable.h"
#include "includes/common_functions.h"
#include "makeCrossSectionMCInputs.C"  // GetAnalysisVariables
#include "plotting_functions.h"

void SetPOT(TFile& fin, CCPi::MacroUtil& util) {
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
void plotCrossSectionFromFile(int signal_definition_int = 6,
                              int plot_errors = 1) {
  // Infiles
  TFile fin(
      //      "DataXSecInputs_20241204_ALL_mix_CCPiWeightSys_p4.root",
      "DataXSecInputs_20241217_ALL_mix_thetapi_CCPiWeighSys_p4.root", "READ");
  // TFile fin("DataXSecInputs_20240819_ME1A_mixed_NewSys_sys_p4.root", "READ");
  TFile fin1("GENIEXSECEXTRACT_AarSignalDefMCME1A_q2.root", "READ");
  TFile fin2("GENIEXSECEXTRACT_AarSignalDefTpichangeMCME1A_q2.root", "READ");
  TFile fin3(
      "/minerva/data/users/abercell/hists/xsec/"
      "xsec_new_jeffrey_flux_MENU1PI_plastic_MinervaME1ABCDEFGLMNOP.root",
      "READ");
  TFile fin4(
      "/minerva/data/users/abercell/hists/xsec_inputs/"
      "Merge_BkgdSub_Unfold_MENU1PI_POTNorm_plastic_MinervaME1ABCDEFGLMNOP."
      "root",
      "READ");
  //  TFile
  //  fin5("/minerva/data/users/abercell/hists/Macro/GridOneLoop_MENU1PI_MinosMatched_plastic_Merged_NewdEdXCal_MinervaME1ABCDEFGLMNOP_Data_Merged_NewdEdXCal_Tracker_MinervaME1ABCDEFGLMNOP_MC.root",
  //  "READ");
  TFile fin5(
      "/minerva/data/users/abercell/hists/xsec_inputs/"
      "Merge_Eff_MENU1PI_POTNorm_plastic_MinervaME1ABCDEFGLMNOP.root",
      "READ");
  TFile fin6("DataXSecInputs_20240318_ALL_AaronSignalDef_nosys_p4.root",
             "READ");
  TFile fin7("DataXSecInputs_20240325_ALL_mixed_NoTpiweight_nosys_p4.root",
             "READ");
  //  TFile
  //  fin7("DataXSecInputs_20240304_ALL_AaronSigDef_plusAaronFidVolCuts_nosys_p4.root",
  //  "READ");
  TFile fin8(
      "/minerva/data/users/abercell/hists/Macro/"
      "GridOneLoop_MENU1PI_MinosMatched_plastic_Merged_NewdEdXCal_"
      "MinervaME1ABCDEFGLMNOP_Data_Merged_NewdEdXCal_Tracker_"
      "MinervaME1ABCDEFGLMNOP_MC.root",
      "READ");
  TFile fin9("DataXSecInputs_20240422_ALL_untracked_NewEstimator_noSys_p4.root",
             "READ");
  TFile fin10("MCXSecInputs_20240417_me1P_AaronSigDef_PEMA_nosys_p4.root",
              "READ");
  // TFile
  // finaux("MCXSecInputs_20240804_ALL_mixed_NoNMichelsCut_nosys_p4_NOMINAL.root",
  // "READ");
  TFile finaux("DataXSecInputs_20241204_ALL_mix_CCPiCVandSys_p4.root", "READ");

  cout << "Reading input from " << fin.GetName() << endl;

  // Set up macro utility object...which gets the list of systematics for us...
  // which we need in order to read in HistWrappers...which we don't need at
  // this point...indeed we only need MnvH1Ds...so that's a TODO: write a
  // function that only loads in MnvH1D's from a file, not HWs.
  // Most of these options aren't used in this script. TODO make a CTOR that
  // doesn't require them.
  // INPUT TUPLES
  // Don't actually use the MC chain, only load it to indirectly access it's
  // systematics
  const std::string plist = "ME1L";
  std::string data_file_list = GetPlaylistFile(plist, false, false);
  std::string mc_file_list = GetPlaylistFile(plist, true, false);
  // std::string data_file_list = GetTestPlaylist(false);
  // std::string mc_file_list = GetTestPlaylist(true);

  // Macro Utility
  const std::string macro("PlotCrossSectionFromFile");
  bool do_truth = false, is_grid = false, do_systematics = true;
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
                       plist, do_truth, is_grid, do_systematics);
  // Get POT from file, not from any chain
  SetPOT(fin, util);
  util.PrintMacroConfiguration(macro);

  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);
  ContainerEraser::erase_if(variables, [](Variable* v) {
    return v->Name() == "tpi_res";  // || v->Name() == "wexp_fit";
  });

  PlotUtils::MnvH1D* AaronXsecGENIE = (PlotUtils::MnvH1D*)fin1.Get("q2_xsec");
  PlotUtils::MnvH1D* BenXsecGENIE = (PlotUtils::MnvH1D*)fin2.Get("q2_xsec");

  PlotUtils::MnvH1D* AaronUnfoldData =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_dataUnfold2");
  PlotUtils::MnvH1D* AaronCorrectedData =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_data_eff_corrected");
  PlotUtils::MnvH1D* AaronEffDen =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_pi_channel_truth_sig");
  PlotUtils::MnvH1D* AaronXsecMC =
      (PlotUtils::MnvH1D*)fin3.Get("h_q2_plastic_pi_channel_mc_xsec_nucleon");
  PlotUtils::MnvH1D* AaronXsecData =
      (PlotUtils::MnvH1D*)fin3.Get("h_q2_plastic_data_xsec_nucleon");
  PlotUtils::MnvH1D* AaronEffNum =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_pi_channel_signal");
  PlotUtils::MnvH1D* AaronEff =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_pi_channel_efficiency");
  PlotUtils::MnvH1D* AaronNewEff =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_pi_channel_efficiency");
  PlotUtils::MnvH1D* AaronDataSel =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_data");
  PlotUtils::MnvH1D* AaronBGSub =
      (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_BkgdSubData");
  PlotUtils::MnvH1D* AaronEffOtherFile =
      (PlotUtils::MnvH1D*)fin5.Get("h_q2_plastic_pi_channel_efficiency");
  PlotUtils::MnvH1D* AaronEffNumOtherFile =
      (PlotUtils::MnvH1D*)fin5.Get("h_q2_plastic_EFF_NUM_pi_channel");
  PlotUtils::MnvH1D* AaronEffDenOtherFile =
      //                 (PlotUtils::MnvH1D*)fin4.Get("h_q2_plastic_pi_channel_truth_sig_ME1A");
      (PlotUtils::MnvH1D*)fin5.Get("h_q2_plastic_EFF_DEN_pi_channel");
  PlotUtils::MnvH1D* AaronFluxNormalizer =
      (PlotUtils::MnvH1D*)fin3.Get("h_q2_plastic_data_xsec_nucleon");
  PlotUtils::MnvH1D* AaronRevUnfoldData =
      (PlotUtils::MnvH1D*)fin3.Get("h_q2_plastic_data_xsec_nucleon");
  PlotUtils::MnvH1D* h_data_POT = (PlotUtils::MnvH1D*)fin3.Get("h_Data_POT");
  PlotUtils::MnvH1D* h_MC_POT = (PlotUtils::MnvH1D*)fin3.Get("h_MC_POT");
  PlotUtils::MnvH1D* AaronNoNormXsecData =
      (PlotUtils::MnvH1D*)fin3.Get("h_q2_plastic_data_xsec_nucleon");

  PlotUtils::MnvH1D* dummyBenXsecData =
      (PlotUtils::MnvH1D*)fin.Get("cross_section_q2");
  PlotUtils::MnvH1D* dummyBenXsecDataCorr =
      (PlotUtils::MnvH1D*)fin.Get("cross_section_q2");
  PlotUtils::MnvH1D* dummyBenEffCorrData =
      (PlotUtils::MnvH1D*)fin.Get("efficiency_corrected_data_q2");
  PlotUtils::MnvH1D* dummyBenEffDen =
      (PlotUtils::MnvH1D*)fin.Get("effden_q2_true");
  PlotUtils::MnvH1D* dummyBenUnfoldData =
      (PlotUtils::MnvH1D*)fin.Get("unfolded_q2");
  PlotUtils::MnvH1D* dummyBenXsecMC =
      (PlotUtils::MnvH1D*)fin.Get("mc_cross_section_q2");
  PlotUtils::MnvH1D* dummyBenEff = (PlotUtils::MnvH1D*)fin.Get("efficiency_q2");
  PlotUtils::MnvH1D* dummyBenEffNum = (PlotUtils::MnvH1D*)fin.Get("effnum_q2");
  PlotUtils::MnvH1D* dummyBenDataSel =
      (PlotUtils::MnvH1D*)fin.Get("selection_data_q2");
  PlotUtils::MnvH1D* dummyBenBGSub =
      (PlotUtils::MnvH1D*)fin.Get("bg_subbed_data_q2");
  PlotUtils::MnvH1D* dummyBenXsecAaronSigData =
      (PlotUtils::MnvH1D*)fin6.Get("cross_section_q2");
  PlotUtils::MnvH1D* dummyBenXsecNoTpiWegihtData =
      (PlotUtils::MnvH1D*)fin7.Get("cross_section_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigCutsEffDen =
      (PlotUtils::MnvH1D*)fin7.Get("effden_q2_true");
  PlotUtils::MnvH1D* dummyBenAaronSigCutsEffNum =
      (PlotUtils::MnvH1D*)fin7.Get("effnum_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigCutsEff =
      (PlotUtils::MnvH1D*)fin7.Get("efficiency_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigCutsEffCorr =
      (PlotUtils::MnvH1D*)fin7.Get("efficiency_corrected_data_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigCutsUnfold =
      (PlotUtils::MnvH1D*)fin7.Get("unfolded_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigCutsDataSel =
      (PlotUtils::MnvH1D*)fin7.Get("selection_data_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigCutsBGSub =
      (PlotUtils::MnvH1D*)fin7.Get("bg_subbed_data_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigEffDen =
      (PlotUtils::MnvH1D*)fin6.Get("effden_q2_true");
  PlotUtils::MnvH1D* dummyBenAaronSigEffNum =
      (PlotUtils::MnvH1D*)fin6.Get("effnum_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigEff =
      (PlotUtils::MnvH1D*)fin6.Get("efficiency_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigEffCorr =
      (PlotUtils::MnvH1D*)fin6.Get("efficiency_corrected_data_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigUnfold =
      (PlotUtils::MnvH1D*)fin6.Get("unfolded_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigDataSel =
      (PlotUtils::MnvH1D*)fin6.Get("selection_data_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigBGSub =
      (PlotUtils::MnvH1D*)fin6.Get("bg_subbed_data_q2");
  PlotUtils::MnvH1D* dummyBenXsecUntracked =
      (PlotUtils::MnvH1D*)fin9.Get("cross_section_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigEffDenOldMacro =
      (PlotUtils::MnvH1D*)fin10.Get("effden_q2_true");
  PlotUtils::MnvH1D* dummyBenAaronSigEffNumOldMacro =
      (PlotUtils::MnvH1D*)fin10.Get("effnum_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigEffOldMacro =
      (PlotUtils::MnvH1D*)fin9.Get("efficiency_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigEffCorrOldMacro =
      (PlotUtils::MnvH1D*)fin9.Get("efficiency_corrected_data_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigUnfoldOldMacro =
      (PlotUtils::MnvH1D*)fin9.Get("unfolded_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigDataSelOldMacro =
      (PlotUtils::MnvH1D*)fin9.Get("selection_data_q2");
  PlotUtils::MnvH1D* dummyBenAaronSigBGSubOldMacro =
      (PlotUtils::MnvH1D*)fin9.Get("bg_subbed_data_q2");
  PlotUtils::MnvH1D* dataPOTBenAaronSig =
      (PlotUtils::MnvH1D*)fin6.Get("data_pot");
  PlotUtils::MnvH1D* MCPOTBenAaronSig = (PlotUtils::MnvH1D*)fin6.Get("mc_pot");

  PlotUtils::MnvH1D* dataPOTBenAaronSigCuts =
      (PlotUtils::MnvH1D*)fin7.Get("data_pot");
  PlotUtils::MnvH1D* MCPOTBenAaronSigCuts =
      (PlotUtils::MnvH1D*)fin7.Get("mc_pot");

  PlotUtils::MnvH1D* dataPOTBenAaronSigOldMacro =
      (PlotUtils::MnvH1D*)fin9.Get("data_pot");
  PlotUtils::MnvH1D* MCPOTBenAaronSigOldMacro =
      (PlotUtils::MnvH1D*)fin10.Get("mc_pot");

  PlotUtils::MnvH1D* AaronNewEffCorr =
      (PlotUtils::MnvH1D*)AaronUnfoldData->Clone("EffCorrAaron");

  PlotUtils::MnvH1D* BenXsecData =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecData");
  BenXsecData->SetTitle("BenXsecData");
  PlotUtils::MnvH1D* BenXsecDataCorr =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecDataCorr");
  PlotUtils::MnvH1D* BenEffCorrData =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecDataCorr");
  BenEffCorrData->SetTitle("BenEffCorr");
  PlotUtils::MnvH1D* BenEffDen =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecData");
  BenEffDen->SetTitle("BenEffDen");
  PlotUtils::MnvH1D* BenUnfoldData =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenUnfold");
  BenUnfoldData->SetTitle("BenUnfold");
  PlotUtils::MnvH1D* BenDataSel =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenDataSel");
  BenDataSel->SetTitle("BenDataSel");
  PlotUtils::MnvH1D* BenBGSub =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenBGSub");
  BenBGSub->SetTitle("BenBGSub");
  PlotUtils::MnvH1D* BenXsecMC =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecDataCorr");
  PlotUtils::MnvH1D* BenEff =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenEff");
  BenEff->SetTitle("BenEff");
  PlotUtils::MnvH1D* BenEffNum =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenEffNum");
  BenEffNum->SetTitle("BenEffNum");
  PlotUtils::MnvH1D* BenXsecAaronSigData =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecAaronSigData");
  BenXsecAaronSigData->SetTitle("BenXsecAaronSigData");
  PlotUtils::MnvH1D* BenAaronSigEffDen =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEffDen");
  BenAaronSigEffDen->SetTitle("BenAaronSigEffDen");
  PlotUtils::MnvH1D* BenAaronSigEffNum =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEffNum");
  BenAaronSigEffNum->SetTitle("Ben w/ Aaron Sig Def");
  PlotUtils::MnvH1D* BenAaronSigEff =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEff");
  BenAaronSigEff->SetTitle("BenAaronSigEff");
  PlotUtils::MnvH1D* BenAaronSigEffCorr =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEffCorr");
  BenAaronSigEffCorr->SetTitle("BenAaronSigEffCorr");
  PlotUtils::MnvH1D* BenAaronSigUnfold =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigUnfold");
  BenAaronSigUnfold->SetTitle("Ben w/ Aaron Sig Def");

  PlotUtils::MnvH1D* BenAaronSigDataSel =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigDataSel");
  BenAaronSigDataSel->SetTitle("BenAaronSigDataSel");
  PlotUtils::MnvH1D* BenAaronSigBGSub =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigBGSub");
  BenAaronSigBGSub->SetTitle("BenAaronSigBGSub");
  PlotUtils::MnvH1D* BenXsecNoTpiWegihtData =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecAaronSigCuts");
  BenXsecNoTpiWegihtData->SetTitle("Ben No T_{#pi} weight");
  PlotUtils::MnvH1D* BenAaronSigCutsEffDen =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigCutsEffDen");
  BenAaronSigCutsEffDen->SetTitle("Ben No T_{#pi} weight");
  PlotUtils::MnvH1D* BenAaronSigCutsEffNum =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigCutsEffNum");
  BenAaronSigCutsEffNum->SetTitle("Ben No T_{#pi} weight");
  PlotUtils::MnvH1D* BenAaronSigCutsEff =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigCutsEff");
  BenAaronSigCutsEff->SetTitle("B/E SDef NotpiW Eff");
  PlotUtils::MnvH1D* BenAaronSigCutsEffCorr =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigCutsEffCorr");
  BenAaronSigCutsEffCorr->SetTitle("BenAaronSigCutsEffCorr");
  PlotUtils::MnvH1D* BenAaronSigCutsUnfold =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigCutsUnfold");
  BenAaronSigCutsUnfold->SetTitle("BenAaronSigCutsUnfold");

  PlotUtils::MnvH1D* BenAaronSigCutsDataSel =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigCutsDataSel");
  BenAaronSigCutsDataSel->SetTitle("BenAaronSigCutsDataSel");
  PlotUtils::MnvH1D* BenAaronSigCutsBGSub =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigCutsBGSub");
  BenAaronSigCutsBGSub->SetTitle("BenAaronSigCutsBGSub");

  PlotUtils::MnvH1D* BenXsecUntracked =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenXsecUntracked");
  BenXsecUntracked->SetTitle("Untracked");
  PlotUtils::MnvH1D* BenAaronSigEffDenOldMacro =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEffDenOldMacro");
  BenAaronSigEffDenOldMacro->SetTitle("PEMA on");
  PlotUtils::MnvH1D* BenAaronSigEffNumOldMacro =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEffNumOldMacro");
  BenAaronSigEffNumOldMacro->SetTitle("PEMA on");
  PlotUtils::MnvH1D* BenAaronSigEffOldMacro =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEffOldMacro");
  BenAaronSigEffOldMacro->SetTitle("BenSigUntrackedEff");
  PlotUtils::MnvH1D* BenAaronSigEffCorrOldMacro =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigEffCorrOldMacro");
  BenAaronSigEffCorrOldMacro->SetTitle("BenSigUntrackedEffCorr");
  PlotUtils::MnvH1D* BenAaronSigUnfoldOldMacro =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigUnfoldOldMacro");
  BenAaronSigUnfoldOldMacro->SetTitle("BenSigUntrackedUnfold");
  PlotUtils::MnvH1D* BenAaronSigDataSelOldMacro =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigDataSelOldMacro");
  BenAaronSigDataSelOldMacro->SetTitle("BenSigUntrackedDataSel");
  PlotUtils::MnvH1D* BenAaronSigBGSubOldMacro =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("BenAaronSigBGSubOldMacro");
  BenAaronSigBGSubOldMacro->SetTitle("BenSigUntrackedBGSub");
  PlotUtils::MnvH1D* EffPEMA =
      (PlotUtils::MnvH1D*)AaronXsecData->Clone("EffPEMA");
  EffPEMA->SetTitle("PEMA Eff");

  BenXsecData->Reset();
  BenXsecDataCorr->Reset();
  BenEffDen->Reset();
  BenUnfoldData->Reset();
  BenEffCorrData->Reset();
  BenXsecMC->Reset();
  BenEff->Reset();
  BenEffNum->Reset();
  BenDataSel->Reset();
  BenBGSub->Reset();
  AaronNewEffCorr->Reset();
  AaronNewEff->Reset();
  AaronFluxNormalizer->Reset();
  BenXsecAaronSigData->Reset();
  BenXsecNoTpiWegihtData->Reset();
  BenAaronSigEffDen->Reset();
  BenAaronSigEffNum->Reset();
  BenAaronSigEff->Reset();
  BenAaronSigEffCorr->Reset();
  BenAaronSigUnfold->Reset();
  BenAaronSigDataSel->Reset();
  BenAaronSigBGSub->Reset();
  BenAaronSigCutsEffDen->Reset();
  BenAaronSigCutsEffNum->Reset();
  BenAaronSigCutsEff->Reset();
  BenAaronSigCutsEffCorr->Reset();
  BenAaronSigCutsUnfold->Reset();
  BenAaronSigCutsDataSel->Reset();
  BenAaronSigCutsBGSub->Reset();
  BenAaronSigEffDenOldMacro->Reset();
  BenAaronSigEffNumOldMacro->Reset();
  BenAaronSigEffOldMacro->Reset();
  BenAaronSigEffCorrOldMacro->Reset();
  BenAaronSigUnfoldOldMacro->Reset();
  BenAaronSigDataSelOldMacro->Reset();
  BenAaronSigBGSubOldMacro->Reset();
  BenXsecUntracked->Reset();
  EffPEMA->Reset();

  int flux_nuPDG = 14;
  std::string flux_playlist = "minervame1d1m1nWeightedAve";
  AaronFluxNormalizer =
      flux_reweighter(flux_playlist, flux_nuPDG, true)
          .GetIntegratedTargetFlux(flux_nuPDG, "tracker", AaronFluxNormalizer,
                                   0., 100., "targets_12345_jointNueIMD");
  //  for (int i = 1; i <= BenXsecData->GetNbinsX(); ++i)
  //    AaronFluxNormalizer->SetBinContent(i, 6.322e-8);
  AaronFluxNormalizer->Scale(1.0e-4);

  PlotUtils::MnvH1D* h_flux_normalization =
      (PlotUtils::MnvH1D*)AaronCorrectedData->Clone("flux_normalization");
  PlotUtils::MnvH1D* h_efficiency_corrected_data =
      (PlotUtils::MnvH1D*)AaronCorrectedData->Clone(
          "efficiency_corrected_data");
  h_flux_normalization->ClearAllErrorBands();
  h_flux_normalization->Reset();

  // Get the flux histo, to be integrated
  static PlotUtils::FluxReweighter* frw = new PlotUtils::FluxReweighter(
      14, CCNuPionIncConsts::kUseNueConstraint, "minervame1D1M1NWeightedAve",
      PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6,
      CCNuPionIncConsts::kNFluxUniverses);

  h_flux_normalization = frw->GetIntegratedFluxReweighted(
      14, h_efficiency_corrected_data, 0., 100.);
  h_flux_normalization->Scale(1.0e-4);
  //  h_flux_normalization->Divide(h_flux_normalization, AaronFluxNormalizer);

  static const double apothemAaron = 850.;
  static const double upstreamAaron = 5991.29;    // ~module 25 plane 1
  static const double downstreamAaron = 8408.91;  // ~module 81 plane 1

  static const double apothemBen = 865.;
  static const double upstreamBen = 5900.;    // ~module 25 plane 1
  static const double downstreamBen = 8430.;  // ~module 81 plane 1

  double n_target_nucleonsAaron =
      PlotUtils::TargetUtils::Get().GetTrackerNNucleons(upstreamAaron,
                                                        downstreamAaron,
                                                        false,  // isMC
                                                        apothemAaron);
  double n_target_nucleonsBen =
      PlotUtils::TargetUtils::Get().GetTrackerNNucleons(upstreamBen,
                                                        downstreamBen,
                                                        false,  // isMC
                                                        apothemBen);
  std::cout << "Number of targets = " << n_target_nucleonsBen << "\n";
  double targetsScale = n_target_nucleonsBen / n_target_nucleonsAaron;
  /*dummyBenXsecDataCorr->Multiply(dummyBenXsecDataCorr, h_flux_normalization);
  dummyBenXsecDataCorr->Scale(targetsScale);
  std::vector<PlotUtils::MnvH1D*> hvec;
  std::vector<PlotUtils::MnvH1D*> hvecinternal;
  std::vector<PlotUtils::MnvH1D*> hvecinternalMixSig;
  std::vector<PlotUtils::MnvH1D*> hvecEffNum;
  std::vector<PlotUtils::MnvH1D*> hvecEffDen;
  std::vector<PlotUtils::MnvH1D*> hvecEff;
  std::vector<PlotUtils::MnvH1D*> hvecEffCorr;
  std::vector<PlotUtils::MnvH1D*> hvecUnfold;
  std::vector<PlotUtils::MnvH1D*> hvecDataSel;
  std::vector<PlotUtils::MnvH1D*> hvecBGSub;

  for (int i = 1; i <= BenXsecData->GetNbinsX(); ++i) {
    BenXsecData->SetBinContent(i, dummyBenXsecData->GetBinContent(i));
    BenXsecDataCorr->SetBinContent(i, dummyBenXsecDataCorr->GetBinContent(i));
    BenEffCorrData->SetBinContent(i, dummyBenEffCorrData->GetBinContent(i));
    BenEffDen->SetBinContent(i, dummyBenEffDen->GetBinContent(i));
    BenDataSel->SetBinContent(i, dummyBenDataSel->GetBinContent(i));
    BenBGSub->SetBinContent(i, dummyBenBGSub->GetBinContent(i));
    BenUnfoldData->SetBinContent(i, dummyBenUnfoldData->GetBinContent(i));
    BenXsecMC->SetBinContent(i, dummyBenXsecMC->GetBinContent(i));
    BenEff->SetBinContent(i, dummyBenEff->GetBinContent(i));
    BenEffNum->SetBinContent(i, dummyBenEffNum->GetBinContent(i));
    AaronNewEffCorr->SetBinContent(
        i, AaronUnfoldData->GetBinContent(i) /
               (AaronEffNum->GetBinContent(i) / AaronEffDen->GetBinContent(i)));
    AaronNewEff->SetBinContent(
        i, AaronEffNum->GetBinContent(i) / AaronEffDen->GetBinContent(i));
    BenXsecAaronSigData->SetBinContent(
        i, dummyBenXsecAaronSigData->GetBinContent(i));
    BenXsecNoTpiWegihtData->SetBinContent(
        i, dummyBenXsecNoTpiWegihtData->GetBinContent(i));
    BenAaronSigEffDen->SetBinContent(i,
                                     dummyBenAaronSigEffDen->GetBinContent(i));
    BenAaronSigEffNum->SetBinContent(i,
                                     dummyBenAaronSigEffNum->GetBinContent(i));
    BenAaronSigEff->SetBinContent(i, dummyBenAaronSigEff->GetBinContent(i));
    BenAaronSigEffCorr->SetBinContent(
        i, dummyBenAaronSigEffCorr->GetBinContent(i));
    BenAaronSigUnfold->SetBinContent(i,
                                     dummyBenAaronSigUnfold->GetBinContent(i));
    BenAaronSigDataSel->SetBinContent(
        i, dummyBenAaronSigDataSel->GetBinContent(i));
    BenAaronSigBGSub->SetBinContent(i, dummyBenAaronSigBGSub->GetBinContent(i));

    BenAaronSigCutsEffDen->SetBinContent(
        i, dummyBenAaronSigCutsEffDen->GetBinContent(i));
    BenAaronSigCutsEffNum->SetBinContent(
        i, dummyBenAaronSigCutsEffNum->GetBinContent(i));
    BenAaronSigCutsEff->SetBinContent(
        i, dummyBenAaronSigCutsEff->GetBinContent(i));
    BenAaronSigCutsEffCorr->SetBinContent(
        i, dummyBenAaronSigCutsEffCorr->GetBinContent(i));
    BenAaronSigCutsUnfold->SetBinContent(
        i, dummyBenAaronSigCutsUnfold->GetBinContent(i));
    BenAaronSigCutsDataSel->SetBinContent(
        i, dummyBenAaronSigCutsDataSel->GetBinContent(i));
    BenXsecUntracked->SetBinContent(i, dummyBenXsecUntracked->GetBinContent(i));
    BenAaronSigCutsBGSub->SetBinContent(
        i, dummyBenAaronSigCutsBGSub->GetBinContent(i));
    BenAaronSigEffDenOldMacro->SetBinContent(
        i, dummyBenAaronSigEffDenOldMacro->GetBinContent(i));
    BenAaronSigEffNumOldMacro->SetBinContent(
        i, dummyBenAaronSigEffNumOldMacro->GetBinContent(i));
    BenAaronSigEffOldMacro->SetBinContent(
        i, dummyBenAaronSigEffOldMacro->GetBinContent(i));
    BenAaronSigEffCorrOldMacro->SetBinContent(
        i, dummyBenAaronSigEffCorrOldMacro->GetBinContent(i));
    BenAaronSigUnfoldOldMacro->SetBinContent(
        i, dummyBenAaronSigUnfoldOldMacro->GetBinContent(i));
    BenAaronSigDataSelOldMacro->SetBinContent(
        i, dummyBenAaronSigDataSelOldMacro->GetBinContent(i));
    BenAaronSigBGSubOldMacro->SetBinContent(
        i, dummyBenAaronSigBGSubOldMacro->GetBinContent(i));
    EffPEMA->SetBinContent(i, BenAaronSigEffNumOldMacro->GetBinContent(i) /
                                  BenAaronSigEffDenOldMacro->GetBinContent(i));
  }

  double AaronDataPOT = h_data_POT->GetBinContent(1);
  double AaronMCPOT = h_MC_POT->GetBinContent(1);
  double DataPOTScale = util.m_data_pot / AaronDataPOT;
  double MCPOTScale = util.m_mc_pot / AaronMCPOT;
  double POTDataScaleASig =
      util.m_data_pot / dataPOTBenAaronSig->GetBinContent(1);
  double POTMCScaleASig = util.m_mc_pot / MCPOTBenAaronSig->GetBinContent(1);
  double POTDataScaleASigCuts =
      util.m_data_pot / dataPOTBenAaronSigCuts->GetBinContent(1);
  double POTMCScaleASigCuts =
      util.m_mc_pot / MCPOTBenAaronSigCuts->GetBinContent(1);
  double POTDataScaleASigOldMacro =
      util.m_data_pot / dataPOTBenAaronSigOldMacro->GetBinContent(1);
  double POTMCScaleASigOldMacro =
      util.m_mc_pot / MCPOTBenAaronSigOldMacro->GetBinContent(1);
  std::cout << "Flux reweighter = " << AaronFluxNormalizer->GetBinContent(1)
            << " Flux scale = " << h_flux_normalization->GetBinContent(1)
            << "\n";
  std::cout << "Number of targets = " << n_target_nucleonsAaron
            << " Target scale = " << targetsScale << "\n";
  std::cout << "Data scale = " << DataPOTScale << " MC Scale = " << DataPOTScale
            << " AaronSigDatascale = " << POTDataScaleASig
            << " AaronSigMCScale = " << POTMCScaleASig << "\n";

  AaronEffDenOtherFile->Scale(AaronMCPOT / AaronDataPOT);
  AaronEffNumOtherFile->Scale(AaronMCPOT / AaronDataPOT);
  AaronEffDenOtherFile->Scale(MCPOTScale);
  AaronEffNumOtherFile->Scale(MCPOTScale);
  BenAaronSigEffDen->Scale(POTMCScaleASig);
  BenAaronSigEffNum->Scale(POTMCScaleASig);
  BenAaronSigCutsEffDen->Scale(POTMCScaleASigCuts);
  BenAaronSigCutsEffNum->Scale(POTMCScaleASigCuts);
  BenAaronSigEffDenOldMacro->Scale(POTMCScaleASigOldMacro);
  BenAaronSigEffNumOldMacro->Scale(POTMCScaleASigOldMacro);

  //  AaronEffOtherFile->Scale(AaronMCPOT/AaronDataPOT);
  AaronUnfoldData->Scale(DataPOTScale);
  AaronCorrectedData->Scale(DataPOTScale);
  AaronDataSel->Scale(DataPOTScale);
  AaronBGSub->Scale(DataPOTScale);

  BenAaronSigUnfold->Scale(POTDataScaleASig);
  BenAaronSigDataSel->Scale(POTDataScaleASig);
  BenAaronSigBGSub->Scale(POTDataScaleASig);
  BenAaronSigEffCorr->Scale(POTDataScaleASig);

  BenAaronSigCutsUnfold->Scale(POTDataScaleASigCuts);
  BenAaronSigCutsDataSel->Scale(POTDataScaleASigCuts);
  BenAaronSigCutsBGSub->Scale(POTDataScaleASigCuts);
  BenAaronSigCutsEffCorr->Scale(POTDataScaleASigCuts);

  BenAaronSigUnfoldOldMacro->Scale(POTDataScaleASigOldMacro);
  BenAaronSigDataSelOldMacro->Scale(POTDataScaleASigOldMacro);
  BenAaronSigBGSubOldMacro->Scale(POTDataScaleASigOldMacro);
  BenAaronSigEffCorrOldMacro->Scale(POTDataScaleASigOldMacro);

  BenEffDen->Scale(MCPOTScale);
  hvec.push_back(BenXsecData);
  //  hvec.push_back(BenXsecDataCorr);
  hvec.push_back(BenXsecAaronSigData);
  hvec.push_back(BenXsecNoTpiWegihtData);
  hvec.push_back(BenXsecUntracked);
  hvecEffNum.push_back(BenEffNum);
  hvecEffNum.push_back(BenAaronSigEffNum);
  //  hvecEffNum.push_back(BenAaronSigCutsEffNum);
  hvecEffNum.push_back(BenAaronSigEffNumOldMacro);
  hvecEffDen.push_back(BenEffDen);
  hvecEffDen.push_back(BenAaronSigEffDen);
  //  hvecEffDen.push_back(BenAaronSigCutsEffDen);
  hvecEffDen.push_back(BenAaronSigEffDenOldMacro);
  hvecinternal.push_back(BenXsecData);
  hvecinternal.push_back(BenXsecNoTpiWegihtData);
  hvecinternal.push_back(BenXsecUntracked);

  hvecinternalMixSig.push_back(BenXsecNoTpiWegihtData);
  hvecinternalMixSig.push_back(BenXsecUntracked);

  hvecEff.push_back(BenEff);
  hvecEff.push_back(BenAaronSigEff);
  hvecEff.push_back(EffPEMA);
  //  hvecEff.push_back(BenAaronSigCutsEff);
  hvecEff.push_back(BenAaronSigEffOldMacro);

  hvecEffCorr.push_back(BenEffCorrData);
  hvecEffCorr.push_back(BenAaronSigEffCorr);
  //  hvecEffCorr.push_back(BenAaronSigCutsEffCorr);
  hvecEffCorr.push_back(BenAaronSigEffCorrOldMacro);

  hvecUnfold.push_back(BenUnfoldData);
  hvecUnfold.push_back(BenAaronSigUnfold);
  //  hvecUnfold.push_back(BenAaronSigCutsUnfold);
  hvecUnfold.push_back(BenAaronSigUnfoldOldMacro);

  hvecDataSel.push_back(BenDataSel);
  hvecDataSel.push_back(BenAaronSigDataSel);
  //  hvecDataSel.push_back(BenAaronSigCutsDataSel);
  hvecDataSel.push_back(BenAaronSigDataSelOldMacro);

  hvecBGSub.push_back(BenBGSub);
  hvecBGSub.push_back(BenAaronSigBGSub);
  //  hvecBGSub.push_back(BenAaronSigCutsBGSub);
  hvecBGSub.push_back(BenAaronSigBGSubOldMacro);

  AaronRevUnfoldData->Multiply(AaronRevUnfoldData, AaronFluxNormalizer);
  AaronRevUnfoldData->Multiply(AaronRevUnfoldData, AaronNewEff);
  AaronRevUnfoldData->Scale(h_data_POT->GetBinContent(1));
  AaronRevUnfoldData->Scale(n_target_nucleonsAaron);

  AaronNoNormXsecData->Multiply(AaronNoNormXsecData, AaronFluxNormalizer);
  AaronNoNormXsecData->Scale(h_data_POT->GetBinContent(1));
  AaronNoNormXsecData->Scale(n_target_nucleonsAaron);

  PlotRatio(BenXsecGENIE, AaronXsecGENIE, "Q2", 1., "BenSigvsAaronsigGENIE",
            false, "BenSigvsAaronsig", "Q^{2}");
  PlotRatio(BenXsecGENIE, AaronXsecGENIE, "Q2", 1., "BenSigvsAaronsigGENIE",
            false, "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatio(BenXsecMC, AaronXsecMC, "Q2", 1., "BenSigvsAaronXSecMC", false,
            "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenXsecData, BenXsecAaronSigData, "Q2", 1.,
            "BenMacroAaronSigXSecData", false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenEffDen, AaronEffDen, "Q2", MCPOTScale, "BenSigvsAaronEffDen",
            false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenEffCorrData, AaronCorrectedData, "Q2", DataPOTScale,
            "BenSigvsAaronEffCorrection", false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenEffCorrData, AaronNewEffCorr, "Q2", DataPOTScale,
            "BenSigvsAaronNewEffCorrection", false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenEffCorrData, AaronNoNormXsecData, "Q2", DataPOTScale,
            "BenSigvsAaronNoNormEffCorr", false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenEff, AaronNewEff, "Q2", 1., "BenSigvsAaronNewEff", false,
            "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(AaronNewEff, AaronEff, "Q2", 1., "AaronEffCheck", false,
            "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenUnfoldData, AaronUnfoldData, "Q2", DataPOTScale,
            "BenSigvsAaronUnfoldData", false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenUnfoldData, AaronRevUnfoldData, "Q2", DataPOTScale,
            "BenSigvsAaronRevUnfoldData", false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(AaronRevUnfoldData, AaronUnfoldData, "Q2", 1.,
            "AaronUnfoldDataCheck", false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(AaronNewEffCorr, AaronCorrectedData, "Q2", 1., "AaronEffCorrCheck",
            false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenEff, AaronEff, "Q2", 1., "BenSigvsAaronEff", false,
            "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(AaronEffOtherFile, AaronEff, "Q2", 1., "AaronEffDiffFiles", false,
            "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(AaronNewEff, AaronEffOtherFile, "Q2", 1., "AaronNewEffvsOtherFile",
            false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(BenEffNum, AaronEffNum, "Q2", MCPOTScale, "BenSigvsAaronEffNum",
            false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(AaronEffNumOtherFile, AaronEffNum, "Q2", 1., "AaronEffNumDiffFiles",
            false, "BenSig/Aaronsig", "Q^{2}");
  PlotRatio(AaronEffDenOtherFile, AaronEffDen, "Q2", 1., "AaronEffDenDiffFiles",
            false, "BenSig/Aaronsig", "Q^{2}");
  // PlotRatio(AaronEffNum, AaronEffDen, "Q2", 1.,"AaronEfficiency", true);
  PlotRatioVec(hvec, AaronXsecData, "Q2", 1., "CrossSectionsData", false,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecEffNum, AaronEffNumOtherFile, "Q2", 1., "EffNum", false,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecEffDen, AaronEffDenOtherFile, "Q2", 1., "EffDen", true,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecEff, AaronEffOtherFile, "Q2", 1., "Eff", false,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecEffCorr, AaronCorrectedData, "Q2", 1., "EffCorr", false,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecUnfold, AaronUnfoldData, "Q2", 1., "Unfold", false,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecDataSel, AaronDataSel, "Q2", 1., "DataSelection", false,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecBGSub, AaronBGSub, "Q2", 1., "BGSub", false,
               "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecinternal, BenXsecAaronSigData, "Q2", 1.,
               "InternalCrossSectionsData", false, "BenSig/Aaronsig",
               "Q^{2} (GeV^{2})");
  PlotRatioVec(hvecinternalMixSig, BenXsecData, "Q2", 1.,
               "InternalCrossSectionsDataTrackerSig", false, "",
               "Q^{2} (GeV^{2})");

  PlotRatio(BenEffNum, BenAaronSigEffNum, "Q2", 1, "OurMasterBranchEffNum",
            false, "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatio(BenEffDen, BenAaronSigEffDen, "Q2", 1, "OurMasterBranchEffDen",
            false, "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatio(BenEff, BenAaronSigEff, "Q2", 1, "OurMasterBranchEff", false,
            "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatio(BenXsecData, BenXsecAaronSigData, "Q2", 1, "OurMasterBranchXSec",
            false, "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatio(BenEffCorrData, BenAaronSigEffCorr, "Q2", 1,
            "OurMasterBranchEffCorr", false, "BenSig/Aaronsig",
            "Q^{2} (GeV^{2})");
  PlotRatio(BenUnfoldData, BenAaronSigUnfold, "Q2", 1, "OurMasterBranchUnfold",
            false, "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatio(BenDataSel, BenAaronSigDataSel, "Q2", 1, "OurMasterBranchDataSel",
            false, "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
  PlotRatio(BenBGSub, BenAaronSigBGSub, "Q2", 1, "OurMasterBranchBGSub", false,
            "BenSig/Aaronsig", "Q^{2} (GeV^{2})");
*/
  for (auto var : variables) {
    if (var->Name() == sidebands::kFitVarString) {
      var->m_hists.m_stacked_wsideband = StackedHistogram<WSidebandType>(
          var->m_hists.m_label, var->m_hists.m_xlabel,
          var->m_hists.m_bins_array, kNWSidebandTypes,
          sidebands::kWSideband_ColorScheme);
    }
    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    var->LoadDataHistsFromFile(fin);

    // if(var->Name() == "ptmu")
    //  PrintUniverseContent(var->m_hists.m_cross_section);
  }

  {  // remove unwanted variables
    ContainerEraser::erase_if(variables, [](Variable* v) {
      return v->Name() == "tpi_mbr";  // || v->Name() == "wexp_fit";
    });

    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "tpi" || v->Name() == "enu"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "thetapi_deg" || v->Name() == "thetamu_deg"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "q2" || v->Name() == "wexp"; });
    // ContainerEraser::erase_if(variables, [](Variable* v) {
    //    return v->Name() == "ptmu" || v->Name() == "pzmu"; });
  }

  // PLOT Event Selection, BGs (error)
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    bool do_cov_area_norm = false;
    for (auto var : variables) {
      if (var->Name() == "wexp_fit") continue;
      std::cout << var->Name() << "\n";
      do_cov_area_norm = false;
      Plotter plot_info(var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                        do_cov_area_norm, include_stat,
                        util.m_signal_definition);

      bool do_bin_width_norm = true, do_log_scale = false, do_bg = true;
      bool do_tuned_bg = false;

      // selection and BG error before tuning
      if (!var->m_is_true) {
        PlotVar_Selection(plot_info, -1., do_log_scale, do_bg, do_tuned_bg,
                          do_bin_width_norm);
        if (plot_errors) PlotVar_ErrorSummary(plot_info);
        if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

        // selection and BG error after tuning
        do_tuned_bg = true;
        PlotVar_Selection(plot_info, -1., do_log_scale, do_bg, do_tuned_bg,
                          do_bin_width_norm);
        if (plot_errors) PlotVar_ErrorSummary(plot_info);
        if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);
      }
    }
  }
  if (true) {
    for (auto var : variables) {
      if (var->Name() == "wexp_fit") continue;
      PlotUtils::MnvH1D* h_bkdtrackedtpi = nullptr;
      //        (PlotUtils::MnvH1D*)finaux.Get(Form("selection_data_tracked_%s",
      //        var->Name().c_str()));
      PlotUtils::MnvH1D* h_bkdtracklesstpi = nullptr;
      //        (PlotUtils::MnvH1D*)finaux.Get(Form("selection_data_untracked_%s",
      //        var->Name().c_str()));
      PlotUtils::MnvH1D* h_bkdmixtpi = nullptr;
      //        (PlotUtils::MnvH1D*)finaux.Get(Form("selection_data_mixed_%s",
      //        var->Name().c_str()));
      PlotUtils::MnvH1D* h_data = nullptr;
      //        (PlotUtils::MnvH1D*)fin.Get(Form("selection_data_%s",
      //        var->Name().c_str()));
      PlotUtils::MnvH1D* h_mc_POT = (PlotUtils::MnvH1D*)finaux.Get("mc_pot");

      PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
      TCanvas cE("c1", "c1");

      if (var->Name() == "q2") {
        h_bkdtrackedtpi = RebinQ2Plot(*(PlotUtils::MnvH1D*)fin.Get(
            Form("selection_data_tracked_%s", var->Name().c_str())));
        h_bkdtracklesstpi = RebinQ2Plot(*(PlotUtils::MnvH1D*)fin.Get(
            Form("selection_data_untracked_%s", var->Name().c_str())));
        h_bkdmixtpi = RebinQ2Plot(*(PlotUtils::MnvH1D*)fin.Get(
            Form("selection_data_mixed_%s", var->Name().c_str())));
        h_data = RebinQ2Plot(*(PlotUtils::MnvH1D*)fin.Get(
            Form("selection_data_%s", var->Name().c_str())));
        cE.SetLogx();
      } else if (var->Name() == "enu" || var->Name() == "pmu" ||
                 var->Name() == "ptmu" || var->Name() == "pzmu") {
        // hist = RebinningtoGeV(*doomyhist, var->Name());
        h_bkdtrackedtpi = RebinningtoGeV(
            *(PlotUtils::MnvH1D*)finaux.Get(
                Form("selection_data_tracked_%s", var->Name().c_str())),
            var->Name());
        h_bkdtracklesstpi = RebinningtoGeV(
            *(PlotUtils::MnvH1D*)finaux.Get(
                Form("selection_data_untracked_%s", var->Name().c_str())),
            var->Name());
        h_bkdmixtpi =
            RebinningtoGeV(*(PlotUtils::MnvH1D*)finaux.Get(Form(
                               "selection_data_mixed_%s", var->Name().c_str())),
                           var->Name());
        h_data = RebinningtoGeV(*(PlotUtils::MnvH1D*)fin.Get(Form(
                                    "selection_data_%s", var->Name().c_str())),
                                var->Name());
      } else {
        // hist = (PlotUtils::MnvH1D*)doomyhist->Clone("hist");
        h_bkdtrackedtpi = (PlotUtils::MnvH1D*)finaux.Get(
            Form("selection_data_tracked_%s", var->Name().c_str()));
        h_bkdtracklesstpi = (PlotUtils::MnvH1D*)finaux.Get(
            Form("selection_data_untracked_%s", var->Name().c_str()));
        h_bkdmixtpi = (PlotUtils::MnvH1D*)finaux.Get(
            Form("selection_data_mixed_%s", var->Name().c_str()));
        h_data = (PlotUtils::MnvH1D*)fin.Get(
            Form("selection_data_%s", var->Name().c_str()));
      }

      SetHistColorScheme(h_bkdtrackedtpi, 1, 5);
      SetHistColorScheme(h_bkdtracklesstpi, 2, 5);
      SetHistColorScheme(h_bkdmixtpi, 3, 5);

      double scl = util.m_data_pot / h_mc_POT->GetBinContent(1);
      TObjArray* stack = new TObjArray();
      h_bkdtrackedtpi->SetTitle("Tracked");
      h_bkdtracklesstpi->SetTitle("Untracked");
      h_bkdmixtpi->SetTitle("Overlap");
      h_bkdtrackedtpi->GetYaxis()->SetTitle(
          Form("Events/%s", var->m_units.c_str()));
      h_data->GetYaxis()->SetTitle(Form("Events/%s", var->m_units.c_str()));
      h_bkdtrackedtpi->Scale(1, "width");
      h_bkdtracklesstpi->Scale(1, "width");
      h_bkdmixtpi->Scale(1, "width");
      h_data->Scale(1., "width");

      double ymax = h_data->GetBinContent(h_data->GetMaximumBin()) * 1.1;

      stack->Add(h_bkdtrackedtpi);
      stack->Add(h_bkdmixtpi);
      stack->Add(h_bkdtracklesstpi);
      mnvPlotter.axis_maximum = ymax;
      mnvPlotter.axis_title_offset_y = 1.;
      // mnvPlotter.DrawStackedMC(stack, 1.0, "TR", -1, -1, 3000, "T_{#pi}
      // (MeV)");
      std::string lbel = var->m_hists.m_xlabel + " " + var->m_units;
      mnvPlotter.DrawDataStackedMC(h_data, stack, 1, "TR", "Data", -1, -1, 3002,
                                   lbel.c_str());
      // mnvPlotter.AddHistoTitle("T_{#pi} Breakdown", 0.05);
      mnvPlotter.WriteNorm(Form("Data POT = %.2E", util.m_data_pot), 0.3, 0.88,
                           0.03);
      mnvPlotter.AddPlotLabel("Stat-error-only", 0.74, 0.73, 0.03, 9);
      // std::string plotname = "Stacked_" + var->Name() +
      // "_thetapisample_tpiweight";
      std::string plotname = "Stacked_" + var->Name() + "_tpiweight";

      mnvPlotter.MultiPrint(&cE, plotname, "pdf");
    }
  }
  // PLOT Efficiency & Migration
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    for (auto var : variables) {
      if (var->Name() == "wexp_fit") continue;
      // var->LoadMCHistsFromFile(fin, util.m_error_bands);

      const Plotter plot_info(var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                              do_cov_area_norm, include_stat,
                              util.m_signal_definition);

      // Efficiency
      if (var->m_is_true) {
        PlotUtils::MnvH1D* eff =
            (PlotUtils::MnvH1D*)var->m_hists.m_effnum.hist->Clone(uniq());
        PlotUtils::MnvH1D* effnum =
            (PlotUtils::MnvH1D*)var->m_hists.m_effnum.hist->Clone(uniq());
        PlotUtils::MnvH1D* effden =
            (PlotUtils::MnvH1D*)var->m_hists.m_effden.hist->Clone(uniq());
        eff->Divide(effnum, effden);
        var->m_hists.m_efficiency = (PlotUtils::MnvH1D*)eff->Clone(uniq());

        double ymax = -1.;
        bool do_log_scale = false;
        bool do_bg = true;
        bool do_tuned_bg = true;
        PlotMC(eff, plot_info, Form("Efficiency_%s", var->Name().c_str()), -1.,
               "Efficiency");
        if (plot_errors) PlotEfficiency_ErrorSummary(plot_info);
      }

      // Migration
      if (!var->m_is_true) {
        PlotUtils::MnvH2D* mig =
            (PlotUtils::MnvH2D*)var->m_hists.m_migration.hist->Clone(uniq());
        PlotMigration_AbsoluteBins(mig, var->Name());
        PlotMigration_VariableBins(mig, var->Name());
      }
    }
  }

  // PLOT Background Subtraction
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    for (auto var : variables) {
      if (var->Name() == "wexp_fit") continue;
      if (var->m_is_true) continue;

      Plotter plot_info(var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                        do_cov_area_norm, include_stat,
                        util.m_signal_definition);

      double ymax = -1.;
      bool do_log_scale = false;
      bool do_bg = true;
      bool do_tuned_bg = true;
      bool do_bin_width_norm = true;
      if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

      do_tuned_bg = false;
      if (plot_errors) PlotBG_ErrorSummary(plot_info, do_tuned_bg);

      Plot_BGSub(plot_info, ".", ymax, do_log_scale, do_bin_width_norm);

      if (plot_errors) PlotBGSub_ErrorSummary(plot_info);
    }
  }

  // PLOT W Sideband Fit
  if (false) {
    const bool do_frac_unc = true;
    const bool do_cov_area_norm = false;
    const bool include_stat = true;

    Plotter plot_info(util.m_mc_pot, util.m_data_pot, do_frac_unc,
                      do_cov_area_norm, include_stat, util.m_signal_definition);

    PlotUtils::MnvH1D* loW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("loW_fit_wgt");
    PlotUtils::MnvH1D* midW_fit_wgt =
        (PlotUtils::MnvH1D*)fin.Get("midW_fit_wgt");
    PlotUtils::MnvH1D* hiW_fit_wgt = (PlotUtils::MnvH1D*)fin.Get("hiW_fit_wgt");

    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, loW_fit_wgt, "WSidebandFit_loW");

    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, midW_fit_wgt,
                                    "WSidebandFit_midW");
    if (plot_errors)
      PlotWSidebandFit_ErrorSummary(plot_info, hiW_fit_wgt, "WSidebandFit_hiW");

    // Wexp distribution, stacked, all cuts except for W, pre-fit
    std::string tag;
    double ymax = 440000;
    double do_BWN = true;
    double do_postfit = false;
    Variable* var = GetVar(variables, sidebands::kFitVarString);
    PlotWSidebandStacked(
        var, var->m_hists.m_wsideband_data, loW_fit_wgt, midW_fit_wgt,
        hiW_fit_wgt, var->GetStackArray(static_cast<WSidebandType>(0)),
        util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag, ymax,
        do_BWN, do_postfit);

    do_postfit = true;
    PlotWSidebandStacked(
        var, var->m_hists.m_wsideband_data, loW_fit_wgt, midW_fit_wgt,
        hiW_fit_wgt, var->GetStackArray(static_cast<WSidebandType>(0)),
        util.m_data_pot, util.m_mc_pot, util.m_signal_definition, tag, ymax,
        do_BWN, do_postfit);
    /*for (auto var : variables) {
      std::string name = var->Name();
      if (var->Name() == "wexp_fit") continue;
      if (var->m_is_true) continue;
      PlotUtils::MnvH1D* h_sig = (PlotUtils::MnvH1D*)fin.Get(
          Form("wsidebandfit_sig_%s", name.c_str()));
      PlotUtils::MnvH1D* h_loW = (PlotUtils::MnvH1D*)fin.Get(
          Form("wsidebandfit_loW_%s", name.c_str()));
      PlotUtils::MnvH1D* h_midW = (PlotUtils::MnvH1D*)fin.Get(
          Form("wsidebandfit_midW_%s", name.c_str()));
      PlotUtils::MnvH1D* h_hiW = (PlotUtils::MnvH1D*)fin.Get(
          Form("wsidebandfit_hiW_%s", name.c_str()));

      PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
      TCanvas cE("c1", "c1");
      TObjArray* stack = new TObjArray();
      h_sig->SetTitle("Signal");
      h_loW->SetTitle("Low W");
      h_midW->SetTitle("Mid W");
      h_hiW->SetTitle("High W");

      h_sig->GetYaxis()->SetTitle(Form("Events/%s", var->m_units.c_str()));
      h_sig->Scale(1., "width");
      h_loW->Scale(1., "width");
      h_midW->Scale(1., "width");
      h_hiW->Scale(1., "width");
      //     if(name == "q2")  cE.SetLogx();

      stack->Add(h_sig);
      stack->Add(h_loW);
      stack->Add(h_midW);
      stack->Add(h_hiW);
      mnvPlotter.DrawStackedMC(
          stack, 1.0, "TR", 2, 1, 3001,
          Form("%s %s", var->m_hists.m_xlabel.c_str(), var->m_units.c_str()));
      mnvPlotter.AddHistoTitle("SidebandRegion", 0.05);

      std::string plotname = "SidebandRegion_" + name;
      mnvPlotter.MultiPrint(&cE, plotname, "png");
    }*/
    // TODO plot pre/postfit
    ymax = -1;
    for (auto var : variables) {
      // if (var->Name() != "wexp_fit") continue;
      if (var->Name() == "wexp") continue;
      if (var->m_is_true) continue;
      tag = "SidebandRegion";
      bool do_prefit = true;
      bool do_bin_width_norm = true;
      CVUniverse* universe = util.m_error_bands.at("cv").at(0);
      PlotFittedW(var, *universe, loW_fit_wgt, midW_fit_wgt, hiW_fit_wgt,
                  util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                  do_prefit, tag, ymax, do_bin_width_norm);
      do_prefit = false;
      PlotFittedW(var, *universe, loW_fit_wgt, midW_fit_wgt, hiW_fit_wgt,
                  util.m_data_pot, util.m_mc_pot, util.m_signal_definition,
                  do_prefit, tag, ymax, do_bin_width_norm);
    }
  }

  // PLOT unfolded
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    const double ymax = -1.;
    const bool do_log_scale = false;
    const bool do_bin_width_norm = true;
    for (auto reco_var : variables) {
      if (reco_var->Name() == "wexp_fit") continue;
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                        do_cov_area_norm, include_stat,
                        util.m_signal_definition);

      Plot_Unfolded(plot_info, reco_var->m_hists.m_unfolded,
                    true_var->m_hists.m_effnum.hist);
      if (plot_errors) PlotUnfolded_ErrorSummary(plot_info);
    }
  }

  // PLOT cross section
  if (false) {
    const bool do_frac_unc = true;
    const bool include_stat = true;
    const bool do_cov_area_norm = false;
    const double ymax = -1.;
    const bool do_log_scale = false;
    bool do_bin_width_norm = true;
    PlotUtils::MnvH1D* h_dummyflux =
        (PlotUtils::MnvH1D*)fin.Get("cross_section_enu");
    PlotUtils::MnvH1D* h_flux =
        (PlotUtils::MnvH1D*)h_dummyflux->Clone("enu_clone");
    h_flux->ClearAllErrorBands();
    h_flux->Reset();

    // Get the flux histo, to be integrated
    static PlotUtils::FluxReweighter* frw = new PlotUtils::FluxReweighter(
        14, CCNuPionIncConsts::kUseNueConstraint, "minervame1D1M1NWeightedAve",
        PlotUtils::FluxReweighter::gen2thin,
        PlotUtils::FluxReweighter::g4numiv6,
        CCNuPionIncConsts::kNFluxUniverses);
    PlotUtils::MnvH1D* flux_normalization =
        (PlotUtils::MnvH1D*)h_dummyflux->Clone("h_flux_normalization");
    flux_normalization->ClearAllErrorBands();
    flux_normalization->Reset();
    flux_normalization =
        frw->GetIntegratedFluxReweighted(14, h_dummyflux, 0., 100.);
    flux_normalization->Scale(1.0e-4);
    h_dummyflux->AddMissingErrorBandsAndFillWithCV(*flux_normalization);
    std::vector<std::string> fluxnor_bands =
        flux_normalization->GetVertErrorBandNames();
    PlotUtils::MnvH1D* GeVflux = nullptr;
    GeVflux = RebinningtoGeV(*h_flux, "flux");
    PlotUtils::MnvH1D* flux =
        UndoBWN(frw->GetRebinnedFluxReweighted(14, GeVflux));
    double lowbin = 0, hibin = 0;
    flux->Scale(1.0e-4);

    for (int i = 1; i <= h_flux->GetNbinsX(); i++) {
      lowbin = h_flux->GetBinLowEdge(i) / 1000;
      hibin = h_flux->GetBinLowEdge(i + 1) / 1000;
      std::cout << "Bin = " << i << " Low binedge = " << lowbin
                << " High binedge = " << hibin << " " << flux->GetBinContent(i)
                << "\n";
      h_flux->SetBinContent(i, flux->GetBinContent(i));
    }

    //  h_dummyflux->AddMissingErrorBandsAndFillWithCV(*h_flux);
    TH1* h_flux_aux = (TH1*)h_flux->Clone("TH1Flux");

    /* for (int i = 1; i <= h_flux_aux->GetNbinsX(); i++){
       lowbin = h_flux->GetBinLowEdge(i)/1000;
       hibin = h_flux->GetBinLowEdge(i+1)/1000;
       std::cout << "Bin = " << i << " Low binedge = " << lowbin <<
                    " High binedge = " << hibin << " " <<
                     flux->GetBinContent(i) << "\n";
       h_flux->SetBinContent(i, flux->GetBinContent(i));
     }*/

    for (auto reco_var : variables) {
      do_bin_width_norm = true;
      if (reco_var->Name() == "wexp_fit") continue;
      if (reco_var->m_is_true) continue;
      Variable* true_var =
          GetVar(variables, reco_var->Name() + std::string("_true"));

      Plotter plot_info(reco_var, util.m_mc_pot, util.m_data_pot, do_frac_unc,
                        do_cov_area_norm, include_stat,
                        util.m_signal_definition);

      PlotUtils::MnvH1D* m_mc_cross_section = (PlotUtils::MnvH1D*)fin.Get(
          Form("mc_cross_section_%s", reco_var->Name().c_str()));

      std::vector<std::string> x_bands =
          reco_var->m_hists.m_cross_section->GetVertErrorBandNames();
      std::vector<std::string> flux_bands = h_flux->GetVertErrorBandNames();
      if (reco_var->Name() == "enu") {
        for (auto a : flux_bands) {
          for (int i = 1; i <= h_flux->GetNbinsX(); i++) {
            h_flux->GetVertErrorBand(a)
                ->GetErrorBand(false, false)
                .SetBinContent(i, h_flux->GetBinContent(i));
          }
        }
        for (auto s : x_bands)
          for (auto a : flux_bands) {
            TH1D* hErr =
                dynamic_cast<TH1D*>(h_flux->GetVertErrorBand(a)
                                        ->GetErrorBand(false, false)
                                        .Clone(Form("Flux%s", a.c_str())));
          }
      }
      // std::cout << reco_var->Name() << "\n";

      TObjArray* stack = new TObjArray();
      static PlotUtils::FluxReweighter* frw1 = new PlotUtils::FluxReweighter(
          14, CCNuPionIncConsts::kUseNueConstraint,
          "minervame1D1M1NWeightedAve", PlotUtils::FluxReweighter::gen2thin,
          PlotUtils::FluxReweighter::g4numiv6,
          CCNuPionIncConsts::kNFluxUniverses);
      PlotUtils::MnvH1D* h_flux_norm =
          reco_var->m_hists.m_cross_section->Clone("FLUX");
      PlotUtils::MnvH1D* h_mc_POT = (PlotUtils::MnvH1D*)finaux.Get("mc_pot");
      double mcpotaux = h_mc_POT->GetBinContent(1);
      h_flux_norm->ClearAllErrorBands();
      h_flux_norm->Reset();
      h_flux_norm =
          frw1->GetIntegratedFluxReweighted(14, h_flux_norm, 0., 100.);
      h_flux_norm->Scale(1.0e-4);

      // double mc_scale = 1.0 / (n_target_nucleonsBen * util.m_mc_pot);
      double mc_scale = 1.0 / (n_target_nucleonsBen * mcpotaux);
      /*for (int i = 0; i < 4; i++) {
        TObject* obj = finaux.Get(
            Form("%s_true_Int_%d", reco_var->Name().c_str(), i));
        TH1* h = dynamic_cast<TH1*>(obj);
        h->Divide(h,h_flux_norm);
        h->Scale(mc_scale);
        stack->Add(h);
      }*/

      if (reco_var->Name() == "enu") {
        do_bin_width_norm = false;
        m_mc_cross_section->Multiply(m_mc_cross_section, flux_normalization);
        reco_var->m_hists.m_cross_section->Multiply(
            reco_var->m_hists.m_cross_section, flux_normalization);
        m_mc_cross_section->DivideSingle(m_mc_cross_section, h_flux_aux);
        reco_var->m_hists.m_cross_section->DivideSingle(
            reco_var->m_hists.m_cross_section, h_flux_aux);
      }

      //      PlotStackedXSec(reco_var, reco_var->m_hists.m_cross_section,
      //      *stack,
      //                    util.m_data_pot, mcpotaux,
      //                    util.m_signal_definition, "Int");

      Plot_CrossSection(plot_info, reco_var->m_hists.m_cross_section,
                        m_mc_cross_section, ".", -1, false, do_bin_width_norm);
      if (plot_errors)
        PlotCrossSection_ErrorSummary(
            plot_info);  // Adds chi2 label and prints out assumed binning.

      // Various error matrices
      // TMatrixD
      // scaled_mtx(reco_var->m_hists.m_cross_section->GetTotalErrorMatrix());
      // TMatrixD scaled_mtx(m_mc_cross_section->GetStatErrorMatrix());

      // TMatrixD
      // data_stat_err_mtx(reco_var->m_hists.m_cross_section->GetStatErrorMatrix());
      // TMatrixD
      // data_tot_err_mtx(reco_var->m_hists.m_cross_section->GetTotalErrorMatrix());
      // if (plot_errors) PlotMatrix(data_stat_err_mtx, reco_var->Name(),
      // "Data_Stat_Err_Mtx"); if (plot_errors) PlotMatrix(data_tot_err_mtx,
      // reco_var->Name(), "Data_Tot_Err_Mtx");

      // TMatrixD mc_stat_err_mtx(m_mc_cross_section->GetStatErrorMatrix());
      // TMatrixD mc_tot_err_mtx(m_mc_cross_section->GetTotalErrorMatrix());
      // if (plot_errors) PlotMatrix(mc_stat_err_mtx,   reco_var->Name(),
      // "MC_Stat_Err_Mtx"); if (plot_errors) PlotMatrix(mc_tot_err_mtx,
      // reco_var->Name(), "MC_Tot_Err_Mtx");

      // TMatrixD
      // data_sys_err_mtx(reco_var->m_hists.m_cross_section->GetSysErrorMatrix());
      // TMatrixD
      // data_sys_corr_mtx(reco_var->m_hists.m_cross_section->GetSysCorrelationMatrix
      // ()); if (plot_errors) PlotMatrix(data_sys_err_mtx,  reco_var->Name(),
      // "Data_Sys_Err_Mtx"); if (plot_errors) PlotMatrix(data_sys_cor_mtx,
      // reco_var->Name(), "Data_Sys_Cor_Mtx");

      // TMatrixD mc_sys_err_mtx (m_mc_cross_section->GetSysErrorMatrix());
      // TMatrixD mc_sys_cor_mtx(m_mc_cross_section->GetSysCorrelationMatrix());
      // if (plot_errors) PlotMatrix(mc_sys_err_mtx,    reco_var->Name(),
      // "MC_Sys_Err_Mtx"); if (plot_errors) PlotMatrix(mc_sys_cor_mtx,
      // reco_var->Name(), "MC_Sys_Cor_Mtx");

      // if (plot_errors) PlotCovarianceMatrix(scaled_mtx, reco_var->Name(),
      // "Data_StatErrMtx"); for(int c = 0; c < scaled_mtx.GetNcols(); ++c ) {
      //  for(int r = 0; r < scaled_mtx.GetNrows(); ++r ) {
      //    std::cout << scaled_mtx(r, c) << "  ";
      //  }
      //  std::cout << "\n";
      //}

      // for(int c = 0; c < scaled_mtx2.GetNcols(); ++c ) {
      //  for(int r = 0; r < scaled_mtx2.GetNrows(); ++r ) {
      //    std::cout << scaled_mtx2(r, c) << "  ";
      //  }
      //  std::cout << "\n";
      //}

      // std::cout << "data (" << scaled_mtx.GetNrows() << "," <<
      // scaled_mtx.GetNcols() << ")\n"; std::cout << "mc (" <<
      // scaled_mtx2.GetNrows() << "," << scaled_mtx2.GetNcols() << ")\n";

      // PrintChi2Info(plot_info, reco_var->m_hists.m_cross_section,
      // m_mc_cross_section); // this one works
    }
  }

  //============================================================================
}

#endif  // plotCrossSectionFromFile_C
