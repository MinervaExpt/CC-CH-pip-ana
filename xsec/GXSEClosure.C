#ifndef GXSEClosure_C
#define GXSEClosure_C

#include "PlotUtils/MnvH1D.h"
#include "TFile.h"
#include "includes/SignalDefinition.h"
#include "includes/Variable.h"
#include "makeCrossSectionMCInputs.C" // GetAnalysisVariables
#include "includes/CCPiMacroUtil.h"
#include "plotting_functions.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/TargetUtils.h"


//==============================================================================
// Main
// TODO function-ify the real xsec script and use those functions.
// Instead of copy-paste.
//==============================================================================
void GXSEClosure(int signal_definition_int = 0) {
  // In and outfiles
    //TFile fin("rootfiles/MCXSecInputs_20190903.root", "READ");
    TFile fin("rootfiles/MCXSecInputs_20190904_ME1A.root", "READ");
    cout << "Reading input from " << fin.GetName() << endl;

  // Set up macro utility object -- which does the systematics for us
    const char* plist = "ALL";
    bool do_data = false, do_mc = false, do_truth = false;
    bool do_systematics = true, do_grid = false;
    CCPiMacroUtil util(signal_definition_int, plist, do_data, do_mc, do_truth,
                       do_systematics, do_grid);

  // Set POT
  PlotUtils::MnvH1D* h_mc_pot=(PlotUtils::MnvH1D*)fin.Get("mc_pot");
  float mc_pot = h_mc_pot->GetBinContent(1);
  util.m_mc_pot = mc_pot;

  // Variables and their histograms
  const bool do_truth_vars = true;
  std::vector<Variable*> variables = GetAnalysisVariables(util.m_signal_definition, 
                                                          do_truth_vars);

  for (auto var : variables)
    var->LoadMCHistsFromFile(fin, util.m_error_bands);

  for (auto var : variables) {
    var->LoadMCHistsFromFile(fin, util.m_error_bands);
    if (var->m_is_true) continue;
    if (var->Name() == std::string("tpi_mbr"))  continue;
    if (var->Name() == sidebands::kFitVarString) continue; 

    if(var->Name() != "pmu") continue;

    const char* name = var->Name().c_str();

    // We'll be needing the true version of this variable
    Variable* true_var = GetVar(variables, var->Name() + std::string("_true"));

      // Start with "fake efficiency corrected data"
      PlotUtils::MnvH1D* h_mc_cross_section = (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());

      //========================================================================
      // Prepare normalizations -- flux, targets
      //========================================================================
        // Flux normalization
        PlotUtils::MnvH1D* h_flux_normalization = 
            (PlotUtils::MnvH1D*)h_mc_cross_section->Clone("flux_normalization");
        h_flux_normalization->ClearAllErrorBands();
        h_flux_normalization->Reset();

        // Get the flux histo, to be integrated
        const bool use_hundred_universes = true;
        static PlotUtils::FluxReweighter* frw = 
            new PlotUtils::FluxReweighter( 14, false, FluxReweighter::minervame1A,
                                           PlotUtils::FluxReweighter::gen2thin, 
                                           PlotUtils::FluxReweighter::g4numiv6,
                                           use_hundred_universes );

        h_flux_normalization = 
            frw->GetIntegratedFluxReweighted(14, h_mc_cross_section, 0., 100.);

        // remove redundant error bands from flux integral
        h_flux_normalization->PopVertErrorBand("Flux_BeamFocus");
        h_flux_normalization->PopVertErrorBand("ppfx1_Total");   

        // Convert flux units from nu/m^2/POT to nu/cm^2/POT
        // Or could we already by in cm^-2?
        //h_flux_normalization->Scale( 1.0e-4 );

        // targets and POT norm
        static const double apothem    = 865.;
        static const double upstream   = 5900.; // ~module 25 plane 1
        static const double downstream = 8430.; // ~module 81 plane 1

        double n_target_nucleons = 
            PlotUtils::TargetUtils::Get().GetTrackerNNucleons(upstream, downstream,
                                                              /*isMC =*/ false,
                                                              apothem);

        // Summarize scales
        std::cout << "  flux_integral cv = "  << h_flux_normalization->GetBinContent(1) << "\n";
        std::cout << "  N target nucleons = " << n_target_nucleons                      << "\n";
        std::cout << "  mc pot = " << util.m_mc_pot << "\n";
        std::cout << "  flux_integral * N target nucleons = " << n_target_nucleons * h_flux_normalization->GetBinContent(1) << "\n";
        std::cout << "  flux_integral * N target nucleons * mc pot= " << n_target_nucleons * h_flux_normalization->GetBinContent(1) * util.m_mc_pot << "\n";
        //static const double data_scale = 1.0 / 3.529606e+42;

      //========================================================================
      // do the normalizations
      //========================================================================
      // Flux
      h_mc_cross_section->AddMissingErrorBandsAndFillWithCV(*h_flux_normalization);
      h_mc_cross_section->Divide( h_mc_cross_section, h_flux_normalization );

      // Targets & Bin width normalization
      static const double mc_scale = 1.0 / ( n_target_nucleons  * util.m_mc_pot );
      //static const double mc_scale = 1.0 / 3.529606e+42;
      h_mc_cross_section->Scale(mc_scale, "width");

      //========================================================================
      // Compare with GenieXSecExtractor
      //========================================================================
      double mc_integral = h_mc_cross_section->Integral();

      if(var->Name() == "pmu") {
        TFile fin_gxse("/minerva/app/users/bmesserl/cmtuser/Minerva_v21r1p1/GENIEXSecExtract/GenieXSecs_20190904_ME1A_pmu.root");

        // Need to convert GXSE result from GeV --> MeV
        PlotUtils::MnvH1D* pmu_xsec_dummy = (PlotUtils::MnvH1D*)fin_gxse.Get("ds_dpmu_xsec");
        assert(pmu_xsec_dummy );
        PlotUtils::MnvH1D* pmu_xsec = (PlotUtils::MnvH1D*)h_mc_cross_section->Clone(uniq());
        pmu_xsec->Reset();
        for(int i = 0; i < pmu_xsec->GetNbinsX()+1; ++i) {
          std::cout << i << "  " << pmu_xsec_dummy->GetBinLowEdge(i) << "  " << pmu_xsec_dummy->GetBinContent(i) << "\n";
          pmu_xsec->SetBinContent(i, pmu_xsec_dummy->GetBinContent(i));
        }

        // What units is the flux in?
        // Maybe we need to convert flux units from nu/cm^2/POT to nu/m^2/POT?
        //pmu_xsec->Scale(1./10000. );

        // Compare integrals
        double gxse_integral = pmu_xsec->Integral();
        std::cout << "  mc xsec integral = "   << mc_integral   << "\n";
        std::cout << "  gxse xsec integral = " << gxse_integral << "\n";
        std::cout << "  gxse / mc = " << gxse_integral / mc_integral << "\n";

        // Area normalize
        h_mc_cross_section->Scale(1./mc_integral);
        pmu_xsec->Scale(1./gxse_integral);
        
        // plot on top of each other
        PlotTogether(h_mc_cross_section, "mc", pmu_xsec, "gxse", "gxse_compare_pmu");

        // Plot ratio
        PlotRatio(h_mc_cross_section, pmu_xsec, Form("GXSEClosure_%s", name));
      }


      //========================================================================
      // Compare EVENT RATE with GenieXSecExtractor
      //========================================================================
      // Start with "fake efficiency corrected data"
      // I've made sure that this matches what we get when we unfold and effcor
      // the bgsub.
      if(var->Name() == "pmu") {
        PlotUtils::MnvH1D* h_all_signal_true = (PlotUtils::MnvH1D*)true_var->m_hists.m_effden.hist->Clone(uniq());

        TFile fin_gxse("/minerva/app/users/bmesserl/cmtuser/Minerva_v21r1p1/GENIEXSecExtract/GenieXSecs_20190904_ME1A_pmu.root");

        // Need to convert GXSE result from GeV --> MeV
        PlotUtils::MnvH1D* pmu_rate_dummy = (PlotUtils::MnvH1D*)fin_gxse.Get("ds_dpmu_xsec_evRate");
        PlotUtils::MnvH1D* pmu_rate = (PlotUtils::MnvH1D*)h_all_signal_true->Clone(uniq());
        pmu_rate->Reset();
        for(int i = 0; i < pmu_rate->GetNbinsX()+1; ++i) {
          std::cout << i << "  " << pmu_rate_dummy->GetBinLowEdge(i) << "  " << pmu_rate_dummy->GetBinContent(i) << "\n";
          pmu_rate->SetBinContent(i, pmu_rate_dummy->GetBinContent(i));
        }

        // BWN
        pmu_rate->Scale(1.,"width");
        h_all_signal_true->Scale(1., "width");

        // Compare integrals
        double mc_rate_integral = h_all_signal_true->Integral();
        double gxse_rate_integral = pmu_rate->Integral();
        std::cout << "  mc rate integral = "   << mc_rate_integral   << "\n";
        std::cout << "  gxse rate integral = " << gxse_rate_integral << "\n";
        std::cout << "  gxse / mc = " << gxse_rate_integral / mc_rate_integral << "\n";

        // Area normalize
        h_all_signal_true->Scale(1./mc_rate_integral);
        pmu_rate->Scale(1./gxse_rate_integral);
        
        // plot on top of each other
        PlotTogether(h_all_signal_true, "mc", pmu_rate, "gxse", "gxse_ratecompare_pmu");

        // Plot ratio
        PlotRatio(h_all_signal_true, pmu_rate, Form("GXSERateClosure_%s", name));
      }
      
  }

  //============================================================================
}

#endif // GXSEClosure_C