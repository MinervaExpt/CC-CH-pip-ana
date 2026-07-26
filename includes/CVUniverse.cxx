#ifndef CVUniverse_cxx
#define CVUniverse_cxx

#include "CVUniverse.h"

#include <TVector3.h>

#include <algorithm>  // max_element
#include <cmath>      //isfinite

#include "PlotUtils/MnvTuneSystematics.h"
#include "utilities.h"  // FixAngle

//==============================================================================
// Constructor
//==============================================================================
CVUniverse::CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
    : PlotUtils::MinervaUniverse(chw, nsigma) {}

//==============================================================================
// Print arachne link
//==============================================================================
void CVUniverse::PrintArachneLink() const {
  int link_size = 200;
  char link[link_size];
  int run = GetInt("mc_run");
  int subrun = GetInt("mc_subrun");
  int gate = GetInt("mc_nthEvtInFile") + 1;
  int slice = GetVecElem("slice_numbers", 0);
  sprintf(link,
          "http://minerva05.fnal.gov/Arachne/"
          "arachne.html\?det=SIM_minerva&recoVer=v21r1p1&run=%d&subrun=%d&gate="
          "%d&slice=%d",
          run, subrun, gate, slice);
  // strncpy(); // FAIL
  // memcpy();  // FAIL
  std::cout << link << std::endl;
}

//==============================================================================
// Dummy access for variable constructors
//==============================================================================
double CVUniverse::GetDummyVar() const { return -999.; }
double CVUniverse::GetDummyHadVar(const int x) const { return -999.; }

//==============================================================================
// Get and set pion candidates
//==============================================================================
TruePionIdx CVUniverse::GetHighestEnergyTruePionIndex() const {
  std::vector<double> tpi_vec = GetTpiTrueVec();  // pip and pim
  const int n_true_pions = GetNChargedPionsTrue();
  // if (n_true_pions != tpi_vec.size()) {
  //  std::cerr << "GetHighestEnergyTruePionIndex Error: True pion size
  //  mismatch\n"; std::exit(1);
  //}

  TruePionIdx reigning_idx = -1;
  double reigning_tpi = 0;
  for (TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
    if (tpi_vec[idx] > reigning_tpi && GetPiChargeTrue(idx) > 0.) {
      reigning_idx = idx;
      reigning_tpi = tpi_vec[idx];
    }
  }

  // This solution below doesn't account for pion charge.
  // see cache/AlgTests for some better attempts at std::alg/lambda solutions
  // TruePionIdx reigning_idx = distance(tpi_vec, max_element(tpi_vec.begin(),
  // tpi_vec.end()));

  return reigning_idx;
}

// -2 --> empty input pion candidates vector
// -1 --> vertex pion
// >= 0 --> index of tracked hadron
// ** No pions with KE > 0 --> fail
// ** TODO what happens when MBR and original method disagree
int CVUniverse::GetHighestEnergyPionCandidateIndex(
    const std::vector<int>& pion_candidate_idxs) const {
  if (pion_candidate_idxs.empty()) {
    return CCNuPionIncConsts::kEmptyPionCandidateVector;  // == -2
  }

  const int dummy_idx = -99;
  int largest_tpi_idx = dummy_idx;
  double largest_tpi = 0.;
  for (uint iter = 0; iter < pion_candidate_idxs.size(); ++iter) {
    int current_idx = pion_candidate_idxs[iter];
    double current_tpi = -997.;
    if (current_idx == CCNuPionIncConsts::kIsVertexPion) {
      // std::cerr << "GetHighestEnergyPionCandidateIndex: pion_idx = -1.\n"
      //             "In the future this will be the code for a vertex pion.\n";
      // std::exit(2);
      double bogus_range = 1.e10;
      current_tpi = GetTpiUntracked(bogus_range);
    } else {
      current_tpi = GetTpi(current_idx);
      if (current_tpi > largest_tpi) {
        largest_tpi_idx = current_idx;
        largest_tpi = current_tpi;
      }
    }
  }
  if (largest_tpi_idx == dummy_idx) {
    // return pion_candidate_idxs[0];
#ifndef NDEBUG
    std::cerr << "GetHighestEnergyPionCandidateIndex: no pion with KE > 0!\n";
#endif
    return -3;
  }
  return largest_tpi_idx;
}

std::vector<RecoPionIdx> CVUniverse::GetPionCandidates() const {
  return m_pion_candidates;
}

void CVUniverse::SetPionCandidates(std::vector<RecoPionIdx> c) {
  m_pion_candidates = c;
  SetNonCalIndices(c);  // for part response syst -- particle(s) that we've
                        // reco-ed by tracking and not by calorimetry
}

void CVUniverse::SetPassesTrakedTracklessCuts(
    bool passesTrackedCuts, bool passesTracklessCuts, bool tracked_sideband,
    bool trackless_sideband, bool tracked_all_ex_w, bool trackless_all_ex_w) {
  m_passesTrackedCuts = passesTrackedCuts;
  m_passesTracklessCuts = passesTracklessCuts;
  m_passesTrackedSideband = tracked_sideband;
  m_passesTracklessSideband = trackless_sideband;
  m_passesTrackedExceptW = tracked_all_ex_w;
  m_passesTracklessExceptW = trackless_all_ex_w;
}

//==============================================================================
// Analysis Variables
//==============================================================================
// muon
double CVUniverse::GetPTmu() const {
  return sqrt(pow(GetPXmu(), 2.0) + pow(GetPYmu(), 2.0));
}

double CVUniverse::GetPXmu() const { return GetMuon4V().Px(); }

double CVUniverse::GetPYmu() const { return GetMuon4V().Py(); }

double CVUniverse::GetPZmu() const { return GetMuon4V().Pz(); }

double CVUniverse::GetThetamuDeg() const {
  return ConvertRadToDeg(GetThetamu());
}

// event-wide
double CVUniverse::GetEhad() const {
  if (m_passesTrackedCuts) {
    return GetCalRecoilEnergy() + GetTrackRecoilEnergy();
  } else if (m_passesTracklessCuts)
    return GetEavail() + GetNMichels() * MinervaUnits::M_pion;
  else if (m_passesTrackedSideband)
    return GetCalRecoilEnergy() + GetTrackRecoilEnergy();
  else if (m_passesTracklessSideband)
    return GetEavail() + GetNMichels() * MinervaUnits::M_pion;
  else if (m_passesTrackedExceptW)
    return GetCalRecoilEnergy() + GetTrackRecoilEnergy();
  else if (m_passesTracklessExceptW)
    return GetEavail() + GetNMichels() * MinervaUnits::M_pion;
  else {
    std::cout << "CVUniverse::GetEhad It is not passing the correct there "
                 "something wrong with the cuts True\n";
    std::exit(1);
  }

  //  return GetEavail() + GetNMichels() * MinervaUnits::M_pion;
  //  return GetCalRecoilEnergy() + GetTrackRecoilEnergy();
}
double CVUniverse::GetEnu() const { return GetEmu() + GetEhad(); }

double CVUniverse::GetQ2() const {
  return CalcQ2(GetEnu(), GetEmu(), GetThetamu());
}

double CVUniverse::GetWexp() const {
  return CalcWexp(GetQ2(), GetEhad());
  /*  if (m_passesTrackedCuts){
      return GetTrackedWexp();
    }
    else if (m_passesTracklessCuts)
      return GetTracklessWexp();
    else if (m_passesTrackedSideband)
      return GetTrackedWexp();
    else if (m_passesTracklessSideband)
      return GetTracklessWexp();
    else {
      std::cout << "CVUniverse::GetMixedWexp It is not passing the correct there
    something wrong with the cuts True\n"; std::exit(1);
    }*/
}

double CVUniverse::GetTracklessWexp() const {
  return CalcWexp(
      CalcQ2(GetEmu() + GetEavail() + GetNMichels() * MinervaUnits::M_pion,
             GetEmu(), GetThetamu()),
      GetEavail() + GetNMichels() * MinervaUnits::M_pion);
}

double CVUniverse::GetTrackedWexp() const {
  return CalcWexp(
      CalcQ2(GetEmu() + GetCalRecoilEnergy() + GetTrackRecoilEnergy(), GetEmu(),
             GetThetamu()),
      GetCalRecoilEnergy() + GetTrackRecoilEnergy());
}

double CVUniverse::Getq0() const { return Calcq0(GetEnu(), GetEmu()); }

double CVUniverse::GetTracklessq3() const {
  return Calcq3(
      CalcQ2(GetEmu() + GetEavail() + GetCalRecoilEnergy() +
                 GetTrackRecoilEnergy(),
             GetEmu(), GetThetamu()),
      GetEmu() + GetEavail() + GetCalRecoilEnergy() + GetTrackRecoilEnergy(),
      GetEmu());
}

// pion
// The output 1.1 means L, the ouput 2.1 means R and the 0 means
// that it is coplanar
double CVUniverse::GetALR(RecoPionIdx hadron) const {
  TVector3 NeuDir(0., 0., 1.);
  TVector3 MuDir(GetPXmu(), GetPYmu(), GetPZmu());
  TVector3 PiDir(GetVecElem("MasterAnaDev_pion_Px", hadron),
                 GetVecElem("MasterAnaDev_pion_Py", hadron),
                 GetVecElem("MasterAnaDev_pion_Pz", hadron));
  TVector3 PlaneDir = NeuDir.Cross(MuDir);
  double proy = PlaneDir.Dot(PiDir);
  if (proy == 0) return 0.1;
  if (proy < 0)
    return 1.1;
  else
    return 2.1;
}

double CVUniverse::GetAdlerCosTheta(RecoPionIdx hadron) const {
  double mumom = GetPmu();
  double pimom = GetVecElem("MasterAnaDev_pion_P", hadron);
  double Enu = GetEnu();
  TVector3 NeuDir(0., 0.057564027, 0.998341817);
  TVector3 MuDir(GetPXmu(), GetPYmu(), GetPZmu());
  MuDir = MuDir.Unit();
  TVector3 PiDir(GetVecElem("MasterAnaDev_pion_Px", hadron),
                 GetVecElem("MasterAnaDev_pion_Py", hadron),
                 GetVecElem("MasterAnaDev_pion_Pz", hadron));
  PiDir = PiDir.Unit();
  TVector3 AdAngle = AdlerAngle(2, mumom /*GeV*/, pimom /*GeV*/, NeuDir, MuDir,
                                PiDir, Enu /*GeV*/);
  return cos(AdAngle[1]);
}

double CVUniverse::GetAdlerPhi(RecoPionIdx hadron) const {
  double mumom = GetPmu();
  double pimom = GetVecElem("MasterAnaDev_pion_P", hadron);
  double Enu = GetEnu();
  TVector3 NeuDir(0., 0.057564027, 0.998341817);
  TVector3 MuDir(GetPXmu(), GetPYmu(), GetPZmu());
  MuDir = MuDir.Unit();
  TVector3 PiDir(GetVecElem("MasterAnaDev_pion_Px", hadron),
                 GetVecElem("MasterAnaDev_pion_Py", hadron),
                 GetVecElem("MasterAnaDev_pion_Pz", hadron));
  PiDir = PiDir.Unit();
  TVector3 AdAngle = AdlerAngle(2, mumom /*GeV*/, pimom /*GeV*/, NeuDir, MuDir,
                                PiDir, Enu /*GeV*/);
  return AdAngle[2];
}

// dEdx tool w assumption that track is pion
double CVUniverse::GetEpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_E", hadron);
}

double CVUniverse::GetPT(RecoPionIdx hadron) const {
  TVector3 pT_mu(GetPXmu(), GetPYmu(), 0);
  TVector3 pT_pi(GetVecElem("MasterAnaDev_pion_Px", hadron),
                 GetVecElem("MasterAnaDev_pion_Py", hadron), 0);
  TVector3 pT = pT_mu + pT_pi;
  return pT.Mag();
}

double CVUniverse::GetPXpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_Px", hadron);
}

double CVUniverse::GetPYpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_Py", hadron);
}

double CVUniverse::GetPZpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_Pz", hadron);
}

double CVUniverse::GetPpi(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_P", hadron);
}

double CVUniverse::GetThetapi(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cerr << "CVU::GetThetapi: pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n"
                 "In that case, this function won't make sense.\n";
    throw hadron;
  }
  return FixAngle(GetVecElem("MasterAnaDev_pion_theta", hadron));
}

double CVUniverse::GetThetapiDeg(RecoPionIdx hadron) const {
  return ConvertRadToDeg(GetThetapi(hadron));
}

double CVUniverse::GetTpi(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cerr << "CVU::GetTpi: pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n";
    throw hadron;
  }
  // return (GetVecElem("MasterAnaDev_hadron_pion_E", hadron)
  //          - CCNuPionIncConsts::CHARGED_PION_MASS)/0.96;
  double Epi = GetVecElem("MasterAnaDev_pion_E", hadron);
  if (Epi == 0) {
    return GetTpiMBR(hadron);  // TODO maybe do from momentum instead
                               // Not sure what this fail mode means.
    // std::cout << "CVUniverse::GetTpi: warning Epi < 0\n";
    // return 0;
  }
  return Epi - CCNuPionIncConsts::CHARGED_PION_MASS;
}

double CVUniverse::GetTpiMBR(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cerr << "CVU::GetTpiMBR: pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n"
                 "In that case, this function won't make sense.\n";
    throw hadron;
  }
  double TLA = GetVecElem("hadron_track_length_area", hadron);
  return 2.3112 * TLA + 37.03;
}

double CVUniverse::GetpimuAngle(
    RecoPionIdx hadron) const {  // Angle beetwen P_pi and P_mu (degrees)
  TVector3 p_mu(GetPXmu(), GetPYmu(), GetPZmu());
  TVector3 p_pi(GetVecElem("MasterAnaDev_pion_Px", hadron),
                GetVecElem("MasterAnaDev_pion_Py", hadron),
                GetVecElem("MasterAnaDev_pion_Pz", hadron));
  double PidotMu = p_pi.Dot(p_mu);
  double Pmu = p_mu.Mag(), Ppi = p_pi.Mag();
  return ConvertRadToDeg(acos((PidotMu) / (Pmu * Ppi)));
}

double CVUniverse::Gett(RecoPionIdx h) const {
  if (h < 0) {
    std::cerr << "CVU::Gett: pion_idx < 0.\n"
                 "This indicates no pion track, so you shouldn't be here.\n";
    throw h;
  }
  ROOT::Math::PxPyPzEVector mu4v = GetMuon4V();
  return Calct(GetPXpi(h), GetPYpi(h), GetPZpi(h), GetEpi(h), mu4v.Px(),
               mu4v.Py(), mu4v.Pz(), GetEmu());
}

double CVUniverse::GetTpiTrackless() const {
  return GetTpiUntracked(m_vtx_michels.m_bestdist);
}

double CVUniverse::GetBestDistance() const { return m_vtx_michels.m_bestdist; }

double CVUniverse::GetMixedTpi(RecoPionIdx idx) const {
  if (m_passesTrackedCuts) {
    return GetTpi(idx);
  } else if (m_passesTracklessCuts)
    return GetTpiTrackless();
  else if (m_passesTrackedSideband)
    return GetTpi(idx);
  else if (m_passesTracklessSideband)
    return GetTpiTrackless();
  else if (m_passesTrackedExceptW)
    return GetTpi(idx);
  else if (m_passesTracklessExceptW)
    return GetTpiTrackless();
  else {
    std::cout << "CVUniverse::GetMixedTpi It is not passing the correctly the "
                 "cuts info, there something wrong with the cuts Reco\n";
    std::exit(1);
  }
}

double CVUniverse::GetThetapitrackless() const {
  return m_vtx_michels.m_bestthetaangle;
}

double CVUniverse::GetThetapitracklessDeg() const {
  return ConvertRadToDeg(GetThetapitrackless());
}

double CVUniverse::GetMixedThetapiDeg(RecoPionIdx Idx) const {
  if (m_passesTrackedCuts) {
    return GetThetapiDeg(Idx);
  } else if (m_passesTracklessCuts)
    return GetThetapitracklessDeg();
  else if (m_passesTrackedSideband)
    return GetThetapiDeg(Idx);
  else if (m_passesTracklessSideband)
    return GetThetapitracklessDeg();
  else if (m_passesTrackedExceptW)
    return GetThetapiDeg(Idx);
  else if (m_passesTracklessExceptW)
    return GetThetapitracklessDeg();
  else {
    std::cout
        << "CVUniverse::GetMixedThetapideg It is not passing the correctly the "
           "cuts info, there something wrong with the cuts Reco\n";
    std::exit(1);
  }
}

double CVUniverse::GetEavail() const {
  double recoiltracker = GetDouble("blob_recoil_E_tracker") -
                         GetTrackerECALMuFuzz()[0];  // GetTrackerMuFuzz();
  double recoilEcal = GetDouble("blob_recoil_E_ecal") -
                      GetTrackerECALMuFuzz()[1];  // GetECALMuFuzz();
  const double Eavailable_scale = 1.17;
  double eavail = recoiltracker + recoilEcal;
  return eavail * Eavailable_scale;
}

//==============================================================================
// Truth
//==============================================================================
// The output 1 means L, the ouput 2 means R and the 0 means
// that it is coplanar
double CVUniverse::GetALRTrue(TruePionIdx idx) const {
  TVector3 NeuDir(GetVecElem("mc_incomingPartVec", 0),
                  GetVecElem("mc_incomingPartVec", 1),
                  GetVecElem("mc_incomingPartVec", 2));
  TVector3 MuDir(GetPXmuTrue(), GetPYmuTrue(), GetPZmuTrue());
  TVector3 PiDir(GetVecElem("truth_pi_px", idx), GetVecElem("truth_pi_py", idx),
                 GetVecElem("truth_pi_pz", idx));
  TVector3 PlaneDir = NeuDir.Cross(MuDir);
  double proy = PlaneDir.Dot(PiDir);
  if (proy == 0) return 0.1;
  if (proy < 0)
    return 1.1;
  else
    return 2.1;
}

double CVUniverse::GetAdlerCosThetaTrue(TruePionIdx idx) const {
  double mumom = GetPlepTrue();
  double Enu = GetEnuTrue();
  TVector3 NeuDir(GetVecElem("mc_incomingPartVec", 0),
                  GetVecElem("mc_incomingPartVec", 1),
                  GetVecElem("mc_incomingPartVec", 2));
  NeuDir = NeuDir.Unit();
  TVector3 MuDir(GetPXmuTrue(), GetPYmuTrue(), GetPZmuTrue());
  MuDir = MuDir.Unit();
  TVector3 PiDir(GetVecElem("truth_pi_px", idx), GetVecElem("truth_pi_py", idx),
                 GetVecElem("truth_pi_pz", idx));
  double pimom = PiDir.Mag();
  PiDir = PiDir.Unit();
  TVector3 AdAngle = AdlerAngle(2, mumom /*GeV*/, pimom /*GeV*/, NeuDir, MuDir,
                                PiDir, Enu /*GeV*/);
  return cos(AdAngle[1]);
}

double CVUniverse::GetAdlerPhiTrue(TruePionIdx idx) const {
  double mumom = GetPlepTrue();
  double Enu = GetEnuTrue();
  TVector3 NeuDir(GetVecElem("mc_incomingPartVec", 0),
                  GetVecElem("mc_incomingPartVec", 1),
                  GetVecElem("mc_incomingPartVec", 2));
  NeuDir = NeuDir.Unit();
  TVector3 MuDir(GetPXmuTrue(), GetPYmuTrue(), GetPZmuTrue());
  MuDir = MuDir.Unit();
  TVector3 PiDir(GetVecElem("truth_pi_px", idx), GetVecElem("truth_pi_py", idx),
                 GetVecElem("truth_pi_pz", idx));
  double pimom = PiDir.Mag();
  PiDir = PiDir.Unit();
  TVector3 AdAngle = AdlerAngle(2, mumom /*GeV*/, pimom /*GeV*/, NeuDir, MuDir,
                                PiDir, Enu /*GeV*/);
  return AdAngle[2];
}

// Need a function that gets truth tracked energy, total E for pions, just KE
// for protons.
double CVUniverse::GetAllTrackEnergyTrue() const {
  double etracks = 0;
  for (const auto& pi_idx : GetPionCandidates()) {
    etracks += GetVecElem("MasterAnaDev_hadron_tm_beginKE", pi_idx);
    // std::cout << GetVecElem("MasterAnaDev_hadron_tm_PDGCode", pi_idx) << " ";
    if (abs(GetVecElem("MasterAnaDev_hadron_tm_PDGCode", pi_idx)) ==
        211)  // TODO may want to only not add pion mass when is proton or
              // neutron
      etracks += CCNuPionIncConsts::CHARGED_PION_MASS;
  }
  // std::cout << "\n";
  return etracks;
}

double CVUniverse::GetEmuTrue() const { return GetElepTrue(); }

double CVUniverse::GetIntVtxXTrue() const { return GetVecElem("mc_vtx", 0); }

double CVUniverse::GetIntVtxYTrue() const { return GetVecElem("mc_vtx", 1); }

double CVUniverse::GetIntVtxZTrue() const { return GetVecElem("mc_vtx", 2); }

double CVUniverse::GetPTTrue(TruePionIdx idx) const {
  TVector3 pT_mu(GetPXmuTrue(), GetPYmuTrue(), 0);
  TVector3 pT_pi(GetVecElem("truth_pi_px", idx), GetVecElem("truth_pi_py", idx),
                 0);
  TVector3 pT = pT_mu + pT_pi;
  return pT.Mag();
}

double CVUniverse::GetPTmuTrue() const {
  return GetPlepTrue() * sin(GetThetalepTrue());
}

double CVUniverse::GetPXmuTrue() const {
  return GetVecElem("mc_primFSLepton", 0);
}

double CVUniverse::GetPYmuTrue() const {
  return GetVecElem("mc_primFSLepton", 1);
}

double CVUniverse::GetPZmuTrue() const {
  return GetPlepTrue() * cos(GetThetalepTrue());
}

double CVUniverse::GetPmuTrue() const { return GetPlepTrue(); }

double CVUniverse::GetThetamuTrue() const {
  return FixAngle(GetThetalepTrue());
}

double CVUniverse::GetThetamuTrueDeg() const {
  return ConvertRadToDeg(GetThetamuTrue());
}

// input: truth index. output: may be pip or pim.
double CVUniverse::GetThetapiTrue(TruePionIdx idx) const {
  double t_pi_theta = GetVecElem("truth_pi_theta_wrtbeam", idx);
  if (t_pi_theta == -9.0) {
    std::cerr << "CVU::GetThetapiTrue: Default angle.\n"
                 "Tried to access truth pion angle for a nonexistent "
                 "truth pion trajectory.\n";
    throw t_pi_theta;
  }
  return FixAngle(t_pi_theta);
}

double CVUniverse::GetThetapiTrueDeg(TruePionIdx idx) const {
  return ConvertRadToDeg(GetThetapiTrue(idx));
}

// Get TRUE Hadron Quantities (always MeV, radians)
// NOTE: These 'truth_pi_*' containers all have 20 elements
// Their non-default elements correspond to all of and only the CHARGED
// (N.B. + & -) pions from the event's TG4Trajectory primary
// trajectories.
// Pion indices here do not correspond to reco hadron indices.
// Be careful with these functions. Make sure you understand the input
// and output.
// input: truth index. output: may be pip or pim.
double CVUniverse::GetTpiTrue(TruePionIdx idx) const {
  double t_pi_E = GetVecElem("truth_pi_E", idx);
  if (t_pi_E == -1.) {
    std::cerr << "CVU::GetTpiTrue: Default energy.\n"
                 "Tried to access truth pion energy for a nonexistent "
                 "truth pion trajectory.\n";
    throw t_pi_E;
  }
  return t_pi_E - CCNuPionIncConsts::CHARGED_PION_MASS;
}

double CVUniverse::GetWexpTrue() const {
  return CalcWexp(GetQ2True(), GetEhadTrue());
}

double CVUniverse::GetWgenie() const { return GetDouble("mc_w"); }

double CVUniverse::GetpimuAngleTrue(
    TruePionIdx idx) const {  // Angle beetwen P_pi and P_mu (degrees)
  TVector3 p_mu(GetPXmuTrue(), GetPYmuTrue(), GetPZmuTrue());
  TVector3 p_pi(GetVecElem("truth_pi_px", idx), GetVecElem("truth_pi_py", idx),
                GetVecElem("truth_pi_pz", idx));
  double PidotMu = p_pi.Dot(p_mu);
  double Pmu = p_mu.Mag(), Ppi = p_pi.Mag();
  return ConvertRadToDeg(acos((PidotMu) / (Pmu * Ppi)));
}

int CVUniverse::GetNChargedPionsTrue() const {
  return GetInt("truth_N_pip") + GetInt("truth_N_pim");
}

// input: truth index. output: may be pip or pim.
int CVUniverse::GetPiChargeTrue(TruePionIdx idx) const {
  int t_pi_charge = GetVecElem("truth_pi_charge", idx);
  if (t_pi_charge == 0) {
    std::cerr << "CVU::GetPiChargeTrue: Default charge.\n"
                 "Tried to access truth pion charge for a nonexistent "
                 "truth pion trajectory.\n";
    throw t_pi_charge;
  }
  return t_pi_charge;
}

// !!!!!! Vector filled with pip and pim !!!!!!!
std::vector<double> CVUniverse::GetTpiTrueVec() const {
  std::vector<double> ret;
  const int n_true_pions = GetNChargedPionsTrue();
  for (TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
    ret.push_back(GetTpiTrue(idx));
  }
  return ret;
}

double CVUniverse::GetMixedTpiTrue(TruePionIdx idx) const {
  if (m_passesTrackedCuts) {
    return GetTpiTrue(idx);
  } else if (m_passesTracklessCuts)
    return GetTrueTpi();
  else if (m_passesTrackedSideband)
    return GetTpiTrue(idx);
  else if (m_passesTracklessSideband)
    return GetTrueTpi();
  else if (m_passesTrackedExceptW)
    return GetTpiTrue(idx);
  else if (m_passesTracklessExceptW)
    return GetTrueTpi();
  else {
    std::cout << "CVUniverse::GetMixedTpiTrue It is not passing the correct "
                 "there something wrong with the cuts True\n";
    std::exit(1);
  }
}

double CVUniverse::GetThetapitracklessTrueDeg() const {
  return ConvertRadToDeg(GetThetapitracklessTrue());
}

double CVUniverse::GetThetapitracklessTrue() const {
  std::vector<double> tpi_vec = GetTpiTrueVec();  // pip and pim
  const int n_true_pions = GetNChargedPionsTrue();
  TruePionIdx reigning_idx = -1;
  double reigning_tpi = 9999;
  for (TruePionIdx idx = 0; idx < n_true_pions; ++idx) {
    if (tpi_vec[idx] < reigning_tpi && GetPiChargeTrue(idx) > 0.) {
      reigning_idx = idx;
      reigning_tpi = tpi_vec[idx];
    }
  }
  return GetThetapiTrue(reigning_idx);
}

double CVUniverse::GetMixedThetapiTrueDeg(TruePionIdx idx) const {
  if (m_passesTrackedCuts) {
    return GetThetapiTrueDeg(idx);
  } else if (m_passesTracklessCuts)
    return GetThetapitracklessTrueDeg();
  else if (m_passesTrackedSideband)
    return GetThetapiTrueDeg(idx);
  else if (m_passesTracklessSideband)
    return GetThetapitracklessTrueDeg();
  else if (m_passesTrackedExceptW)
    return GetThetapiTrueDeg(idx);
  else if (m_passesTracklessExceptW)
    return GetThetapitracklessTrueDeg();
  else {
    std::cout << "CVUniverse::GetMixedThetapiTrueDeg It is not passing the "
                 "correct there something wrong with the cuts True\n";
    std::exit(1);
  }
}
//==============================
// Ehad (GetErecoil) Variables
//==============================
// If the calorimetric energy is too great, abandon trying to calculate tracked
// energy separately. This requires a coordinated effort from both functions.

// Passive-corrected pion energy
// NukeCCPion_pion_recoilE_passive =
// m_caloUtils->applyCalConsts(hadronProng,"Default",true,true);
double CVUniverse::GetCalEpi(int iProng) const {
  return GetVecElem("MasterAnaDev_hadron_pion_E_recoil_corr", iProng);
}

// Untracked recoil energy
double CVUniverse::GetCalRecoilEnergy() const {
  const double ecal_nopi = GetCalRecoilEnergyNoPi_DefaultSpline();
  if (ecal_nopi > 1000)
    return GetCalRecoilEnergy_CCPiSpline();
  else
    return GetCalRecoilEnergyNoPi_Corrected(ecal_nopi);
}

// Apply an additive, ad hoc correction to the CalRecoilENoPi
double CVUniverse::GetCalRecoilEnergyNoPi_Corrected(
    const double ecal_nopi) const {
  if (GetPionCandidates().size() == 0) return ecal_nopi;

  // I've shown that low-t (likely coherent) events don't need the
  // correction. 20210102_ErecStudies, slides 46-48.
  RecoPionIdx best_pion =
      GetHighestEnergyPionCandidateIndex(GetPionCandidates());

  if (best_pion >= 0 && Gett(best_pion) < 125.e3) return ecal_nopi;

  // Otherwise, do make the correction
  double ecal_nopi_corrected = ecal_nopi;

  TArrayD ecal_nopi_bins = CCPi::GetBinning("ecal_nopi");

  const std::vector<double> corrections = {-0.060, -0.050, -0.210, -0.180,
                                           -0.165, -0.180, -0.180, -0.180,
                                           -0.195, -0.360, -0.400};

  // Find the correction corresponding to the nominal ehad-nopi.
  // Set the corrected value to be the (nominal - correction) (where the
  // corrections are all negative, so in effect we're shifting everything
  // UP)
  for (int i_bin = 0; i_bin < ecal_nopi_bins.GetSize() - 1; ++i_bin) {
    if (((1e3) * ecal_nopi_bins[i_bin] < ecal_nopi) &&
        (ecal_nopi < (1e3) * ecal_nopi_bins[i_bin + 1])) {
      if (ecal_nopi < 1000) {
        ecal_nopi_corrected = ecal_nopi - (1e3) * corrections[i_bin];
        break;
      }
    }
  }

  return ecal_nopi_corrected;
}

// Cal recoil energy minus calorimetrically-measured pion energy.
// Used to determined whether we should try to use the correction or not.
double CVUniverse::GetCalRecoilEnergyNoPi_DefaultSpline() const {
  double nopi_recoilE = GetCalRecoilEnergy_DefaultSpline();
  for (const auto& pi_idx : GetPionCandidates()) {
    nopi_recoilE -= GetCalEpi(pi_idx);
  }
  return nopi_recoilE;
}

// Total recoil with CCPi spline correction.
// Spline measured from: CC, 1pi+, 1mu, NBaryons,
// True W exp < 1.4 GeV, thmu_true < 20 deg, 35 < tpi < 350 MeV, Minos match
// RecoilUtils->calcRecoilEFromClusters(event, muonProng,
// "NukeCCPion_TwoTrack_Nu_Tracker");
double CVUniverse::GetCalRecoilEnergy_CCPiSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_two_track");
}

// RecoilUtils->calcRecoilEFromClusters(event, muonProng, "Default" );
double CVUniverse::GetCalRecoilEnergy_DefaultSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_default");
}

// This is what the response universe calls our tracked recoil energy
double CVUniverse::GetNonCalRecoilEnergy() const {
  if (GetPionCandidates().empty()) {
#ifndef NDEBUG
//    std::cout << "CVU::GetNonCalRecoilEnergy WARNING: no pion candidates!\n";
#endif
    /*    int nMichels = GetNMichels();
        double pimass = 0;
        if (nMichels > 0){
          pimass = nMichels*MinervaUnits::M_pion;
          return pimass;
        }
        else*/
    return 0.;
  }

  double etracks = 0.;

  const double ecal_nopi = GetCalRecoilEnergyNoPi_DefaultSpline();
  if (ecal_nopi > 1000) {
    etracks = 0;
  } else {
    for (const auto& pi_idx : GetPionCandidates()) {
      // std::cout << GetVecElem("MasterAnaDev_hadron_tm_PDGCode", pi_idx) << "
      // ";
      etracks += GetEpi(pi_idx);
    }
    // std::cout << "\n";
  }

  return etracks;
}

// (Tracked) recoil energy, not determined from calorimetry
double CVUniverse::GetTrackRecoilEnergy() const {
  return GetNonCalRecoilEnergy();
}

/*double CVUniverse::GetCalRecoilEnergy() const { //Using the Anezka definition
  return GetDouble("part_response_total_recoil_passive_allNonMuonClusters_id")+
GetDouble("part_response_total_recoil_passive_allNonMuonClusters_od"); // in MeV
}*/
//==============================
// ehad old variables
//==============================
// Old-school attempt at Ehad - Epi
// Subtract dEdx-tool-pion-energy from Ehad-from-calcRecoilEFromClusters
// ("CCInclusive" splines)
double CVUniverse::GetCalRecoilEnergyNoPi_CCIncSpline() const {
  double ecal_nopi = GetCalRecoilEnergy_CCIncSpline();
  for (const auto& pi_idx : GetPionCandidates()) ecal_nopi -= GetEpi(pi_idx);
  return ecal_nopi;
}

// Ehad CCInclusive Spline Variables
// Ehad ccinclusive splines -- doesn't account for pion
double CVUniverse::GetCalRecoilEnergy_CCIncSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_CCInc");
}

//==============================
// ehad truth variables
//==============================
// Given reco pions – subtract their true matched energy from truth ehad = (enu
// - emu)
double CVUniverse::GetCalRecoilEnergyNoPiTrue() const {
  double nopi_recoilE = GetEhadTrue();
  for (const auto& pi_idx : GetPionCandidates())
    nopi_recoilE -= GetEpiTrueMatched(pi_idx);
  return nopi_recoilE;
}

double CVUniverse::GetEhadTrue() const { return GetEnuTrue() - GetElepTrue(); }

// Given a reco track, get the track's truth E, assuming it's a pion
double CVUniverse::GetEpiTrueMatched(RecoPionIdx hadron) const {
  return GetTpiTrueMatched(hadron) + CCNuPionIncConsts::CHARGED_PION_MASS;
}

// Given a reco track, get the track's truth KE
double CVUniverse::GetTpiTrueMatched(RecoPionIdx hadron) const {
  // TruePionIdx true_index = -1;
  // TruePionIdx true_index = GetVecElem("MasterAnaDev_hadron_tm_trackID",
  // hadron); int true_pdg   = GetVecElem("MasterAnaDev_hadron_tm_PDGCode",
  // hadron); if(true_pdg != 211) std::cout << "pion mis-identified as a " <<
  // true_pdg << "!\n";
  return GetVecElem("MasterAnaDev_hadron_tm_beginKE", hadron);
}

//==============================================================================
// Cuts, Systematics, Studies
//==============================================================================
bool CVUniverse::IsInHexagon(double x, double y, double apothem) const {
  double lenOfSide = apothem * (2 / sqrt(3));
  double slope = (lenOfSide / 2.0) / apothem;
  double xp = fabs(x);
  double yp = fabs(y);

  if ((xp * xp + yp * yp) < apothem * apothem)
    return true;
  else if (xp <= apothem && yp * yp < lenOfSide / 2.0)
    return true;
  else if (xp <= apothem && yp < lenOfSide - xp * slope)
    return true;

  return false;
}

bool CVUniverse::IsInPlastic() const {
  if (!IsInHexagon(GetVecElem("mc_vtx", 0), GetVecElem("mc_vtx", 1), 1000.0))
    return false;  // This is in the calorimeters

  double mc_vtx_z = GetVecElem("mc_vtx", 2);
  if (mc_vtx_z > 8467.0) return false;  // Ditto

  int mc_nuclei = GetInt("mc_targetZ");
  // In the carbon target?  The z is gotten from NukeBinningUtils
  if (fabs(mc_vtx_z - 4945.92) <=
          PlotUtils::TargetProp::ThicknessMC::Tgt3::C / 2 &&
      mc_nuclei == 6)
    return false;

  // In the water target?
  if (5200 < mc_vtx_z && mc_vtx_z < 5420 && (mc_nuclei == 8 || mc_nuclei == 1))
    return false;

  // Finally, do you have target material?  I'm going to say lead/iron isn't a
  // big consideration elsewhere in the detector
  if (mc_nuclei == 26 || mc_nuclei == 82) return false;

  return true;
}

double CVUniverse::GetEmichel(RecoPionIdx hadron) const {
  return GetVecElem("has_michel_cal_energy", hadron);
}

double CVUniverse::GetEnode0(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q0", hadron);
}

double CVUniverse::GetEnode1(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q1", hadron);
}

double CVUniverse::GetEnode01(RecoPionIdx hadron) const {
  return GetEnode0(hadron) + GetEnode1(hadron);
}

double CVUniverse::GetEnode2(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q2", hadron);
}

double CVUniverse::GetEnode3(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q3", hadron);
}

double CVUniverse::GetEnode4(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q4", hadron);
}

double CVUniverse::GetEnode5(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_lastnode_Q5", hadron);
}

double CVUniverse::GetPpionCorr(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_hadron_pion_p_corr", hadron);
}

double CVUniverse::GetFitVtxX() const {
  return GetVecElem("MasterAnaDev_vtx", 0);
}  // cm?

double CVUniverse::GetFitVtxY() const {
  return GetVecElem("MasterAnaDev_vtx", 1);
}  // cm?

double CVUniverse::GetFitVtxZ() const {
  return GetVecElem("MasterAnaDev_vtx", 2);
}  // cm?

double CVUniverse::GetLLRScore(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cout << "CVU::GetLLRScore pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n"
                 "In that case, this function won't make sense.\n";
    throw hadron;
  } else if (hadron < -1) {
#ifndef NDEBUG
    std::cerr << "CVU::GetLLRScore bogus pion_idx." << hadron << "\n";
#endif
    return -1;
    // throw hadron;
  }
  return GetVecElem("MasterAnaDev_hadron_piFit_scoreLLR", hadron);
}

double CVUniverse::GetLargestIsoProngSep() const {
  std::vector<double> sep_vec = GetVec<double>("iso_prong_separation");
  return *std::max_element(sep_vec.begin(), sep_vec.end());
}

double CVUniverse::GetLargestPrimProngSep() const {
  std::vector<double> sep_vec = GetVec<double>("prong_separation");
  return *std::max_element(sep_vec.begin(), sep_vec.end());
}

// double CVUniverse::GetTpiFResidual(const int hadron, const bool MBR) const {
double CVUniverse::GetTpiFResidual(const int hadron) const {
  // ALERT 2023-03-10
  // Christian reports that in fact
  // true_index = GetVecElem("MasterAnaDev_hadron_tm_trackID", hadron) - 1;
  // I don't use true_index here, but if I do in the future I need to check it!

  // double T_reco = !MBR ? GetMixedTpi(hadron) : GetTpiMBR(hadron);
  double T_reco = GetMixedTpi(hadron);
  int true_index = -1;
  true_index = GetVecElem("MasterAnaDev_hadron_tm_trackID", hadron);
  double T_true = GetVecElem("MasterAnaDev_hadron_tm_beginKE", hadron);
  double fresid =
      (!std::isfinite(T_reco / T_true) || (std::abs(T_reco / T_true)) > 10000)
          ? 0.96
          : T_reco / T_true - 1.;
  int true_pdg = GetVecElem("MasterAnaDev_hadron_tm_PDGCode", hadron);
  return fresid;
}

double CVUniverse::GetWexpFResidual() const {
  return GetWexp() / GetWexpTrue() - 1.;
  // return GetWexp()/GetWgenie() - 1.;
}

double CVUniverse::GetdEdxScore(RecoPionIdx hadron) const {
  if (hadron == -1) {
    std::cout << "CVU::GetdEdxScore: pion_idx = -1.\n"
                 "In the future this will be the code for a vertex pion.\n"
                 "In that case, this function won't make sense.\n";
    throw hadron;
  }
  return GetVecElem("MasterAnaDev_hadron_piFit_score1", hadron);
}

int CVUniverse::GetNAnchoredLongTracks() const {
  return GetInt("n_anchored_long_trk_prongs");
}

int CVUniverse::GetNAnchoredShortTracks() const {
  return GetInt("n_anchored_short_trk_prongs");
}

int CVUniverse::GetNIsoProngs() const {
  return GetDouble("n_nonvtx_iso_blobs_all");
  //  return GetDouble("iso_prongs_count");// branch for p3
}

int CVUniverse::GetNNodes(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_pion_nNodes", hadron);
}

int CVUniverse::GetNhadrons() const {
  return GetInt("MasterAnaDev_hadron_number");
}

int CVUniverse::GetTrackReconstructionMethod(RecoPionIdx hadron) const {
  return GetVecElem("MasterAnaDev_hadron_1stTrackPatRec", hadron);
}

//==============================================================================
// Weights
//==============================================================================
double CVUniverse::GetAnisoDeltaDecayWarpWeight() const {
  return GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 4);
}

// Note, this assumes you're not using the diffractive model in GENIE
// As of 03/2021, we don't really trust our diffractive model, so
// as a rough approximation, weight every coherent event
// (diffractive is coherent on hydrogen) by 1.4368.
// Coherent xsec scales by A^(1/3), and 1/(12^(1/3)) = 0.4368
double CVUniverse::GetDiffractiveWeight() const {
  if (GetInt("mc_intType") != 4) return 1.;
  // Note: diffractive should be applied only to plastic. This approximates that
  // if( PlotUtils::TargetUtils::Get().InCarbon3VolMC( GetVecElem("mc_vtx",0),
  //                                 GetVecElem("mc_vtx",1),
  //                                 GetVecElem("mc_vtx",2) ) ) return 1.;
  // if( GetInt("mc_nucleiZ") != 6 ) 1.;
  if (!IsInPlastic() && !PlotUtils::TargetUtils::Get().InWaterTargetMC(
                            GetIntVtxXTrue(), GetIntVtxYTrue(),
                            GetIntVtxZTrue(), GetInt("mc_targetZ")))
    return 1.;

  return 1.4368;
}

double CVUniverse::GetGenieWarpWeight(
    double p) const {  // p is the factor that
                       // the weight is modified
  double wgt = GetVecElem("truth_genie_wgt_MaRES", 4);
  wgt = 1 + (wgt - 1) *
                p;  // double the size of the shift from 1 (e.g. 1.1 --> 1.2)
  return wgt;
}

// This is a systematic shift version of the lowq2 weight.
// We do not use this currently. Instead we use the MAT version, which gives us
// a "CV" weight. Keeping it around in case we want it for a warp study or
// something.
double CVUniverse::GetLowQ2PiWeight(double q2, std::string channel) const {
  if (!PlotUtils::IsCCRes(*this))
    return 1.;
  else
    return PlotUtils::weight_lowq2pi().getWeight(q2, channel, +1);
}

double CVUniverse::GetWeight() const {
  // Warping strategy is to only turn on one of these at a time.
  const bool do_genie_warping = false;  // we are not using this werp anymore
  const bool do_aniso_warping = false;
  const bool do_mk_warping = false;
  const bool do_MARESFrac_warping = false;
  const bool do_tpi_warping = false;
  const bool closureTest = false;

  double wgt_flux = 1., wgt_2p2h = 1.;
  double wgt_rpa = 1., wgt_lowq2 = 1.;
  double wgt_genie = 1., wgt_mueff = 1.;
  double wgt_anisodd = 1.;
  double wgt_michel = 1.;
  double wgt_diffractive = 1.;
  double wgt_mk = 1.;
  double wgt_target = 1.;
  double wgt_CCPiWegiht = 1.;
  double wgt_fsi = 1., wgt_coh = 1., wgt_geant = 1.,
         wgt_sbfit = 1. /* This weight depends of the sidebands, Will we applay
                           this weight?*/
      ;
  double wgt_pionReweight = 1.;
  double wgt_tpi_warp = 1.;
  // genie
  wgt_genie = GetGenieWeight();
  if (do_genie_warping) wgt_genie = GetGenieWarpWeight(2.);
  if (do_MARESFrac_warping) wgt_genie = GetGenieWarpWeight(0.2);

  // flux

  wgt_flux = GetFluxAndCVWeight();

  // rpa
  wgt_rpa = GetRPAWeight();
  if (!closureTest) {
    // MINOS efficiency
    if (!m_is_truth && GetInt("isMinosMatchTrack") == 1 &&
        GetInt("MasterAnaDev_nuHelicity") == 1) {  // Remove for closure test
      wgt_mueff = GetMinosEfficiencyWeight();
    }
    // 2p2h
    wgt_2p2h = GetLowRecoil2p2hWeight();

    // low q2
    wgt_lowq2 = (GetQ2True() > 0)  // remove for closure test
                    ? GetLowQ2PiWeight(CCNuPionIncShifts::kLowQ2PiChannel)
                    : 1.;

    // aniso delta decay weight -- currently being used for warping
    if (do_aniso_warping)
      wgt_anisodd = GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 4);

    // Michel efficiency
    wgt_michel = GetMichelEfficiencyWeight();  // remove for closure test

    // Diffractive
    wgt_diffractive = GetDiffractiveWeight();  // remove for closure test

    // MK Weight
    if (do_mk_warping) wgt_mk = GetMKWeight();

    // Target Mass
    wgt_target = GetTargetMassWeight();  // remove for closure test

    // Tpi Mehreen's weight
    wgt_pionReweight = GetUntrackedPionWeight();  // remove for closure test

    if (do_tpi_warping) wgt_tpi_warp = 0.8 + 0.2 * wgt_pionReweight;

    // New Weights added taking as reference Aaron's weights

    wgt_fsi = GetFSIWeight(0);  // Remove for closure test

    if (GetInt("mc_intType") == 4) {  // Remove for closure test
      int idx = (int)GetHighestEnergyTruePionIndex();
      if (GetNChargedPionsTrue() > 1)
        std::cout << " More that one charge pion in Coherent events "
                  << GetNChargedPionsTrue() << "Index "
                  << GetHighestEnergyTruePionIndex() << "\n";
      double deg_theta_pi = GetThetapiTrueDeg(idx);
      if (GetTpiTrue(idx) > 0.)
        wgt_coh *= GetCoherentPiWeight(
            deg_theta_pi, (GetTpiTrue(idx) + MinervaUnits::M_pion) / 1000);
    }

    wgt_geant = GetGeantHadronWeight();  // Remove for closure test
  }

  if (m_is_signal) {
    int idx = (int)GetHighestEnergyTruePionIndex();
    wgt_CCPiWegiht = GetChargedPionTuneWeight(GetTpiTrue(idx), GetQ2True());
  }

  /*  std::cout << "GENIE " << wgt_genie << " Flux " <<  wgt_flux << " 2p2h " <<
              wgt_2p2h << " RPA " << wgt_rpa << " LowQ2 " <<  wgt_lowq2 <<
              " MuEff " << wgt_mueff << " Michel " << wgt_michel << "
     Diffractive "
              << wgt_diffractive << " Target " << wgt_target << " FSI " <<
     wgt_fsi
              << " Coherent " << wgt_coh << " GEANT " << wgt_geant << " SBband "
     << wgt_sbfit << " PionReweght " << wgt_pionReweight << " Anisodd " <<
              wgt_anisodd << " MK " << wgt_mk << "\n";
  */

  return wgt_genie * wgt_flux * wgt_2p2h * wgt_rpa * wgt_lowq2 * wgt_mueff *
         wgt_anisodd * wgt_michel * wgt_diffractive * wgt_mk * wgt_target *
         wgt_fsi * wgt_coh * wgt_geant * wgt_sbfit * wgt_pionReweight *
         wgt_tpi_warp * wgt_CCPiWegiht;
}

//==============================================================================
// Physics Calculations
//==============================================================================
TVector3 CVUniverse::AdlerAngle(int RefSystemDef, double dmumom /*GeV*/,
                                double dpimom /*GeV*/, TVector3 NeuDir,
                                TVector3 MuDir, TVector3 PiDir,
                                double Enu /*GeV*/) const {
  if (Enu < 0 || dmumom < 0 || dpimom < 0) return TVector3(-9, -9, -9);

  // double MUON_MASS = 105.658;
  // double PION_MASS = 139.570;
  // double PROTON_MASS = 938.272;
  // double bindE = 25.;

  double Epion = sqrt(CCNuPionIncConsts::CHARGED_PION_MASS *
                          CCNuPionIncConsts::CHARGED_PION_MASS +
                      dpimom * dpimom);  // ??
  double Emuon =
      sqrt(CCNuPionIncConsts::MUON_MASS * CCNuPionIncConsts::MUON_MASS +
           dmumom * dmumom);  // ??
  double Eproton =
      CCNuPionIncConsts::PROTON_MASS - CCNuPionIncConsts::bindE;  // MeV

  double Edelta = (CCNuPionIncConsts::PROTON_MASS - CCNuPionIncConsts::bindE) +
                  Enu - Emuon;  // Assume target nucleon at rest.

  double delta[3];  /// this is just Q_3
  for (int i = 0; i < 3; i++)
    delta[i] = Enu * NeuDir[i] -
               dmumom * MuDir[i];  ////vector directors must be normalized

  double beta[3];  /// proper of the boost to the Delta rest frame

  double b2 = 0.;

  for (int i = 0; i < 3; i++) {
    beta[i] = delta[i] / Edelta;
    b2 += beta[i] * beta[i];
  }

  double gamma = 1. / sqrt(1. - b2);

  double b = sqrt(b2);

  TVector3 AdlerSyst;

  if (b > 1.) {
    AdlerSyst[0] = -1000;
    AdlerSyst[1] = -1000;
    AdlerSyst[2] = -1000;
    std::cout << " |beta| should be less than one ----->    " << b
              << std::endl;  /// this canot happen by definition

    return AdlerSyst;
  }

  double piparallel[3];
  double piperpend[3];
  double piboost[3];
  double piboostabs = 0.;

  double nuparallel[3];
  double nuperpend[3];
  double nuboost[3];
  double nuboostabs = 0.;

  double muparallel[3];
  double muperpend[3];
  double muboost[3];
  double muboostabs = 0.;

  double prparallel[3];
  double prperpend[3];
  double prboost[3];
  double prboostabs = 0.;

  for (int i = 0; i < 3; i++) {
    //    piparallel[i] =
    //    beta[0]/b*dpimom*PiDir[0]+beta[1]/b*dpimom*PiDir[1]+beta[2]/b*dpimom*PiDir[2];
    //    // Pre boost
    piparallel[i] =
        (beta[i] / b2) *
        (beta[0] * PiDir[0] + beta[1] * PiDir[1] + beta[2] * PiDir[2]) * dpimom;
    piperpend[i] = dpimom * PiDir[i] - piparallel[i];
    piparallel[i] = gamma * (piparallel[i] - beta[i] * Epion);  // After boost
    piboost[i] = piparallel[i] + piperpend[i];
    piboostabs += piboost[i] * piboost[i];

    prparallel[i] = 0.;  // Pre boost
    prperpend[i] = 0.;
    prparallel[i] = gamma * (prparallel[i] - beta[i] * Eproton);  // After boost
    prboost[i] = prparallel[i] + prperpend[i];
    prboostabs += prboost[i] * prboost[i];

    // nuparallel[i] =
    // beta[i]/b*Enu/Abspv[ineut]*(beta[0]/b*Pmomv[ineut][0]+beta[1]/b*Pmomv[ineut][1]+beta[2]/b*Pmomv[ineut][2]);
    // // Pre boost
    nuparallel[i] = (beta[i] / b2) *
                    (beta[0] * Enu * NeuDir[0] + beta[1] * Enu * NeuDir[1] +
                     beta[2] * Enu * NeuDir[2]);  // Pre boost
    nuperpend[i] = Enu * NeuDir[i] - nuparallel[i];
    nuparallel[i] = gamma * (nuparallel[i] - beta[i] * Enu);  // After boost
    nuboost[i] = nuparallel[i] + nuperpend[i];
    nuboostabs += nuboost[i] * nuboost[i];

    muparallel[i] = (beta[i] / b2) *
                    (beta[0] * dmumom * MuDir[0] + beta[1] * dmumom * MuDir[1] +
                     beta[2] * dmumom * MuDir[2]);  // Pre boost
    muperpend[i] = dmumom * MuDir[i] - muparallel[i];
    muparallel[i] = gamma * (muparallel[i] - beta[i] * Emuon);  // After boost
    muboost[i] = muparallel[i] + muperpend[i];
    muboostabs += muboost[i] * muboost[i];
  }

  piboostabs = sqrt(piboostabs);
  nuboostabs = sqrt(nuboostabs);
  muboostabs = sqrt(muboostabs);

  double angtest = 0.;

  for (int i = 0; i < 3; i++)
    angtest += piboost[i] * beta[i] / b / piboostabs;  // Z = beta

  angtest = acos(angtest);

  double x[3];
  double y[3];
  double z[3];

  if (RefSystemDef ==
      0) {  // Closer to Nuclear angle defintion ? --->> Z = beta, Y = Z x mu

    z[0] = beta[0] / b;
    z[1] = beta[1] / b;
    z[2] = beta[2] / b;

    y[0] = muboost[1] * z[2] - muboost[2] * z[1];
    y[1] = muboost[2] * z[0] - muboost[0] * z[2];
    y[2] = muboost[0] * z[1] - muboost[1] * z[0];

  } else if (RefSystemDef ==
             1) {  // Radecky et al. (Wrong?) //// ---> Z = nu-mu (after boost),
                   // y = nu x mu (after boost)

    z[0] = nuboost[0] / nuboostabs - muboost[0] / muboostabs;
    z[1] = nuboost[1] / nuboostabs - muboost[1] / muboostabs;
    z[2] = nuboost[2] / nuboostabs - muboost[2] / muboostabs;

    double norm = sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);

    z[0] /= norm;
    z[1] /= norm;
    z[2] /= norm;

    y[0] = nuboost[1] * muboost[2] - nuboost[2] * muboost[1];
    y[1] = nuboost[2] * muboost[0] - nuboost[0] * muboost[2];
    y[2] = nuboost[0] * muboost[1] - nuboost[1] * muboost[0];

  } else {  // P.Allen et al. //// Z = nu-mu (after boost), Y = Z x mu (after
            // boost)

    z[0] = nuboost[0] - muboost[0];
    z[1] = nuboost[1] - muboost[1];
    z[2] = nuboost[2] - muboost[2];

    double norm = sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);

    z[0] /= norm;
    z[1] /= norm;
    z[2] /= norm;

    // std::cout << z[0] << "  " << z[1] << "   " << z[2] << std::endl;
    y[0] = z[1] * muboost[2] - z[2] * muboost[1];
    y[1] = z[2] * muboost[0] - z[0] * muboost[2];
    y[2] = z[0] * muboost[1] - z[1] * muboost[0];
  }

  double normy = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);

  for (int i = 0; i < 3; i++) y[i] = y[i] / normy;

  x[0] = y[1] * z[2] - y[2] * z[1];
  x[1] = y[2] * z[0] - y[0] * z[2];
  x[2] = y[0] * z[1] - y[1] * z[0];

  double ang = 0.;

  for (int i = 0; i < 3; i++) ang += piboost[i] * z[i] / piboostabs;

  ang = acos(ang);

  double px = piboost[0] * x[0] + piboost[1] * x[1] + piboost[2] * x[2];
  double py = piboost[0] * y[0] + piboost[1] * y[1] + piboost[2] * y[2];

  double phi = atan2(py / piboostabs, px / piboostabs);

  // if ((px < 0 && py < 0) || (px < 0 && py > 0)) phi = phi +
  // CCNuPionIncConsts::PI;
  // if (px > 0 && py < 0) phi = phi + 2*CCNuPionIncConsts::PI;

  AdlerSyst[0] = angtest;
  AdlerSyst[1] = ang;
  AdlerSyst[2] = phi;

  return AdlerSyst;
}

double CVUniverse::CalcQ2(const double Enu, const double Emu,
                          const double Thetamu) const {
  double Q2 =
      2.0 * Enu *
          (Emu - sqrt(pow(Emu, 2.0) - pow(CCNuPionIncConsts::MUON_MASS, 2.0)) *
                     cos(Thetamu)) -
      pow(CCNuPionIncConsts::MUON_MASS, 2.0);
  if (Q2 < 0.) Q2 = 0.0;
  return Q2;
}

double CVUniverse::CalcWexp(const double Q2, const double Ehad) const {
  double mass_hadron =
      (CCNuPionIncConsts::PROTON_MASS + CCNuPionIncConsts::NEUTRON_MASS) / 2;
  double W = pow(mass_hadron, 2.0) - Q2 + 2.0 * (mass_hadron)*Ehad;
  W = W > 0 ? sqrt(W) : 0.0;
  return W;
}

double CVUniverse::Calcq0(const double Enu, const double Emu) const {
  return Enu - Emu;
}

double CVUniverse::Calcq3(const double Q2, const double Enu,
                          const double Emu) const {
  return sqrt(Q2 + pow(Enu - Emu, 2.0));
}

double CVUniverse::Calct(const double pxpi, const double pypi,
                         const double pzpi, const double epi, const double pxmu,
                         const double pymu, const double pzmu,
                         const double emu) const {
  return pow((epi + emu - pzmu - pzpi), 2.0) + pow(pxpi + pxmu, 2.0) +
         pow(pypi + pymu, 2.0);
}

#endif  // CVUniverse_cxx
