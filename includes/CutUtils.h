#ifndef CutUtils_h
#define CutUtils_h

#include "Constants.h"  // enum ECuts, CCNuPionIncConsts

// Analysis Cuts - default vector
const std::vector<ECuts> kCutsVector = {kNoCuts,
                                        kPrecuts,
                                        kVtx,
                                        kMinosMuon,
//  					kHelicity, //Aaron cut
//                                        khasIntVtx,//AaronCut
//                                        kMultiplicityCut,//AaronCut
//                                        kExitingMuon,//Aaron Cut
//                                        kBrokenRockMuonCut,// Aaron Cut
//                                        kHadronContainment,// Aaron Cut
                                        kAtLeastOnePionCandidateTrack,
                                        kAtLeastOneMichel,
//  					kGoodMomentum, // Aarons Cut
                                        kLLR,
                                        kNode,
                                        kWexp,
                                        kIsoProngs,
                                        kPionMult,
                                        kThetamu, 
                                        kPmu};

const std::vector<ECuts> kDefCutsVector = {
    kNoCuts, kPrecuts, kVtx, kMinosMuon, kPmu, kThetamu, kIsoProngs};

const std::vector<ECuts> kTrackedCutsVector = {//kHelicity,
 	//				       khasIntVtx,
	//				       kMultiplicityCut,
        //                                       kExitingMuon,
  	//				       kBrokenRockMuonCut,
  	//				       kHadronContainment,
					       kAtLeastOnePionCandidateTrack,
                                               kAtLeastOneMichel,
  	//				       kGoodMomentum,
                                               kLLR,
                                               kNode,
                                               kPionMult,
                                               kWexp};

const std::vector<ECuts> kUntrackedCutsVector = {//kOneMichel,
    kHasMichel, kBestMichelDistance, kClosestMichel, kOneMichel,
    kTpi, kPTmu,       kUntrackedWexp};

const std::vector<ECuts> kHasPionCut = {kHasPion};
// Remove W cut from cuts vector
const std::vector<ECuts> GetWSidebandCuts() {
  std::vector<ECuts> w_sideband_cuts = kCutsVector;
  w_sideband_cuts.erase(
      std::remove(w_sideband_cuts.begin(), w_sideband_cuts.end(), kWexp),
      w_sideband_cuts.end());
  return w_sideband_cuts;
}

// Gaudi tool cuts - only work when checking truth tuple
bool IsPrecut(ECuts c) {
  if (c == kNoCuts || c == kGoodObjects || c == kGoodVertex ||
      c == kFiducialVolume || c == kMinosActivity || c == kPrecuts)
    return true;
  else
    return false;
}

// Cut Names
std::string GetCutName(ECuts cut) {
  switch (cut) {
    case kNoCuts:
      return "No Cuts";

    case kGoodObjects:
      return "Good Objects";

    case kGoodVertex:
      return "Good Vertex";

    case kFiducialVolume:
      return "Fiducial Volume";

    case kMinosActivity:
      return "MINOS Activity";

    case kPrecuts:
      return "Anatool Precuts";

    case kVtx:
      return "vertex position Cut";

    case kMinosMatch:
      return "MINOS Muon";

    case kMinosCharge:
      return "MINOS Charge";

    case kMinosMuon:
      return "MINOS Muon";

    case kWexp:
      return "$Tracked W_{exp}$";

    case kUntrackedWexp:
      return "$Untracked W_{exp}$";

    case kIsoProngs:
      return "$<$2 Isolated Prongs";

    case kNPionCandidates:
      return "$\\pi$ candidate";

    case kAtLeastOneMichel:
      return "$>$= 1 Michel";

    case kAtLeastOnePionCandidateTrack:
      return "$>$= 1 Hadron Track";

    case kNode:
      return "Node";

    case kPionMult:
      return "Pion Multiplicity";

    case kLLR:
      return "LLR PID";

    case kAllCuts:
      return "Total";

    case kTrackQuality:
      return "General Track Quality";

    case kPmu:
      return "1.5 GeV $<$ Pmu $<$ 20 GeV";

    case kPTmu:
      return "Ptmu $<$ 1.8 GeV";

    case kThetamu:
      return "$\\theta_{\\mu}$ $<$ 20 degrees";

    case kAtLeastOnePionCandidate:
      return "At Least One Pion";

    case kHasMichel:
      return "Untracked Has Michel";

    case kBestMichelDistance:
      return "Best Michel Distance";

    case kClosestMichel:
      return "Closest Michel";

    case kOneMichel:
      return "One michel";

    case kUntrackedMpi:
      return "Michel idx";

    case kTrackedMpi:
      return "Michel idx";

    case kHasPion:
      return "Has Pion";

    case kTpi:
      return "$T_\\pi<$ 350 MeV";

    case kHelicity:
      return "Helicity";
 
    case khasIntVtx:
      return "Has Interaction Vertex";
    case kMultiplicityCut:
      return "Multiplicity";
    case kExitingMuon:
      return "Exiting Muon";
    case kBrokenRockMuonCut:
      return "Rock muon cut";
    case kHadronContainment:
      return "Hadron Containment";
    case kGoodMomentum:
      return "Good Momentum";

    default:
      std::cout << "ERROR: GetCutName unknown cut!" << std::endl;
      return "";
  };
}

#endif
