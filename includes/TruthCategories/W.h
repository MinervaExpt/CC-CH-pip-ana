#ifndef W_H
#define W_H

#include "../CVUniverse.h"
#include "../SignalDefinition.h"

enum WType
{
  kHighW, kLowW, kNWTypes
};

WType GetTruthWType(const CVUniverse& universe, SignalDefinition signal_definition) {
  return universe.GetWexpTrue() < signal_definition.m_w_max ? kLowW : kHighW;
}

std::string GetTruthClassification_LegendLabel(WType w_category){
  switch(w_category){
    case kLowW:
      return "Low W_{exp}";
    case kHighW:
      return "High W_{exp}";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(WType w_category){
  switch(w_category){
    case kLowW:
      return "TrueLowW";
    case kHighW:
      return "TrueHighW";
    default:
      return "ERROR";
  }
}

#endif
