#ifndef FinalStateAnalysis_DataFormats_SingleFwd_h
#define FinalStateAnalysis_DataFormats_SingleFwd_h
#include "FinalStateAnalysis/DataFormats/interface/FinalStateT.h"
#include "FinalStateAnalysis/DataFormats/interface/FwdIncludes.h"
#include "FinalStateAnalysis/DataFormats/interface/Macros.h"

// Forward declarations

namespace pat {
  class Electron;
  class Muon;
  class Tau;
  class Photon;
  class Jet;
}

typedef FinalStateT<pat::Electron> ElecFinalState;
typedef FinalStateT<pat::Muon> MuFinalState;
typedef FinalStateT<pat::Tau> TauFinalState;
typedef FinalStateT<pat::Photon> PhoFinalState;
typedef FinalStateT<pat::Jet> JetFinalState;

FWD_TYPEDEFS(ElecFinalState)
FWD_TYPEDEFS(MuFinalState)
FWD_TYPEDEFS(TauFinalState)
FWD_TYPEDEFS(PhoFinalState)
FWD_TYPEDEFS(JetFinalState)

#endif
