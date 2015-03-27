#ifndef FinalStateAnalysis_DataFormats_LeptonPairsFwd_h
#define FinalStateAnalysis_DataFormats_LeptonPairsFwd_h
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

typedef FinalStateT<pat::Electron, pat::Electron> ElecElecFinalState;
typedef FinalStateT<pat::Electron, pat::Muon> ElecMuFinalState;
typedef FinalStateT<pat::Electron, pat::Tau> ElecTauFinalState;
typedef FinalStateT<pat::Electron, pat::Photon> ElecPhoFinalState;
typedef FinalStateT<pat::Muon, pat::Muon> MuMuFinalState;
typedef FinalStateT<pat::Muon, pat::Tau> MuTauFinalState;
typedef FinalStateT<pat::Muon, pat::Photon> MuPhoFinalState;
typedef FinalStateT<pat::Tau, pat::Tau> TauTauFinalState;
typedef FinalStateT<pat::Tau, pat::Photon> TauPhoFinalState;
typedef FinalStateT<pat::Photon, pat::Photon> PhoPhoFinalState;
typedef FinalStateT<pat::Muon, pat::Jet> MuJetFinalState;
typedef FinalStateT<pat::Electron, pat::Jet> ElecJetFinalState;

FWD_TYPEDEFS(ElecElecFinalState)
FWD_TYPEDEFS(ElecMuFinalState)
FWD_TYPEDEFS(ElecTauFinalState)
FWD_TYPEDEFS(ElecPhoFinalState)
FWD_TYPEDEFS(MuMuFinalState)
FWD_TYPEDEFS(MuTauFinalState)
FWD_TYPEDEFS(MuPhoFinalState)
FWD_TYPEDEFS(TauTauFinalState)
FWD_TYPEDEFS(TauPhoFinalState)
FWD_TYPEDEFS(PhoPhoFinalState)
FWD_TYPEDEFS(MuJetFinalState)
FWD_TYPEDEFS(ElecJetFinalState)

#endif
