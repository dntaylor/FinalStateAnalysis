#ifndef FinalStateAnalysis_DataFormats_LeptonTripletsFwd_h
#define FinalStateAnalysis_DataFormats_LeptonTripletsFwd_h
#include "FinalStateAnalysis/DataFormats/interface/FinalStateT.h"

// Forward declarations

namespace pat {
  class Electron;
  class Muon;
  class Tau;
  class Photon;
}

#include "FinalStateAnalysis/DataFormats/interface/FwdIncludes.h"
#include "FinalStateAnalysis/DataFormats/interface/Macros.h"

typedef FinalStateT<pat::Electron, pat::Electron, pat::Electron> ElecElecElecFinalState;
FWD_TYPEDEFS(ElecElecElecFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Muon> ElecElecMuFinalState;
FWD_TYPEDEFS(ElecElecMuFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Tau> ElecElecTauFinalState;
FWD_TYPEDEFS(ElecElecTauFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Photon> ElecElecPhoFinalState;
FWD_TYPEDEFS(ElecElecPhoFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Muon> ElecMuMuFinalState;
FWD_TYPEDEFS(ElecMuMuFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Tau> ElecMuTauFinalState;
FWD_TYPEDEFS(ElecMuTauFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Photon> ElecMuPhoFinalState;
FWD_TYPEDEFS(ElecMuPhoFinalState)
typedef FinalStateT<pat::Electron, pat::Tau, pat::Tau> ElecTauTauFinalState;
FWD_TYPEDEFS(ElecTauTauFinalState)
typedef FinalStateT<pat::Electron, pat::Photon, pat::Photon> ElecPhoPhoFinalState;
FWD_TYPEDEFS(ElecPhoPhoFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Muon> MuMuMuFinalState;
FWD_TYPEDEFS(MuMuMuFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Tau> MuMuTauFinalState;
FWD_TYPEDEFS(MuMuTauFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Photon> MuMuPhoFinalState;
FWD_TYPEDEFS(MuMuPhoFinalState)
typedef FinalStateT<pat::Muon, pat::Tau, pat::Tau> MuTauTauFinalState;
FWD_TYPEDEFS(MuTauTauFinalState)
typedef FinalStateT<pat::Muon, pat::Photon, pat::Photon> MuPhoPhoFinalState;
FWD_TYPEDEFS(MuPhoPhoFinalState)
typedef FinalStateT<pat::Muon, pat::Jet, pat::Jet> MuJetJetFinalState;
FWD_TYPEDEFS(MuJetJetFinalState)
typedef FinalStateT<pat::Tau, pat::Tau, pat::Tau> TauTauTauFinalState;
FWD_TYPEDEFS(TauTauTauFinalState)
#endif
