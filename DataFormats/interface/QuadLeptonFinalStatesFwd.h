#ifndef FinalStateAnalysis_DataFormats_LeptonQuadsFwd_h
#define FinalStateAnalysis_DataFormats_LeptonQuadsFwd_h
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

/*

   Use LeptonFinalStatesFwd_gen.py to generate statements.

*/

typedef FinalStateT<pat::Electron, pat::Electron, pat::Electron, pat::Electron> ElecElecElecElecFinalState;
FWD_TYPEDEFS(ElecElecElecElecFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Electron, pat::Muon> ElecElecElecMuFinalState;
FWD_TYPEDEFS(ElecElecElecMuFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Electron, pat::Tau> ElecElecElecTauFinalState;
FWD_TYPEDEFS(ElecElecElecTauFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Electron, pat::Photon> ElecElecElecPhoFinalState;
FWD_TYPEDEFS(ElecElecElecPhoFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Muon, pat::Muon> ElecElecMuMuFinalState;
FWD_TYPEDEFS(ElecElecMuMuFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Muon, pat::Tau> ElecElecMuTauFinalState;
FWD_TYPEDEFS(ElecElecMuTauFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Muon, pat::Photon> ElecElecMuPhoFinalState;
FWD_TYPEDEFS(ElecElecMuPhoFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Tau, pat::Tau> ElecElecTauTauFinalState;
FWD_TYPEDEFS(ElecElecTauTauFinalState)
typedef FinalStateT<pat::Electron, pat::Electron, pat::Photon, pat::Photon> ElecElecPhoPhoFinalState;
FWD_TYPEDEFS(ElecElecPhoPhoFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Muon, pat::Muon> ElecMuMuMuFinalState;
FWD_TYPEDEFS(ElecMuMuMuFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Muon, pat::Tau> ElecMuMuTauFinalState;
FWD_TYPEDEFS(ElecMuMuTauFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Muon, pat::Photon> ElecMuMuPhoFinalState;
FWD_TYPEDEFS(ElecMuMuPhoFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Tau, pat::Tau> ElecMuTauTauFinalState;
FWD_TYPEDEFS(ElecMuTauTauFinalState)
typedef FinalStateT<pat::Electron, pat::Muon, pat::Photon, pat::Photon> ElecMuPhoPhoFinalState;
FWD_TYPEDEFS(ElecMuPhoPhoFinalState)
typedef FinalStateT<pat::Electron, pat::Tau, pat::Tau, pat::Tau> ElecTauTauTauFinalState;
FWD_TYPEDEFS(ElecTauTauTauFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Muon, pat::Muon> MuMuMuMuFinalState;
FWD_TYPEDEFS(MuMuMuMuFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Muon, pat::Tau> MuMuMuTauFinalState;
FWD_TYPEDEFS(MuMuMuTauFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Muon, pat::Photon> MuMuMuPhoFinalState;
FWD_TYPEDEFS(MuMuMuPhoFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Tau, pat::Tau> MuMuTauTauFinalState;
FWD_TYPEDEFS(MuMuTauTauFinalState)
typedef FinalStateT<pat::Muon, pat::Muon, pat::Photon, pat::Photon> MuMuPhoPhoFinalState;
FWD_TYPEDEFS(MuMuPhoPhoFinalState)
typedef FinalStateT<pat::Muon, pat::Tau, pat::Tau, pat::Tau> MuTauTauTauFinalState;
FWD_TYPEDEFS(MuTauTauTauFinalState)
typedef FinalStateT<pat::Tau, pat::Tau, pat::Tau, pat::Tau> TauTauTauTauFinalState;
FWD_TYPEDEFS(TauTauTauTauFinalState)

#endif
