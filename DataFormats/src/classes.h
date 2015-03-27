#include "FinalStateAnalysis/DataFormats/interface/FinalState.h"
#include "FinalStateAnalysis/DataFormats/interface/FinalStateFwd.h"

#include "FinalStateAnalysis/DataFormats/interface/FinalStateProxy.h"

#include "FinalStateAnalysis/DataFormats/interface/MultiCandFinalState.h"
#include "FinalStateAnalysis/DataFormats/interface/MultiCandFinalStateFwd.h"

#include "FinalStateAnalysis/DataFormats/interface/FinalStateEvent.h"

#include "FinalStateAnalysis/DataFormats/interface/FinalStateLS.h"

#include "FinalStateAnalysis/DataFormats/interface/SingleFinalStates.h"

#include "FinalStateAnalysis/DataFormats/interface/DiLeptonFinalStates.h"

#include "FinalStateAnalysis/DataFormats/interface/TriLeptonFinalStates.h"

#include "FinalStateAnalysis/DataFormats/interface/QuadLeptonFinalStates.h"

#include "FinalStateAnalysis/DataAlgos/interface/VBFVariables.h"

#include "FinalStateAnalysis/DataFormats/interface/Macros.h"

namespace {
  struct FinalStateAnalysis_DataFormats_dicts {
    // General missing dictionaries
    edm::Ptr<reco::Vertex> dummyVertexPtr;
    edm::PtrVector<reco::Vertex> dummyVertexPtrVector;

    edm::Ptr<pat::Photon> dummyPhotonPtr;

    edm::RefProd<pat::ElectronCollection> dummyElectronRefProd;
    edm::RefProd<pat::MuonCollection> dummyMuonRefProd;
    edm::RefProd<pat::TauCollection> dummyTauRefProd;
    edm::RefProd<pat::PhotonCollection> dummyPhotonRefProd;
    edm::RefProd<pat::JetCollection> dummyJetRefProd;

    std::map<std::string, float> dummyFloatMap;
    std::map<std::string, int> dummyIntMap;
    std::pair<std::string, float> dummyFloatPair;
    std::pair<std::string, int> dummyIntPair;
    std::map<std::string, edm::Ptr<pat::MET> > dummyMETMap;

    // For the VBF variables
    VBFVariables dummyVBFVars;

    // shared pointer wrapper class
    FinalStateProxy proxyDummy;

    // base classes
    FWD_ABS_CLASSDECL(FinalState)
    FWD_CLASSDECL(FinalStateEvent)
    FWD_CLASSDECL(FinalStateLS)

    // n-cand state
    FWD_CLASSDECL(MultiCandFinalState)

    // single final states
    FWD_MIN_CLASSDECL(ElecFinalState)
    FWD_MIN_CLASSDECL(MuFinalState)
    FWD_MIN_CLASSDECL(TauFinalState)
    FWD_MIN_CLASSDECL(PhoFinalState)
    FWD_MIN_CLASSDECL(JetFinalState)

    // pair final states
    FWD_MIN_CLASSDECL(ElecElecFinalState)
    FWD_MIN_CLASSDECL(ElecMuFinalState)
    FWD_MIN_CLASSDECL(ElecTauFinalState)
    FWD_MIN_CLASSDECL(ElecPhoFinalState)
    FWD_MIN_CLASSDECL(MuMuFinalState)
    FWD_MIN_CLASSDECL(MuTauFinalState)
    FWD_MIN_CLASSDECL(MuPhoFinalState)
    FWD_MIN_CLASSDECL(TauTauFinalState)
    FWD_MIN_CLASSDECL(TauPhoFinalState)
    FWD_MIN_CLASSDECL(PhoPhoFinalState)
    FWD_MIN_CLASSDECL(MuJetFinalState)
    FWD_MIN_CLASSDECL(ElecJetFinalState)

    // triplet final states
    FWD_MIN_CLASSDECL(ElecElecElecFinalState)
    FWD_MIN_CLASSDECL(ElecElecMuFinalState)
    FWD_MIN_CLASSDECL(ElecElecTauFinalState)
    FWD_MIN_CLASSDECL(ElecElecPhoFinalState)
    FWD_MIN_CLASSDECL(ElecMuMuFinalState)
    FWD_MIN_CLASSDECL(ElecMuTauFinalState)
    FWD_MIN_CLASSDECL(ElecMuPhoFinalState)
    FWD_MIN_CLASSDECL(ElecTauTauFinalState)
    FWD_MIN_CLASSDECL(ElecPhoPhoFinalState)
    FWD_MIN_CLASSDECL(MuMuMuFinalState)
    FWD_MIN_CLASSDECL(MuMuTauFinalState)
    FWD_MIN_CLASSDECL(MuMuPhoFinalState)
    FWD_MIN_CLASSDECL(MuTauTauFinalState)
    FWD_MIN_CLASSDECL(MuPhoPhoFinalState)
    FWD_MIN_CLASSDECL(MuJetJetFinalState)
    FWD_MIN_CLASSDECL(TauTauTauFinalState)

    // quad final states
    FWD_MIN_CLASSDECL(ElecElecElecElecFinalState)
    FWD_MIN_CLASSDECL(ElecElecElecMuFinalState)
    FWD_MIN_CLASSDECL(ElecElecElecTauFinalState)
    FWD_MIN_CLASSDECL(ElecElecElecPhoFinalState)
    FWD_MIN_CLASSDECL(ElecElecMuMuFinalState)
    FWD_MIN_CLASSDECL(ElecElecMuTauFinalState)
    FWD_MIN_CLASSDECL(ElecElecMuPhoFinalState)
    FWD_MIN_CLASSDECL(ElecElecTauTauFinalState)
    FWD_MIN_CLASSDECL(ElecElecPhoPhoFinalState)
    FWD_MIN_CLASSDECL(ElecMuMuMuFinalState)
    FWD_MIN_CLASSDECL(ElecMuMuTauFinalState)
    FWD_MIN_CLASSDECL(ElecMuMuPhoFinalState)
    FWD_MIN_CLASSDECL(ElecMuTauTauFinalState)
    FWD_MIN_CLASSDECL(ElecTauTauTauFinalState)
    FWD_MIN_CLASSDECL(ElecMuPhoPhoFinalState)
    FWD_MIN_CLASSDECL(MuMuMuMuFinalState)
    FWD_MIN_CLASSDECL(MuMuMuTauFinalState)
    FWD_MIN_CLASSDECL(MuMuMuPhoFinalState)
    FWD_MIN_CLASSDECL(MuMuTauTauFinalState)
    FWD_MIN_CLASSDECL(MuTauTauTauFinalState)
    FWD_MIN_CLASSDECL(TauTauTauTauFinalState)
    FWD_MIN_CLASSDECL(MuMuPhoPhoFinalState)

  };
}
