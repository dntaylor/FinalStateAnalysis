#include "FinalStateAnalysis/PatTools/plugins/PATQuadFinalStateBuilderHzzT.h"
#include "FinalStateAnalysis/DataFormats/interface/PATQuadLeptonFinalStates.h"

typedef PATQuadFinalStateBuilderHzzT<PATElecElecElecElecFinalState> PATElecElecElecElecFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecElecMuFinalState> PATElecElecElecMuFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecElecTauFinalState> PATElecElecElecTauFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecElecPhoFinalState> PATElecElecElecPhoFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecMuMuFinalState> PATElecElecMuMuFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecMuTauFinalState> PATElecElecMuTauFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecMuPhoFinalState> PATElecElecMuPhoFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecTauTauFinalState> PATElecElecTauTauFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecElecPhoPhoFinalState> PATElecElecPhoPhoFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecMuMuMuFinalState> PATElecMuMuMuFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecMuMuTauFinalState> PATElecMuMuTauFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecMuMuPhoFinalState> PATElecMuMuPhoFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecMuTauTauFinalState> PATElecMuTauTauFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATElecMuPhoPhoFinalState> PATElecMuPhoPhoFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATMuMuMuMuFinalState> PATMuMuMuMuFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATMuMuMuTauFinalState> PATMuMuMuTauFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATMuMuMuPhoFinalState> PATMuMuMuPhoFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATMuMuTauTauFinalState> PATMuMuTauTauFinalStateHzzProducer;
typedef PATQuadFinalStateBuilderHzzT<PATMuMuPhoPhoFinalState> PATMuMuPhoPhoFinalStateHzzProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATElecElecElecElecFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecElecMuFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecElecTauFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecElecPhoFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecMuMuFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecMuTauFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecMuPhoFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecTauTauFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecElecPhoPhoFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecMuMuMuFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecMuMuTauFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecMuMuPhoFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecMuTauTauFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATElecMuPhoPhoFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATMuMuMuMuFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATMuMuMuTauFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATMuMuMuPhoFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATMuMuTauTauFinalStateHzzProducer);
DEFINE_FWK_MODULE(PATMuMuPhoPhoFinalStateHzzProducer);
