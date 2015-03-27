#include "FinalStateAnalysis/DataFormats/interface/FinalStateProxy.h"
#include "FinalStateAnalysis/DataFormats/interface/FinalState.h"

FinalStateProxy::FinalStateProxy(FinalState* finalState):
  finalState_(finalState) {}

FinalStateProxy::FinalStateProxy() {}

const FinalState* FinalStateProxy::get() const {
  return finalState_.get();
}

const FinalState* FinalStateProxy::operator->() const {
  return finalState_.get();
}
