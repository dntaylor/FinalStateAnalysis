#include "FinalStateAnalysis/DataFormats/interface/FinalStateLS.h"
#include "FinalStateAnalysis/DataAlgos/interface/SmartTrigger.h"

FinalStateLS::FinalStateLS() {}

FinalStateLS::FinalStateLS(const edm::LuminosityBlockID& id,
        double integratedLuminosity,
        double instantaneousLumi):
  id_(id),
  integratedLumi_(integratedLuminosity),
  instaneousLumi_(instantaneousLumi) {}

const edm::LuminosityBlockID&
FinalStateLS::lsID() const { return id_; }

double FinalStateLS::instantaneousLumi() const {
  return instaneousLumi_;
}

double FinalStateLS::intLumi() const {
  return integratedLumi_;
}
