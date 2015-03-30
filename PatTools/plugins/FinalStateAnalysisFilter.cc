/*
 * EDFilter which saves events which pass all selections of a
 * FinalStateAnalysis.
 *
 * Author: Evan K. Friis UW Madison
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/LuminosityBlockBase.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FinalStateAnalysis/NtupleTools/interface/FinalStateAnalysis.h"

class FinalStateAnalysisFilter : public edm::EDFilter {
  public:
    FinalStateAnalysisFilter(const edm::ParameterSet& pset);
    virtual ~FinalStateAnalysisFilter(){}
    void beginJob();
    void endJob();
    bool filter(edm::Event& evt, const edm::EventSetup& es);
    bool beginLuminosityBlock(
        edm::LuminosityBlock& ls, const edm::EventSetup & es);
    bool endLuminosityBlock(
        edm::LuminosityBlock& ls, const edm::EventSetup & es);
  private:
    std::auto_ptr<FinalStateAnalysis> analysis_;
};

FinalStateAnalysisFilter::FinalStateAnalysisFilter(
    const edm::ParameterSet& pset) {
  edm::Service<TFileService> fs;
  TFileDirectory &fd =  fs->tFileDirectory();
  analysis_.reset(new FinalStateAnalysis(pset, fd));
}

bool FinalStateAnalysisFilter::filter(
    edm::Event& evt, const edm::EventSetup& es) {
  const edm::EventBase& evtBase = evt;
  return analysis_->filter(evtBase);
}

bool FinalStateAnalysisFilter::beginLuminosityBlock(
    edm::LuminosityBlock& ls, const edm::EventSetup& es) {
  const edm::LuminosityBlockBase& lsBase = ls;
  analysis_->beginLuminosityBlock(lsBase);
  return true;
}

bool FinalStateAnalysisFilter::endLuminosityBlock(
    edm::LuminosityBlock& ls, const edm::EventSetup& es) {
  const edm::LuminosityBlockBase& lsBase = ls;
  analysis_->endLuminosityBlock(lsBase);
  return true;
}

void FinalStateAnalysisFilter::beginJob() {
  analysis_->beginJob();
}
void FinalStateAnalysisFilter::endJob() {
  analysis_->endJob();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FinalStateAnalysisFilter);
