/*
 * A workaround that just copies PaTFinalStates into a new collection.
 *
 * For process naming purposes.
 *
 * Author: Evan K. Friis, UW Madison
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include <string>
#include "FinalStateAnalysis/DataFormats/interface/FinalState.h"
#include "FinalStateAnalysis/DataFormats/interface/FinalStateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class FinalStateCopier : public edm::EDProducer {
  public:
    FinalStateCopier(const edm::ParameterSet& pset);
    virtual ~FinalStateCopier(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);
  private:
    edm::InputTag src_;
    std::string name_;
};

FinalStateCopier::FinalStateCopier(
    const edm::ParameterSet& pset) {
  src_ = pset.getParameter<edm::InputTag>("src");
  produces<FinalStateCollection>();
}
void FinalStateCopier::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::auto_ptr<FinalStateCollection> output(new FinalStateCollection);

  edm::Handle<edm::View<FinalState> > finalStatesH;
  evt.getByLabel(src_, finalStatesH);

  for (size_t i = 0; i < finalStatesH->size(); ++i) {
    FinalState* embedInto = finalStatesH->ptrAt(i)->clone();
    output->push_back(embedInto); // takes ownership
  }
  evt.put(output);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FinalStateCopier);
