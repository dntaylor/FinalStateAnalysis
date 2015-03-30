/*
 * Embed float into a PAT final state.  The float is taken from a double in the
 * event.
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

class FinalStateFloatEmbedder : public edm::EDProducer {
  public:
    FinalStateFloatEmbedder(const edm::ParameterSet& pset);
    virtual ~FinalStateFloatEmbedder(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);
  private:
    edm::InputTag src_;
    edm::InputTag toEmbedSrc_;
    std::string name_;
};

FinalStateFloatEmbedder::FinalStateFloatEmbedder(
    const edm::ParameterSet& pset) {
  src_ = pset.getParameter<edm::InputTag>("src");
  toEmbedSrc_ = pset.getParameter<edm::InputTag>("toEmbedSrc");
  name_ = pset.getParameter<std::string>("name");
  produces<FinalStateCollection>();
}
void FinalStateFloatEmbedder::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::auto_ptr<FinalStateCollection> output(new FinalStateCollection);

  edm::Handle<edm::View<FinalState> > finalStatesH;
  evt.getByLabel(src_, finalStatesH);

  edm::Handle<double> toEmbedH;
  evt.getByLabel(toEmbedSrc_, toEmbedH);

  for (size_t i = 0; i < finalStatesH->size(); ++i) {
    FinalState* embedInto = finalStatesH->ptrAt(i)->clone();
    embedInto->addUserFloat(name_, *toEmbedH);
    output->push_back(embedInto); // takes ownership
  }
  evt.put(output);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FinalStateFloatEmbedder);
