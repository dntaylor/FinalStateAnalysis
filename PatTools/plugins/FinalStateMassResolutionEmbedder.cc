/*
 * Embed the the mass resolution (calculated analytically)
 * into the FinalState object. 
 * If the error is not available for daughter type the mass error is set to -1
 *
 * Author: L. Gray, FNAL
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
#include "FinalStateAnalysis/PatTools/interface/FinalStateMassResolution.h"

class FinalStateMassResolutionEmbedder : public edm::EDProducer {
public:
  FinalStateMassResolutionEmbedder(const edm::ParameterSet& pset);
  virtual ~FinalStateMassResolutionEmbedder(){ delete resoCalc_;}
  void produce(edm::Event& evt, const edm::EventSetup& es);
private:
  const bool debug_;
  edm::InputTag src_;    
  FinalStateMassResolution* resoCalc_;  
  
};

FinalStateMassResolutionEmbedder::
FinalStateMassResolutionEmbedder(const edm::ParameterSet& pset):
  debug_(pset.getParameter<bool>("debug")),
  src_(pset.getParameter<edm::InputTag>("src"))
{  
  resoCalc_ = new FinalStateMassResolution();    
  produces<FinalStateCollection>();
}
void 
FinalStateMassResolutionEmbedder::
produce(edm::Event& evt, const edm::EventSetup& es) {
  std::auto_ptr<FinalStateCollection> output(new FinalStateCollection);

  resoCalc_->init(es);

  edm::Handle<edm::View<FinalState> > finalStatesH;
  evt.getByLabel(src_, finalStatesH);
  
  for (size_t i = 0; i < finalStatesH->size(); ++i) {
    FinalState* embedInto = finalStatesH->ptrAt(i)->clone();    
    std::vector<double> components;
    double dM = 0;
    
    try {
      dM = resoCalc_->getMassResolutionWithComponents(*embedInto,
						      components);
    } catch (std::bad_cast& e) {
      dM = -1;
    }
    
    if( debug_ ) {
      std::cout << "Candidate " << i 
		<< " calculated dM = "<< dM  << " GeV" <<std::endl;
      std::cout << "\tComponents: ";
      for( size_t k = 0; k < components.size(); ++k )
	std::cout << components[k] 
		  << ( (k != components.size() - 1) ? " GeV, " : " GeV" );
      std::cout << std::endl;
    }

    embedInto->addUserFloat("cand_dM", dM);
    for( size_t k = 0; k < components.size(); ++k )
      embedInto->addUserFloat(Form("cand_dM_%lu",k), components[k]);
    
    output->push_back(embedInto); // takes ownership
  }
  evt.put(output);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FinalStateMassResolutionEmbedder);
