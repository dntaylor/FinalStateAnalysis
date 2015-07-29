// -*- C++ -*-
//
// Package:    GenLeptonFSRFilter
// Class:      GenLeptonFSRFilter
// 
/**\class GenLeptonFSRFilter GenLeptonFSRFilter.cc FinalStateAnalysis/GenLeptonFSRFilter/src/GenLeptonFSRFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nick Smith
//         Created:  Fri Jun 12 02:15:46 CDT 2015
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class GenLeptonFSRFilter : public edm::EDFilter {
   public:
      explicit GenLeptonFSRFilter(const edm::ParameterSet&);
      ~GenLeptonFSRFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

   edm::InputTag GenParticleTag_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenLeptonFSRFilter::GenLeptonFSRFilter(const edm::ParameterSet& iConfig) :
      GenParticleTag_(iConfig.getUntrackedParameter<edm::InputTag> ("GenTag", edm::InputTag("genParticles")))
{

}


GenLeptonFSRFilter::~GenLeptonFSRFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GenLeptonFSRFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel(GenParticleTag_, genParticles);

   for ( auto& genPhoton : *genParticles ) {
      // Get the photon
      if ( genPhoton.status() == 1 
            && genPhoton.pdgId() == 22 
            && (abs(genPhoton.mother(0)->pdgId()) == 11 || abs(genPhoton.mother(0)->pdgId()) == 13 || abs(genPhoton.mother(0)->pdgId()) == 15 || abs(genPhoton.mother(0)->pdgId()) < 7) 
            //&& ((abs(genPhoton.mother(0)->pdgId()) == 11 || abs(genPhoton.mother(0)->pdgId()) == 13 || abs(genPhoton.mother(0)->pdgId()) == 15 || abs(genPhoton.mother(0)->pdgId()) < 7) 
            //|| (abs(genPhoton.mother(1)->pdgId()) == 11 || abs(genPhoton.mother(1)->pdgId()) == 13 || abs(genPhoton.mother(1)->pdgId()) == 15 || abs(genPhoton.mother(1)->pdgId()) < 7))
            && genPhoton.pt() > 10.
            ) {
         bool keepEvent = false;
         for  ( auto& genLepton : *genParticles ) {
           // remove if all leptons are far from photon
           double deltaR = reco::deltaR(genPhoton.p4(), genLepton.p4());
           if ( genLepton.status() == 1
                && (abs(genLepton.pdgId()) == 11 || abs(genLepton.pdgId()) == 13 || abs(genLepton.pdgId()) == 15)
                && deltaR < 0.4
                ) {
             keepEvent = true;
          }
          if (!keepEvent) {
             std::cout << "killing a gamma pt = " << genPhoton.pt() << std::endl;
             return false;
          }
        }
      }
   }

   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenLeptonFSRFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenLeptonFSRFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
GenLeptonFSRFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
GenLeptonFSRFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
GenLeptonFSRFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
GenLeptonFSRFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenLeptonFSRFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenLeptonFSRFilter);

