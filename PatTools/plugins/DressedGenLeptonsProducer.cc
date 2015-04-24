#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FinalStateAnalysis/DataFormats/interface/DressedLepton.h"
#include "DataFormats/Math/interface/deltaPhi.h"

class DressedGenLeptonsProducer : public edm::EDProducer {
  public:
    DressedGenLeptonsProducer(const edm::ParameterSet& pset);
    virtual ~DressedGenLeptonsProducer(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);

  private:
    edm::InputTag src_;
    int genPDGID_;
    int genSTATUS_;
    double minPT_;
    double maxETA_;
};

using namespace edm;
using namespace std;
using namespace reco;

DressedGenLeptonsProducer::DressedGenLeptonsProducer(
 const edm::ParameterSet& pset):
 src_(pset.getParameter<edm::InputTag>("src")),
 genPDGID_(pset.getUntrackedParameter<int>("genPDGID",0)),
 genSTATUS_(pset.getUntrackedParameter<int>("genSTATUS",1)),
 minPT_(pset.getUntrackedParameter<double>("minPT",0)),
 maxETA_(pset.getUntrackedParameter<double>("maxETA",100))
{
  produces<reco::DressedLeptonCollection>();
}

void DressedGenLeptonsProducer::produce(edm::Event& evt, const edm::EventSetup& es) {

  Handle<GenParticleCollection> GPCollection;
  evt.getByLabel("genParticles", GPCollection);
  const GenParticleCollection & pGenPart = *(GPCollection.product());
 
  auto_ptr<reco::DressedLeptonCollection> out(new reco::DressedLeptonCollection);
 
  for( size_t i = 0; i < pGenPart.size(); ++ i ) {
    reco::DressedLepton genpart(pGenPart[i]);
    if(genPDGID_!=0 && fabs(genpart.pdgId())!=genPDGID_) continue;
    else if (fabs(genpart.pdgId()) != 11 && fabs(genpart.pdgId())!=13) continue;
    if(genpart.status()!=genSTATUS_)continue;
    double originalPt=genpart.pt();
    double originalPx=genpart.px();
    double originalPy=genpart.py();
    double originalPz=genpart.pz();
    double dressedPx=originalPx;
    double dressedPy=originalPy;
    double dressedPz=originalPz;
 
    int particleID=genpart.pdgId();
 
    for (size_t j = 0; j < pGenPart.size(); ++ j ) {
      const reco::GenParticle photoncand=pGenPart[j];
      if (photoncand.status()!=1) continue;
      if (photoncand.pdgId()!=22) continue;
      double dPhi=deltaPhi(photoncand.phi(),genpart.phi());
      double dEta=fabs(photoncand.eta()-genpart.eta());
      double dR=sqrt(dPhi*dPhi+dEta*dEta);
      if(dR>0.1) continue;
      //cout<<"found photon!"<<endl;
      dressedPx+=photoncand.px();
      dressedPy+=photoncand.py();
      dressedPz+=photoncand.pz();
    }
 
    double mass=genpart.mass();
 
    math::XYZTLorentzVector newmomentum(dressedPx,dressedPy,dressedPz,mass);
 
    if(newmomentum.pt()<minPT_ || fabs(newmomentum.eta())>maxETA_) continue;
 
 
    genpart.setP4(newmomentum);
    genpart.setOriginalPt(originalPt);
 
    int type=0;
    int momPDGID=0;
    if(genpart.numberOfMothers()>0){
      reco::Candidate *mother=(reco::Candidate*)genpart.mother();
   
      momPDGID=mother->pdgId();
   
      if(momPDGID==particleID && mother->status()!=3)
      do{
      reco::Candidate* motherB=(reco::Candidate*)mother->mother();
      mother=motherB;
      momPDGID=mother->pdgId();
      }
      while(momPDGID==particleID && mother->status()!=3);
   
      if (abs(momPDGID)==15) { type=3;}
      else if (int(abs(momPDGID)/1000)==5 || int(abs(momPDGID)/100)==5) { type=5;}
      else if (int(abs(momPDGID)/1000)==4 || int(abs(momPDGID)/100)==4) { type=4;}
      else if (int(abs(momPDGID)/100)==1 || int(abs(momPDGID)/100)==2 || int(abs(momPDGID)/100)==3) { type=2;}
      else { type=1; }
      if(momPDGID!=particleID && abs(momPDGID)!=24 && abs(momPDGID)!=23 && abs(momPDGID)!=22){
        cout<<"WARNING"<<endl;
        cout<<"This should have been prompt but it isnt!!! We should not be here"<<endl;
        cout<<"\t D: "<<genpart.pdgId()<<" "<<genpart.status()<<endl;
        cout<<"\t M: "<<mother->pdgId()<<" "<<mother->numberOfMothers()<<" "<<mother->status()<<endl;
      }
    }
 
    genpart.setType(type);
    genpart.setMotherPDGID(momPDGID);
 
    out->push_back(genpart);
  }

  evt.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DressedGenLeptonsProducer);

