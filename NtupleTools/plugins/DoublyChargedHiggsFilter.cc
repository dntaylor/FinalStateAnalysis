// Basic (& silly) program to produce a small genlevel ntuple, and a bunch of gen level plots
// Useful for debugging
// FSA independent

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
class DoublyChargedHiggsFilter : public edm::EDFilter {

 public:
  DoublyChargedHiggsFilter (const edm::ParameterSet &);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
 private:
  edm::InputTag GenParticleTag_;
  bool selectPairProduction_, selectAssociatedProduction_;
  int decaysToEE_, decaysToMM_,decaysToTT_,decaysToEM_,decaysToET_,decaysToMT_;
  bool verbose_;

  ULong_t countPP, countAP;

  ULong_t countEEEE,countMMEE,countTTEE,countEMEE,countETEE,countMTEE;
  ULong_t countEEMM,countMMMM,countTTMM,countEMMM,countETMM,countMTMM;
  ULong_t countEETT,countMMTT,countTTTT,countEMTT,countETTT,countMTTT;
  ULong_t countEEEM,countMMEM,countTTEM,countEMEM,countETEM,countMTEM;
  ULong_t countEEET,countMMET,countTTET,countEMET,countETET,countMTET;
  ULong_t countEEMT,countMMMT,countTTMT,countEMMT,countETMT,countMTMT;
  ULong_t countEEE,countMME,countTTE,countEME,countETE,countMTE;
  ULong_t countEEM,countMMM,countTTM,countEMM,countETM,countMTM;
  ULong_t countEET,countMMT,countTTT,countEMT,countETT,countMTT;
  

  double nall;
  double nsel;
};
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <map>
#include <memory>

using namespace edm;
using namespace std;
using namespace reco;


DoublyChargedHiggsFilter::DoublyChargedHiggsFilter( const ParameterSet & cfg ) :
 GenParticleTag_(cfg.getUntrackedParameter<edm::InputTag> ("GenTag", edm::InputTag("genParticles"))),
 selectPairProduction_(cfg.getParameter<bool>("selectPairProduction")),
 selectAssociatedProduction_(cfg.getParameter<bool>("selectAssociatedProduction")),
 decaysToEE_(cfg.getUntrackedParameter<int>("decaysToEE",-1)),
 decaysToMM_(cfg.getUntrackedParameter<int>("decaysToMM",-1)),
 decaysToTT_(cfg.getUntrackedParameter<int>("decaysToTT",-1)),
 decaysToEM_(cfg.getUntrackedParameter<int>("decaysToEM",-1)),
 decaysToET_(cfg.getUntrackedParameter<int>("decaysToET",-1)),
 decaysToMT_(cfg.getUntrackedParameter<int>("decaysToMT",-1)),
 verbose_(cfg.getUntrackedParameter<bool>("verbose",false))
{
}

void DoublyChargedHiggsFilter::beginJob() {
 nall=0;
 nsel=0;	
 countPP=0;
 countAP=0;

 countEEEE=0,countMMEE=0,countTTEE=0,countEMEE=0,countETEE=0,countMTEE=0;
 countEEMM=0,countMMMM=0,countTTMM=0,countEMMM=0,countETMM=0,countMTMM=0;
 countEETT=0,countMMTT=0,countTTTT=0,countEMTT=0,countETTT=0,countMTTT=0;
 countEEEM=0,countMMEM=0,countTTEM=0,countEMEM=0,countETEM=0,countMTEM=0;
 countEEET=0,countMMET=0,countTTET=0,countEMET=0,countETET=0,countMTET=0;
 countEEMT=0,countMMMT=0,countTTMT=0,countEMMT=0,countETMT=0,countMTMT=0;
 countEEE=0,countEEM=0,countEET=0;
 countEME=0,countEMM=0,countEMT=0;
 countETE=0,countETM=0,countETT=0;
 countMME=0,countMMM=0,countMMT=0;
 countMTE=0,countMTM=0,countMTT=0;
 countTTE=0,countTTM=0,countTTT=0;

 if(selectPairProduction_&&selectAssociatedProduction_) LogError("")<<"This wont work (you cannot select both pair&associated production), check your input parameters!";
}

void DoublyChargedHiggsFilter::endJob() {
 cout<<"********************************************************************"<<endl;
 cout<<"GEN LEVEL FILTERING"<<endl<<endl;
 cout<<" Selecting PairProduction? "<<selectPairProduction_<<endl;
 cout<<" Selecting AssociatedProduction? "<<selectAssociatedProduction_<<endl;

 cout<<"Total Analyzed =   "<<nall<<endl;
 cout<<"Pair Production=   "<<countPP<<endl;
 cout<<"Associated Production=   "<<countAP<<endl;
 cout<<" -  | EE |  MM | TT | EM | ET | MT "<<endl;
 cout<<" EE "<<countEEEE<<"   "<<countMMEE<<"   "<<countTTEE<<"   "<<countEMEE<<"   "<<countETEE<<"   "<<countMTEE<<endl;
 cout<<" MM "<<countEEMM<<"   "<<countMMMM<<"   "<<countTTMM<<"   "<<countEMMM<<"   "<<countETMM<<"   "<<countMTMM<<endl;
 cout<<" TT "<<countEETT<<"   "<<countMMTT<<"   "<<countTTTT<<"   "<<countEMTT<<"   "<<countETTT<<"   "<<countMTTT<<endl; 
 cout<<" EM "<<countEEEM<<"   "<<countMMEM<<"   "<<countTTEM<<"   "<<countEMEM<<"   "<<countETEM<<"   "<<countMTEM<<endl;
 cout<<" ET "<<countEEET<<"   "<<countMMET<<"   "<<countTTET<<"   "<<countEMET<<"   "<<countETET<<"   "<<countMTET<<endl;  
 cout<<" MT "<<countEEMT<<"   "<<countMMMT<<"   "<<countTTMT<<"   "<<countEMMT<<"   "<<countETMT<<"   "<<countMTMT<<endl;
 cout<<" -  | EE |  MM | TT | EM | ET | MT "<<endl;
 cout<<" E  "<<countEEE<<"   "<<countMME<<"   "<<countTTE<<"   "<<countEME<<"   "<<countETE<<"   "<<countMTE<<endl;
 cout<<" M  "<<countEEM<<"   "<<countMMM<<"   "<<countTTM<<"   "<<countEMM<<"   "<<countETM<<"   "<<countMTM<<endl;
 cout<<" T  "<<countEET<<"   "<<countMMT<<"   "<<countTTT<<"   "<<countEMT<<"   "<<countETT<<"   "<<countMTT<<endl;
 cout<<" L  "<<countEEE+countEEM+countEET
     <<"   "<<countMME+countMMM+countMMT
     <<"   "<<countTTE+countTTM+countTTT
     <<"   "<<countEME+countEMM+countEMT
     <<"   "<<countETE+countETM+countETT
     <<"   "<<countMTE+countMTM+countMTT<<endl;

 cout<<"Selection  =   "<<nsel<<endl;
 cout<<"********************************************************************"<<endl;




}

bool DoublyChargedHiggsFilter::filter (Event & ev, const EventSetup &) {
 nall++;

 bool found=true;

 edm::Handle< vector<reco::GenParticle> >pGenPart;
 if(!ev.getByLabel(GenParticleTag_, pGenPart)) {return false; LogError("")<<"Collection not found!";}

 int countHPP=0, countHMM=0, countHM=0, countHP=0;

 bool pairProduction=false, associatedProduction=0; 

 int hpp_ee=0, hpp_mm=0, hpp_tt=0, hpp_em=0, hpp_et=0, hpp_mt=0;
 int hmm_ee=0, hmm_mm=0, hmm_tt=0, hmm_em=0, hmm_et=0, hmm_mt=0;
 int hp_e=0, hp_m=0, hp_t=0;
 int hm_e=0, hm_m=0, hm_t=0;

 for( size_t i = 0; i < pGenPart->size(); ++ i ) {
  const reco::GenParticle& genpart = (*pGenPart)[i];
  if(genpart.status()!=3) continue;
  if(abs(genpart.pdgId())!=9900041 && abs(genpart.pdgId())!=37) continue;
  if(verbose_) cout <<"\t"<<i<<"   "<<genpart.pdgId()<<"   "<<genpart.status()<<"    "<<genpart.pt()<<"   "<<genpart.numberOfDaughters()<<"   "<<genpart.mass()<<"   "<<genpart.charge()<<endl;

  int electron=0, muon=0, tau=0;
  for (size_t j=0; j<genpart.numberOfDaughters(); j++){
   const reco::Candidate* higgsdaughter=genpart.daughter(j);
   if(genpart.status()!=3) continue;
   if(verbose_) cout <<"\t\t daughter"<<j<<"   "<<higgsdaughter->pdgId()<<"   "<<higgsdaughter->status()<<"    "<<higgsdaughter->pt()<<"    "<<higgsdaughter->numberOfDaughters()<<endl;
   if(abs(higgsdaughter->pdgId())==11) electron++;
   if(abs(higgsdaughter->pdgId())==13) muon++;
   if(abs(higgsdaughter->pdgId())==15) tau++;
  }

  if(genpart.pdgId()==9900041) {
   countHPP++; 
   if(electron==2) hpp_ee++;
   else if(muon==2) hpp_mm++;
   else if(tau==2) hpp_tt++;
   else if(electron==1&&muon==1) hpp_em++;
   else if(electron==1&&tau==1) hpp_et++;
   else if(muon==1&&tau==1) hpp_mt++;
   else {cout<<"What is this???"<<endl;}
  }
  if(genpart.pdgId()==-9900041) {
   countHMM++; 
   if(electron==2) hmm_ee++;
   else if(muon==2) hmm_mm++;
   else if(tau==2) hmm_tt++;
   else if(electron==1&&muon==1) hmm_em++;
   else if(electron==1&&tau==1) hmm_et++;
   else if(muon==1&&tau==1) hmm_mt++;
   else {cout<<"What is this???"<<endl;}
  }
  if(genpart.pdgId()==37) {
   countHP++;
   if(electron==1) hp_e++;
   else if(muon==1) hp_m++;
   else if(tau==1) hp_t++;
   else {cout<<"What is this???"<<endl;}
  }
  if(genpart.pdgId()==-37) {
   countHM++;
   if(electron==1) hm_e++;
   else if(muon==1) hm_m++;
   else if(tau==1) hm_t++;
   else {cout<<"What is this???"<<endl;}
  }
 }
 if(countHPP==1 && countHMM==1) {pairProduction=true; countPP++;}
 if( (countHPP+countHMM)==1 && (countHM+countHP)==1 ) {associatedProduction=true; countAP++;}

 if(selectPairProduction_&& !pairProduction) found=false;
 if(selectAssociatedProduction_&& !associatedProduction) found=false;

 if(decaysToEE_!=-1 &&(hpp_ee+hmm_ee)!=decaysToEE_) found=false;
 if(decaysToMM_!=-1 &&(hpp_mm+hmm_mm)!=decaysToMM_) found=false;
 if(decaysToTT_!=-1 &&(hpp_tt+hmm_tt)!=decaysToTT_) found=false;
 if(decaysToEM_!=-1 &&(hpp_em+hmm_em)!=decaysToEM_) found=false;
 if(decaysToET_!=-1 &&(hpp_et+hmm_et)!=decaysToET_) found=false;
 if(decaysToMT_!=-1 &&(hpp_mt+hmm_mt)!=decaysToMT_) found=false;


 // This is a bit ugly, sorry - I'm too sleepy to code it nicely in a matrix
 if(hpp_ee==1&&hmm_ee==1) countEEEE++;  if(hpp_ee==1&&hmm_mm==1) countEEMM++; if(hpp_ee==1&&hmm_tt==1) countEETT++; 
 if(hpp_mm==1&&hmm_ee==1) countMMEE++;  if(hpp_mm==1&&hmm_mm==1) countMMMM++; if(hpp_mm==1&&hmm_tt==1) countMMTT++;
 if(hpp_tt==1&&hmm_ee==1) countTTEE++;  if(hpp_tt==1&&hmm_mm==1) countTTMM++; if(hpp_tt==1&&hmm_tt==1) countTTTT++;
 if(hpp_em==1&&hmm_ee==1) countEMEE++;  if(hpp_em==1&&hmm_mm==1) countEMMM++; if(hpp_em==1&&hmm_tt==1) countEMTT++;
 if(hpp_et==1&&hmm_ee==1) countETEE++;  if(hpp_et==1&&hmm_mm==1) countETMM++; if(hpp_et==1&&hmm_tt==1) countETTT++;
 if(hpp_mt==1&&hmm_ee==1) countMTEE++;  if(hpp_mt==1&&hmm_mm==1) countMTMM++; if(hpp_mt==1&&hmm_tt==1) countMTTT++;

 if(hpp_ee==1&&hmm_em==1) countEEEM++;  if(hpp_ee==1&&hmm_et==1) countEEET++; if(hpp_ee==1&&hmm_mt==1) countEEMT++; 
 if(hpp_mm==1&&hmm_em==1) countMMEM++;  if(hpp_mm==1&&hmm_et==1) countMMET++; if(hpp_mm==1&&hmm_mt==1) countMMMT++;
 if(hpp_tt==1&&hmm_em==1) countTTEM++;  if(hpp_tt==1&&hmm_et==1) countTTET++; if(hpp_tt==1&&hmm_mt==1) countTTMT++;
 if(hpp_em==1&&hmm_em==1) countEMEM++;  if(hpp_em==1&&hmm_et==1) countEMET++; if(hpp_em==1&&hmm_mt==1) countEMMT++;
 if(hpp_et==1&&hmm_em==1) countETEM++;  if(hpp_et==1&&hmm_et==1) countETET++; if(hpp_et==1&&hmm_mt==1) countETMT++;
 if(hpp_mt==1&&hmm_em==1) countMTEM++;  if(hpp_mt==1&&hmm_et==1) countMTET++; if(hpp_mt==1&&hmm_mt==1) countMTMT++;

 if((hpp_ee+hmm_ee)==1&&(hp_e+hm_e)==1) countEEE++; if((hpp_ee+hmm_ee)==1&&(hp_m+hm_m)==1) countEEM++; if((hpp_ee+hmm_ee)==1&&(hp_t+hm_t)==1) countEET++;
 if((hpp_mm+hmm_mm)==1&&(hp_e+hm_e)==1) countMME++; if((hpp_mm+hmm_mm)==1&&(hp_m+hm_m)==1) countMMM++; if((hpp_mm+hmm_mm)==1&&(hp_t+hm_t)==1) countMMT++;
 if((hpp_tt+hmm_tt)==1&&(hp_e+hm_e)==1) countTTE++; if((hpp_tt+hmm_tt)==1&&(hp_m+hm_m)==1) countTTM++; if((hpp_tt+hmm_tt)==1&&(hp_t+hm_t)==1) countTTT++;
 if((hpp_em+hmm_em)==1&&(hp_e+hm_e)==1) countEME++; if((hpp_em+hmm_em)==1&&(hp_m+hm_m)==1) countEMM++; if((hpp_em+hmm_em)==1&&(hp_t+hm_t)==1) countEMT++;
 if((hpp_et+hmm_et)==1&&(hp_e+hm_e)==1) countETE++; if((hpp_et+hmm_et)==1&&(hp_m+hm_m)==1) countETM++; if((hpp_et+hmm_et)==1&&(hp_t+hm_t)==1) countETT++;
 if((hpp_mt+hmm_mt)==1&&(hp_e+hm_e)==1) countMTE++; if((hpp_mt+hmm_mt)==1&&(hp_m+hm_m)==1) countMTM++; if((hpp_mt+hmm_mt)==1&&(hp_t+hm_t)==1) countMTT++;

 if (found) nsel++;
 return found;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DoublyChargedHiggsFilter);
