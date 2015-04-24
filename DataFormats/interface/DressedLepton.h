#ifndef DressedLepton_H
#define DressedLepton_H

#include "DataFormats/TrackReco/interface/Track.h"
#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/LorentzVector.h"

namespace reco {
  class DressedLepton : public GenParticle {
    public:
      DressedLepton(){}
      DressedLepton(const GenParticle& gen) :
        GenParticle(gen){ }
  
      void setOriginalPt(double val){originalPt_ = val;}
      void setType(int val){Type_ = val;}
      void setMotherPDGID(int val){MotherPDGID_ = val;}
      double originalPt() const {return originalPt_;}
      int Type() const {return Type_;} // 0=not filled/error 1= Good! 2=Decay in Flight Pion/Kaon/Etc 3=Tau 4=Decay from D meson 5=Decay from B hadron
      int MotherPDGID() const {return MotherPDGID_;}
      // Add some methods we assume from PATObject
      const CandidatePtrVector & overlaps(const std::string &label) const {
        static const CandidatePtrVector EMPTY;
        return EMPTY;
      }
      CandidatePtr userCand( const std::string & key ) const { return CandidatePtr(); }
  
    private:
      double originalPt_;
      int Type_;
      int MotherPDGID_;
  };

  typedef std::vector<DressedLepton> DressedLeptonCollection;
  typedef edm::Ref<DressedLeptonCollection> DressedLeptonRef;
  typedef edm::RefProd<DressedLeptonCollection> DressedLeptonRefProd;
  typedef edm::RefVector<DressedLeptonCollection> DressedLeptonRefVector;
}

#endif
