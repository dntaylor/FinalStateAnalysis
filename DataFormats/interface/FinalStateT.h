#ifndef FinalStateAnalysis_DataFormats_FinalStateT_h
#define FinalStateAnalysis_DataFormats_FinalStateT_h

#include "FinalStateAnalysis/DataFormats/interface/FinalState.h"

template<>
class FinalStateT : public FinalState {
  public:
    FinalStateT():FinalState(){}

    FinalStateT(const edm::Ptr<FinalStateEvent>& evt)
      :FinalState(0,NULL,evt) {
      }

    virtual FinalStateT<>* clone() const {
      return new FinalStateT<>(*this);
    }

    virtual const reco::Candidate* daughterUnsafe(size_t i) const {
      const reco::Candidate* output = NULL;
      return output;
    }

    virtual const reco::CandidatePtr daughterPtrUnsafe(size_t i) const {
      reco::CandidatePtr output;
      return output;
    }

    size_t numberOfDaughters() const { return 0; }

    virtual reco::CandidatePtr daughterUserCandUnsafe(size_t i,
        const std::string& tag) const {
      reco::CandidatePtr output;
      return output;
    }

    virtual const reco::CandidatePtrVector& daughterOverlaps(
        size_t i, const std::string& label) const {
      throw cms::Exception("NullOverlaps") <<
        "FinalState::daughterOverlaps(" << i << "," << label
        << ") is null!" << std::endl;
    }

};

template<class T>
int charge(T first) {
  return first->charge();
}

template<class T, class... Rest>
int charge(T first, Rest... rest) {
  return first->charge() + charge(rest...);
}

template<class T>
reco::Candidate::LorentzVector vector(T first) {
  return first->p4();
}

template<class T, class... Rest>
reco::Candidate::LorentzVector vector(T first, Rest... rest) {
  return first->p4() + vector(rest...);
}

template<class T, class... Rest>
class FinalStateT : public FinalState {
  public:
    typedef T daughter_type;

    FinalStateT():FinalState(){}

    FinalStateT(const edm::Ptr<T>& p1, const edm::Ptr<Rest>&... pRest,
        const edm::Ptr<FinalStateEvent>& evt)
      :FinalState(charge(p1,pRest), vector(p1,pRest), evt) {
        p1_ = p1;
        rest_ = new FinalStateT(pRest,evt);
      }

    virtual FinalStateT<T, Rest>* clone() const {
      return new FinalStateT<T, Rest>(*this);
    }

    virtual const reco::Candidate* daughterUnsafe(size_t i) const {
      const reco::Candidate* output = NULL;
      if (i == 0)
        output = p1_.get();
      else
        output = rest_.daughterUnsafe(i-1);
      return output;
    }

    virtual const reco::CandidatePtr daughterPtrUnsafe(size_t i) const {
      reco::CandidatePtr output;
      if (i == 0)
        output = p1_;
      else
        output = rest_.daughterPtrUnsafe(i-1);
      return output;
    }

    size_t numberOfDaughters() const { return 1+sizeof...(Rest); }

    virtual reco::CandidatePtr daughterUserCandUnsafe(size_t i,
        const std::string& tag) const {
      reco::CandidatePtr output;
      if (i == 0)
        output = p1_->userCand(tag);
      else
        output = rest_.daughterUserCandUnsafe(i-1,tag);
      return output;
    }

    virtual const reco::CandidatePtrVector& daughterOverlaps(
        size_t i, const std::string& label) const {
      if (i == 0)
        return p1_->overlaps(label);
      else
        return rest_.daughterOverlaps(i-1,label);
      throw cms::Exception("NullOverlaps") <<
        "FinalState::daughterOverlaps(" << i << "," << label
        << ") is null!" << std::endl;
    }

  private:
    edm::Ptr<T> p1_;
    FinalStateT<Rest> rest_;
};


#endif /* end of include guard: FinalStateAnalysis_DataFormats_FinalStateT_h */
