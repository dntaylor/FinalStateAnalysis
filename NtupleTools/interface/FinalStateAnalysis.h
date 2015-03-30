/*
 * Wrapper about FinalStateSelection which handles the getting objects from
 * the event.
 */

#ifndef PATFINALSTATEANALYSIS_FRM3UCVB
#define PATFINALSTATEANALYSIS_FRM3UCVB

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Provenance/interface/RunID.h"

#include "FinalStateAnalysis/DataFormats/interface/FinalStateEventFwd.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"

class TH1;
class TTree;
class FinalStateSelection;
class TFileDirectory;
namespace edm {
  class LuminosityBlockBase;
}

class FinalStateAnalysis : public edm::BasicAnalyzer {
  public:
    FinalStateAnalysis(const edm::ParameterSet& pset, TFileDirectory& fs);
    virtual ~FinalStateAnalysis();
    void beginJob() {}
    void endJob();
    // Alias for filter with no return value
    void analyze(const edm::EventBase& evt);
    bool filter(const edm::EventBase& evt);
    // Do nothing at beginning
    void beginLuminosityBlock(const edm::LuminosityBlockBase& ls){};
    void endLuminosityBlock(const edm::LuminosityBlockBase& ls);

  private:
    edm::InputTag src_;
    std::string name_;
    TFileDirectory& fs_;
    edm::ParameterSet analysisCfg_;
    boost::shared_ptr<FinalStateSelection> analysis_;

    // Tools for applying event weights
    typedef StringObjectFunction<FinalStateEvent> EventFunction;
    edm::InputTag evtSrc_;
    std::vector<EventFunction> evtWeights_;

    // Tool for examining individual runs
    bool splitRuns_;
    std::auto_ptr<TFileDirectory> runDir_;
    typedef std::map<edm::RunNumber_t, boost::shared_ptr<FinalStateSelection> > RunMap;
    RunMap runAnalysis_;

    // For counting events
    TH1* eventCounter_;
    TH1* eventCounterWeighted_;
    TH1* eventWeights_;
    // For keeping track of the skimming
    edm::InputTag skimCounter_;
    TH1* skimEventCounter_;
    // For counting the luminosity
    edm::InputTag lumiProducer_;
    TH1* integratedLumi_;

    // Keep track of the processed events in each lumi in a tree
    TTree* metaTree_;
    Int_t treeRunBranch_;
    Int_t treeLumiBranch_;
    Int_t treeEventsProcessedBranch_;
    Float_t treeIntLumi_; // The estimated integrated luminosity

    bool filter_;
};

#endif /* end of include guard: PATFINALSTATEANALYSIS_FRM3UCVB */
