/*
 * Since FinalStateEvent does not inherit from reco::Candidate,
 * we can't use the CandViewHistoAnalyzer.  This is its replacment.
 *
 * Author: Evan K. Friis, UW Madison
 *
 */

#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/HistoAnalyzer.h"
#include "FinalStateAnalysis/DataFormats/interface/FinalStateEvent.h"
#include "FinalStateAnalysis/DataFormats/interface/FinalStateEventFwd.h"

typedef HistoAnalyzer<FinalStateEventCollection> FinalStateEventHistoAnalyzer;

DEFINE_FWK_MODULE( FinalStateEventHistoAnalyzer );
