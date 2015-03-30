import FWCore.ParameterSet.Config as cms

finalStateMassResolutionEmbedder = cms.EDProducer(
    "FinalStateMassResolutionEmbedder",
    debug = cms.bool(False),
    src = cms.InputTag("fixme")
    )
