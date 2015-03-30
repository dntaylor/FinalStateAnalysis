import FWCore.ParameterSet.Config as cms

'''
Build the FinalStateLS wrapper data class which holds the lumi summary.
'''

finalStateLS = cms.EDProducer(
    "FinalStateLSProducer",
    lumiSrc = cms.InputTag("lumiProducer"),
    trigSrc = cms.InputTag("patTriggerEvent"),
    xSec = cms.double(-1),
)
