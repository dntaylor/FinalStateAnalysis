import FWCore.ParameterSet.Config as cms

electronsWZID = cms.EDProducer(
    'PATElectronWZEmbedder',
    src = cms.InputTag("fixme"),
    srcVertices = cms.InputTag("selectedPrimaryVertex"),
    sigmaEtaEta = cms.vdouble(0.01,0.03,0.01,0.03),
    deltaEta = cms.vdouble(0.007,0.009,0.007,0.009),
    deltaPhi = cms.vdouble(0.15,0.1,0.15,0.1),
    hoE = cms.vdouble(0.12,0.1,0.12,0.1),
    id = cms.string("WZID"),
)
