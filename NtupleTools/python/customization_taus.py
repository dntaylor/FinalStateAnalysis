# Embed IDs for taus
import FWCore.ParameterSet.Config as cms

def preTaus(process, use25ns, tSrc, vSrc,**kwargs):
    postfix = kwargs.pop('postfix','')
    modName = 'genembeddedTaus{0}'.format(postfix)
    mod=cms.EDProducer("PATTauGenInfoEmbedder",
          src=cms.InputTag(tSrc)
    )
    setattr(process,modName,mod)
    tSrc = modName
    modPath = 'embeddedTaus{0}'.format(postfix)
    setattr(process,modPath,cms.Path(getattr(process,modName)))
    
    process.schedule.append(getattr(process,modPath))

    return tSrc

def postTaus(process, use25ns, tSrc, jSrc,**kwargs):
    postfix = kwargs.pop('postfix','')
    modName = 'miniAODTauJetInfoEmbedding{0}'.format(postfix)
    mod = cms.EDProducer(
        "MiniAODTauJetInfoEmbedder",
        src = cms.InputTag(tSrc),
        embedBtags = cms.bool(False),
        suffix = cms.string(''),
        jetSrc = cms.InputTag(jSrc),
        maxDeltaR = cms.double(0.1),
    )
    setattr(process,modName,mod)
    tSrc = modName
    modPath = 'TauJetInfoEmbedding{0}'.format(postfix)
    setattr(process,modPath,cms.Path(getattr(process,modName)))
    process.schedule.append(getattr(process,modPath))

    return tSrc

