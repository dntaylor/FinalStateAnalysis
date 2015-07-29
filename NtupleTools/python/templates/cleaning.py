'''

Ntuple branch template sets for applying cleaning and extra object vetoes

Each string is transformed into an expression on a FinalStateEvent object.

{object} should be replaced by an expression which evaluates to a pat::Muon
i.e. daughter(1) or somesuch.

Author: Evan K. Friis

'''

from FinalStateAnalysis.Utilities.cfgtools import PSet

# Vetos on extra stuff in the event
vetos = PSet(
    #MUON VETOS
    muVetoPt5 = 'vetoMuons(0.4, "pt > 5 & abs(eta) < 2.4").size()',
    muGlbIsoVetoPt10 = 'vetoMuons(0.4, "isGlobalMuon & isTrackerMuon & pt > 10 & abs(eta) < 2.4 & (userIso(0) + max(photonIso + neutralHadronIso - 0.5*puChargedHadronIso, 0))/pt < 0.4").size()',
    muVetoPt5IsoIdVtx = 'vetoMuons(0.4, "pt > 5 & abs(eta) < 2.4 & userInt(\'tightID\') > 0.5 & ((userIso(0) + max(photonIso()+neutralHadronIso()-0.5*puChargedHadronIso,0.0))/pt()) < 0.15 & userFloat(\'dz\') < 0.2").size()',
    muVetoPt15IsoIdVtx = 'vetoMuons(0.4, "pt > 15 & abs(eta) < 2.4 & userInt(\'tightID\') > 0.5 & ((userIso(0) + max(photonIso()+neutralHadronIso()-0.5*puChargedHadronIso,0.0))/pt()) < 0.15 & userFloat(\'dz\') < 0.2").size()',
    muonVeto4l = 'vetoMuons(0.4, "(isGlobalMuon | isTrackerMuon) & pt > 10 & abs(eta) < 2.4 & (userIso(0) + max(photonIso + neutralHadronIso - 0.5*puChargedHadronIso, 0))/pt < 0.4 & userFloat(\'dz\') < 1.0 & userFloat(\'ip3DS\') < 4.0 & userFloat(\'ipDXY\') < 0.5").size()',
    muonVetoWZLooseWideCone = 'vetoMuons(0.4, "userInt(\'tightID\') & pt > 10 & abs(eta) < 2.4 & (userIso(0) + max(photonIso + neutralHadronIso - 0.5*puChargedHadronIso, 0))/pt < 0.2").size()',
    muonVetoWZTightWideCone = 'vetoMuons(0.4, "userInt(\'tightID\') & pt > 10 & abs(eta) < 2.4 & (userIso(0) + max(photonIso + neutralHadronIso - 0.5*puChargedHadronIso, 0))/pt < 0.12").size()',
    muonVetoWZLoose = 'vetoMuons(0.1, "userInt(\'tightID\') & pt > 10 & abs(eta) < 2.4 & (userIso(0) + max(photonIso + neutralHadronIso - 0.5*puChargedHadronIso, 0))/pt < 0.2").size()',
    muonVetoWZTight = 'vetoMuons(0.1, "userInt(\'tightID\') & pt > 10 & abs(eta) < 2.4 & (userIso(0) + max(photonIso + neutralHadronIso - 0.5*puChargedHadronIso, 0))/pt < 0.12").size()',
    
    #TAU VETOS
    #OLD DMs
    tauVetoPt20Loose3HitsVtx = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFinding\') & tauID(\'byLooseCombinedIsolationDeltaBetaCorr3Hits\') & userFloat(\'dz\') < 0.2").size()',
    tauVetoPt20VLooseHPSVtx  = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFinding\') & tauID(\'byVLooseCombinedIsolationDeltaBetaCorr\') & userFloat(\'dz\') < 0.2").size()',
    tauVetoPt20TightMVALTVtx = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFinding\') & tauID(\'byTightIsolationMVA3oldDMwLT\') & userFloat(\'dz\') < 0.2").size()',
    tauVetoPt20TightMVAVtx   = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFinding\') & tauID(\'byTightIsolationMVA3oldDMwoLT\') & userFloat(\'dz\') < 0.2").size()',

    #NEW DMs
    tauVetoPt20Loose3HitsNewDMVtx = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFindingNewDMs\') & tauID(\'byLooseCombinedIsolationDeltaBetaCorr3Hits\') & userFloat(\'dz\') < 0.2").size()',
    tauVetoPt20VLooseHPSNewDMVtx  = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFindingNewDMs\') & tauID(\'byVLooseCombinedIsolationDeltaBetaCorr\') & userFloat(\'dz\') < 0.2").size()',
    tauVetoPt20TightMVALTNewDMVtx = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFinding\') & tauID(\'byTightIsolationMVA3newDMwLT\') & userFloat(\'dz\') < 0.2").size()',
    tauVetoPt20TightMVANewDMVtx   = 'vetoTaus(0.4, "pt > 20 & abs(eta) < 2.5 & tauID(\'decayModeFinding\') & tauID(\'byTightIsolationMVA3newDMwoLT\') & userFloat(\'dz\') < 0.2").size()',
    
    #ELECTRON VETOS
    eVetoMVAIsoVtx = 'vetoElectrons(0.4, "pt > 10 & abs(eta) < 2.5 & userInt(\'mvaidwp\') > 0.5 & ((userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt) < 0.3 & userFloat(\'dz\') < 0.2").size()',
    eVetoMVAIso = 'vetoElectrons(0.4, "pt > 10 & abs(eta) < 2.5 & userInt(\'mvaidwp\') > 0.5 & (userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt < 0.3").size()',
    eVetoCicTightIso = 'vetoElectrons(0.4, "pt > 10 & abs(eta) < 2.5 &  test_bit(electronID(\'cicTight\'), 0) > 0.5 & (userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt < 0.3").size()',
    elecVeto4l = 'vetoElectrons(0.4, "pt > 10 & abs(eta) < 2.5 & userInt(\'mvaidwp\') > 0.5 & ((userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt) < 0.4 & userFloat(\'dz\') < 1.0 & userFloat(\'ip3DS\') < 4.0 & userFloat(\'ipDXY\') < 0.5 & gsfTrack.trackerExpectedHitsInner.numberOfHits() <= 1").size()',
    elecVetoWZLooseWideCone = 'vetoElectrons(0.4, "userFloat(\'WZID\') & !userFloat(\'hasConversion\') & pt > 10 & abs(eta) < 2.5 & (userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt < 0.2").size()',
    elecVetoWZTightWideCone = 'vetoElectrons(0.4, "userFloat(\'WZID\') & !userFloat(\'hasConversion\') & pt > 10 & abs(eta) < 2.5 & userFloat(\'ipDXY\')<0.02 & userFloat(\'dz\')<0.1 & ((abs(superCluster().eta())<0.8 & electronID(\'mvaTrigV0\')>0.94) | (0.8<abs(superCluster().eta())<1.479 & electronID(\'mvaTrigV0\')>0.85) | (abs(superCluster().eta())>1.479 & electronID(\'mvaTrigV0\')>0.92)) & (userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt < 0.15").size()',
    elecVetoWZLoose = 'vetoElectrons(0.1, "userFloat(\'WZID\') & !userFloat(\'hasConversion\') & pt > 10 & abs(eta) < 2.5 & (userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt < 0.2").size()',
    elecVetoWZTight = 'vetoElectrons(0.1, "userFloat(\'WZID\') & !userFloat(\'hasConversion\') & pt > 10 & abs(eta) < 2.5 & userFloat(\'ipDXY\')<0.02 & userFloat(\'dz\')<0.1 & ((abs(superCluster().eta())<0.8 & electronID(\'mvaTrigV0\')>0.94) | (0.8<abs(superCluster().eta())<1.479 & electronID(\'mvaTrigV0\')>0.85) | (abs(superCluster().eta())>1.479 & electronID(\'mvaTrigV0\')>0.92)) & (userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt < 0.15").size()',
    
    #B-JET Vetos
    bjetVeto = 'vetoJets(0.4, "pt > 20 & abs(eta) < 2.4  & userInt(\'fullIdLoose\') & bDiscriminator(\'\') > 3.3").size()',
    bjetCSVVeto = 'vetoJets(0.4, "pt > 20 & abs(eta) < 2.4 & userInt(\'fullIdLoose\') & bDiscriminator(\'combinedSecondaryVertexBJetTags\') > 0.679").size()',
    bjetCSVVetoZHLike = 'vetoJets(0.4, "pt > 20 & abs(eta) < 2.4 & userInt(\'fullIdLoose\') & bDiscriminator(\'combinedSecondaryVertexBJetTags\') > 0.898").size()',
    bjetCSVVetoZHLikeNoJetId = 'vetoJets(0.4, "pt > 20 & abs(eta) < 2.4 & bDiscriminator(\'combinedSecondaryVertexBJetTags\') > 0.898").size()',
    bjetCSVVeto30 = 'vetoJets(0.4, "pt > 30 & abs(eta) < 2.4 & userInt(\'fullIdLoose\') & bDiscriminator(\'combinedSecondaryVertexBJetTags\') > 0.679").size()',

    #JET VETOS
    jetVeto20 = 'vetoJets(0.4, "pt > 20 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto20jes_plus = 'vetoJets(0.4, "userCand(\'jes+\').pt() > 20 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto20jes_minus = 'vetoJets(0.4, "userCand(\'jes-\').pt() > 20 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto20_DR05 = 'vetoJets(0.5, "pt > 20 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto30 = 'vetoJets(0.4, "pt > 30 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto30jes_plus = 'vetoJets(0.4, "userCand(\'jes+\').pt() > 30 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto30jes_minus = 'vetoJets(0.4, "userCand(\'jes-\').pt() > 30 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto30_DR05 = 'vetoJets(0.5, "pt > 30 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto40 = 'vetoJets(0.4, "pt > 40 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    jetVeto40_DR05 = 'vetoJets(0.5, "pt > 40 & abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size()',
    #leadingJetPt = '? (vetoJets(0.4, "abs(eta) < 5.0 & userInt(\'fullIdLoose\')").size() > 0) ? vetoJets(0.4, "abs(eta) < 5.0 & userInt(\'fullIdLoose\')").at(0).pt() : -1.',
)

overlaps = PSet(
    objectElectronPt15IdIsoVtxOverlap = 'overlapElectrons({object_idx}, 0.4, "pt > 15 & abs(eta) < 2.5 & userInt(\'mvaidwp\') > 0.5 & ((userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt) < 0.3 & abs(userFloat(\'dz\')) < 0.2").size()', 
    objectElectronPt10IdIsoVtxOverlap = 'overlapElectrons({object_idx}, 0.4, "pt > 10 & abs(eta) < 2.5 & userInt(\'mvaidwp\') > 0.5 & ((userIso(0) + max(userIso(1) + neutralHadronIso - 0.5*userIso(2), 0))/pt) < 0.3 & abs(userFloat(\'dz\')) < 0.2").size()', 
    objectElectronPt15IdVtxOverlap = 'overlapElectrons({object_idx}, 0.4, "pt > 15 & abs(eta) < 2.5 & userInt(\'mvaidwp\') > 0.5 & abs(userFloat(\'dz\')) < 0.2").size()',
    objectElectronPt10IdVtxOverlap = 'overlapElectrons({object_idx}, 0.4, "pt > 10 & abs(eta) < 2.5 & userInt(\'mvaidwp\') > 0.5 & abs(userFloat(\'dz\')) < 0.2").size()',
   
    objectMuonIdIsoVtxOverlap = 'overlapMuons({object_idx}, 0.4, "pt > 10 & abs(eta) < 2.4 & userInt(\'tightID\') > 0.5 & ((userIso(0) + max(photonIso()+neutralHadronIso()-0.5*puChargedHadronIso,0.0))/pt()) < 0.15 & userFloat(\'dz\') < 0.2").size()',
    objectMuonIdVtxOverlap = 'overlapMuons({object_idx}, 0.4, "pt > 10 & abs(eta) < 2.4 & userInt(\'tightID\') > 0.5 & userFloat(\'dz\') < 0.2").size()',
    objectMuonIdIsoStdVtxOverlap = 'overlapMuons({object_idx}, 0.4, "pt > 10 & abs(eta) < 2.4 & userInt(\'tightID\') > 0.5 & ((chargedHadronIso() + max(photonIso()+neutralHadronIso()-0.5*puChargedHadronIso,0.0))/pt()) < 0.15 & userFloat(\'dz\') < 0.2").size()',
    objectGlobalMuonVtxOverlap = 'overlapMuons({object_idx}, 0.4, "pt > 10 & abs(eta) < 2.4 & isGlobalMuon & userFloat(\'dz\') < 0.2").size()',

    objectMuOverlap = 'overlapMuons({object_idx}, 0.4, "pt > 5").size()',
    objectElecOverlap = 'overlapElectrons({object_idx}, 0.4, "pt > 10").size()',
    objectCiCTightElecOverlap = 'overlapElectrons({object_idx}, 0.4, "pt > 10 & test_bit(electronID(\'cicTight\'), 0)").size()',
)
