#!/usr/bin/env cmsRun
'''

Ntuple Generation
=================

Generates the ntuples for a given list of final state generation.

Usage:

./make_ntuples_cfg.py channels="mt,em,mm,eemm" [options]

There are some additional pre-defined groups of channels which are expanded
for your convenience::

zh = eeem, eeet, eemt, eett, emmm, emmt, mmmt, mmtt,
zz = eeee, eemm, mmmm,
zgg = eegg, mmgg
llt = emt, mmt, eet, mmm, emm
zg = mmg,eeg
zgxtra = mgg, emg, egg,

The available options (which are set to zero or one) are::

skipEvents=0   - events to skip (for debugging)
maxEvents=-1   - events to run on
rerunMCMatch=0 - rerun MC matching
eventView=0    - make a row in the ntuple correspond to an event
                 instead of a final state in an event.
passThru=0     - turn off any preselection/skim
dump=0         - if one, dump process python to stdout
verbose=0      - print out timing information
noPhotons=0    - don't build things which depend on photons.
rerunMVAMET=0  - rerun the MVAMET algorithm
svFit=1        - run the SVfit on appropriate pairs
rerunQGJetID=0 - rerun the quark-gluon JetID
rerunJets=0    - rerun with new jet energy corrections
runMetFilter=0 - apply met filters
use25ns=1      - run on 25 ns miniAOD (0 -> 50ns)
runDQM=0       - run over single object final states to test all 
                 object properties (wont check diobject properties)
hzz=0          - Include FSR contribution a la HZZ4l group, 
                 include all ZZ candidates (including alternative lepton pairings).
isSync=0       - (with eCalib=1 and isMC=1) Apply electron energy 
                 resolution corrections as a 1sigma shift instead of smearing 
                 for synchronization purposes.
eCalib=0       - Apply electron energy scale and resolution corrections.
nExtraJets=0   - Include basic info about this many jets (ordered by pt). 
                 Ignored if final state involves jets.
paramFile=''   - custom parameter file for ntuple production
keepPat=0      - Instead of making flat ntuples, write high level 
                 physics objects including the PATFinalState objects
skipMET=0      - Don't do MET corrections and systematics (good way to reduce 
                 memory use if you don't use them)

'''

import FWCore.ParameterSet.Config as cms
import os
import copy
from FinalStateAnalysis.NtupleTools.hzg_sync_mod import set_passthru
from FinalStateAnalysis.NtupleTools.ntuple_builder import \
    make_ntuple, add_ntuple, _producer_translation
from FinalStateAnalysis.Utilities.version import cmssw_major_version, \
    cmssw_minor_version
import PhysicsTools.PatAlgos.tools.helpers as helpers

process = cms.Process("Ntuples")

# if you want to debug in the future, uncomment this
#process.ProfilerService = cms.Service (
#      "ProfilerService",
#       firstEvent = cms.untracked.int32(2),
#       lastEvent = cms.untracked.int32(500),
#       paths = cms.untracked.vstring('schedule') 
#)
#
#process.SimpleMemoryCheck = cms.Service(
#    "SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)


process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
)

import FinalStateAnalysis.Utilities.TauVarParsing as TauVarParsing
options = TauVarParsing.TauVarParsing(
    skipEvents=0,  # Start at an event offset (for debugging)
    reportEvery=100,
    channels='mt,em,et,mm,ee',
    rerunMCMatch=False,
    eventView=0,  # Switch between final state view (0) and event view (1)
    passThru=0,  # Turn off preselections
    dump=0,  # If one, dump process python to stdout
    verbose=0,  # If one print out the TimeReport
    noPhotons=0,  # If one, don't assume that photons are in the PAT tuples.
    svFit=0,  # If one, SVfit appropriate lepton pairs.
    rochCor="",
    eleCor="",
    rerunQGJetID=0,  # If one reruns the quark-gluon JetID
    runMVAMET=0,  # If one, (re)build the MVA MET
    runMETNoHF=0,  # If one, use get metnohf (needs to be recalculated in miniaodv1)
    usePUPPI=0,
    rerunJets=0,
    dblhMode=False, # For double-charged Higgs analysis
    runTauSpinner=0,
    GlobalTag="",
    runMetFilter=0,
    use25ns=1,
    runDQM=0,
    hzz=0,
    paramFile='',
    skipGhost=0,
    runWZ=0,
    runMetUncertainties=0,
    metShift='',
    runFSRFilter=0, # 1 = filter for ZG, -1 inverts filter for DY
    eventsToSkip='',
)

options.register(
    'skimCuts',
    '',
    TauVarParsing.TauVarParsing.multiplicity.list,
    TauVarParsing.TauVarParsing.varType.string,
    'additional cuts to impose on the NTuple'
)
options.register(
    'nExtraJets',
    0,
    TauVarParsing.TauVarParsing.multiplicity.singleton,
    TauVarParsing.TauVarParsing.varType.int,
    'Number of pt-ordered jets to keep some info about. Ignored if final state involves jets.',
)
options.register(
    'keepPat',
    0,
    TauVarParsing.TauVarParsing.multiplicity.singleton,
    TauVarParsing.TauVarParsing.varType.int,
    'If 1, write final PAT objects rather than a flat ntuple. '
    'If 2, also keep packedGenParticles. If 3 or more, also keep '
    'packedPFCands (increases size significantly).',
)
options.register(
    'eCalib',
    0,
    TauVarParsing.TauVarParsing.multiplicity.singleton,
    TauVarParsing.TauVarParsing.varType.int,
    'Apply electron energy scale and resolution corrections. '
    'For data, this is a correction; for MC, it is a smearing'
)
options.register(
    'isSync',
    0,
    TauVarParsing.TauVarParsing.multiplicity.singleton,
    TauVarParsing.TauVarParsing.varType.int,
    'Apply electron energy correction as a 1-sigma shift instead of a '
    'smearing. Only used if eCalib=0 and isMC=1.'
)
options.register(
    'skipMET',
    0,
    TauVarParsing.TauVarParsing.multiplicity.singleton,
    TauVarParsing.TauVarParsing.varType.int,
    "Skip MET corrections and systematics (good way to reduce memory "
    "use if you don't need them)."
)

options.outputFile = "ntuplize.root"
options.parseArguments()

#########################
### Customize the job ###
#########################

# list of filters to apply
filters = []

# SV Fit requires MVA MET
options.runMVAMET = (options.runMVAMET or options.svFit)

eventsToSkip = cms.untracked.VEventRange()
if options.eventsToSkip:
    # if it is a file
    if os.path.isfile(options.eventsToSkip):
        with open(options.eventsToSkip,'r') as f:
            for e in f:
                eventsToSkip.append(e.rstrip())
    else:
        print 'Warning: event file {0} does not exist.'.format(options.eventsToSkip)


process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
    skipEvents=cms.untracked.uint32(options.skipEvents),
    eventsToSkip=eventsToSkip,
)

from FinalStateAnalysis.NtupleTools.parameters.default import parameters
if options.paramFile:
    # add custom parameters
    if os.path.isfile(options.paramFile):
        print 'Using custom parameter file %s' % os.path.abspath(options.paramFile)
        import imp
        custParamModule = imp.load_source('custParamModule',options.paramFile)
        from custParamModule import parameters as custParams
        parameters.update(custParams)
    else:
        print 'Failed to load custom parameters, using default.'
    pass
if options.runWZ:
    from FinalStateAnalysis.NtupleTools.parameters.wz import parameters as wzParams
    parameters.update(wzParams)
    options.eCalib=1
    

if options.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(
        options.eventsToProcess)

# If desired, apply a luminosity mask
if options.lumiMask:
    print "Applying LumiMask from", options.lumiMask
    process.source.lumisToProcess = options.buildPoolSourceLumiMask()

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string(options.outputFile)
)

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents))

process.schedule = cms.Schedule()

### To use IgProf's neat memory profiling tools, uncomment the following 
### lines then run this cfg with igprof like so:
###      $ igprof -d -mp -z -o igprof.mp.gz cmsRun ... 
### this will create a memory profile every 250 events so you can track use
### Turn the profile into text with
###      $ igprof-analyse -d -v -g -r MEM_LIVE igprof.yourOutputFile.gz > igreport_live.res
### To do a performance profile instead of a memory profile, change -mp to -pp
### in the first command and remove  -r MEM_LIVE from the second
### For interpretation of the output, see http://igprof.org/text-output-format.html

#from IgTools.IgProf.IgProfTrigger import igprof
#process.load("IgTools.IgProf.IgProfTrigger")
#process.igprofPath = cms.Path(process.igprof)
#process.igprof.reportEventInterval     = cms.untracked.int32(250)
#process.igprof.reportToFileAtBeginJob  = cms.untracked.string("|gzip -c>igprof.begin-job.gz")
#process.igprof.reportToFileAtEvent = cms.untracked.string("|gzip -c>igprof.%I.%E.%L.%R.event.gz")
#process.schedule.append(process.igprofPath)

#load magfield and geometry (for mass resolution)
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')

# Need the global tag for geometry etc.
envvar = 'mcgt' if options.isMC else 'datagt'
GT = {'mcgt': '74X_mcRun2_asymptotic_v4', 'datagt': '74X_dataRun2_reMiniAOD_v1'}

process.GlobalTag.globaltag = cms.string(GT[envvar])

print 'Using globalTag: %s' % process.GlobalTag.globaltag

# Count events at the beginning of the tuplization
process.load("FinalStateAnalysis.RecoTools.eventCount_cfi")
process.load("FinalStateAnalysis.PatTools.finalStates.patFinalStateLSProducer_cfi")
process.generateMetaInfo = cms.Path(process.eventCount *
                                    process.summedWeight *
                                    process.finalStateLS
                                    )
process.schedule.append(process.generateMetaInfo)

######################################################################
### Uncomment if a hard process tau skim is desired
# process.load("FinalStateAnalysis.NtupleTools.genTauFilter_cfi")
# filters.append(process.filterForGenTaus)
######################################################################

# Drop the input ones, just to make sure we aren't screwing anything up
process.buildFSASeq = cms.Sequence()
from FinalStateAnalysis.PatTools.patFinalStateProducers \
    import produce_final_states
# Which collections are used to build the final states
fs_daughter_inputs = {
    'electrons': 'slimmedElectrons',
    'muons': 'slimmedMuons',
    'taus': 'slimmedTaus',
    'photons': 'slimmedPhotons',
    'jets': 'slimmedJets',
    'pfmet': 'slimmedMETs',         # slimmedMETs, slimmedMETsNoHF (miniaodv2), slimmmedMETsPuppi (not correct in miniaodv1)
    'mvamet': 'fixme',              # produced later
    'vertices': 'offlineSlimmedPrimaryVertices',
}

# add additional final states to ntuples with different parameters... work in progress... difficult with VID
additional_fs = {}



###############################
### Customize object inputs ###
###############################

# use metNoHF
if options.runMETNoHF:
    fs_daughter_inputs['pfmet'] = 'slimmedMETsNoHF'

# use puppi instead
if options.usePUPPI:
    fs_daughter_inputs['pfmet'] = 'slimmedMETsPuppi'
    fs_daughter_inputs['jets'] = 'slimmedJetsPuppi'


##############
### mvamet ###
##############

if options.runMVAMET:
    process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    process.ak4PFJets.src = cms.InputTag("packedPFCandidates")
    process.ak4PFJets.doAreaFastjet = cms.bool(True)
    
    from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3
    
    process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
    process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
    process.pfMVAMEt.srcVertices = cms.InputTag(fs_daughter_inputs['vertices'])
    process.pfMVAMEt.inputFileNames.U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_25NS_July2015.root')
    process.pfMVAMEt.inputFileNames.DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_25NS_July2015.root')
    process.pfMVAMEt.inputFileNames.CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_25NS_July2015.root')
    process.pfMVAMEt.inputFileNames.CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_25NS_July2015.root')
    if not options.use25ns:
        process.pfMVAMEt.inputFileNames.U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_50NS_July2015.root')
        process.pfMVAMEt.inputFileNames.DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_50NS_July2015.root')
        process.pfMVAMEt.inputFileNames.CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_50NS_July2015.root')
        process.pfMVAMEt.inputFileNames.CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_50NS_July2015.root')
    
    process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
    process.puJetIdForPFMVAMEt.vertexes = cms.InputTag(fs_daughter_inputs['vertices'])
    process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")
    
    from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
    
    process.miniAODMVAMEt = patMETs.clone(
        metSource=cms.InputTag("pfMVAMEt"),
        addMuonCorrections = cms.bool(False),
        addGenMET = cms.bool(False)
    )
    fs_daughter_inputs['mvamet'] = 'miniAODMVAMEt'
    
    process.mvaMetSequence = cms.Path(
        process.ak4PFJets *
        process.pfMVAMEtSequence *
        process.miniAODMVAMEt
    )



################################################
### add filters (that wont make it into fsa) ###
################################################

# add met filters
if options.runMetFilter:
    # HBHE Loose
    # flag in miniaod wrong, must rerun
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
    process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
    
    process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
       inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
       reverseDecision = cms.bool(False)
    )
    
    process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
       inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
       reverseDecision = cms.bool(False)
    )
    filters += [process.HBHENoiseFilterResultProducer, process.ApplyBaselineHBHENoiseFilter]

    # CSC Tight Halo and EE Bad SC 4
    eventListFile = 'FinalStateAnalysis/NtupleTools/data/csc2015_Dec01.txt'
    process.MiniAODCSCTightHaloFilterProducer = cms.EDProducer('MiniAODEventListProducer',
        label = cms.string('CSCTightHaloFilterResult'),
        eventList = cms.FileInPath(eventListFile),
    )
    process.ApplyCSCTightHaloFilter = cms.EDFilter('BooleanFlagFilter',
        inputLabel = cms.InputTag('MiniAODCSCTightHaloFilterProducer','CSCTightHaloFilterResult'),
        reverseDecision = cms.bool(False),
    )
    filters += [process.MiniAODCSCTightHaloFilterProducer, process.ApplyCSCTightHaloFilter]

    eventListFile = 'FinalStateAnalysis/NtupleTools/data/ecalscn1043093_Dec01.txt'
    process.MiniAODBadSCEE4FilterProducer = cms.EDProducer('MiniAODEventListProducer',
        label = cms.string('BadSCEE4FilterResult'),
        eventList = cms.FileInPath(eventListFile),
    )
    process.ApplyBadSCEE4Filter = cms.EDFilter('BooleanFlagFilter',
        inputLabel = cms.InputTag('MiniAODBadSCEE4FilterProducer','BadSCEE4FilterResult'),
        reverseDecision = cms.bool(False),
    )
    filters += [process.MiniAODBadSCEE4FilterProducer, process.ApplyBadSCEE4Filter]

    # good vertices and ee bad sc filter
    # flag in miniaod, so just filter on that
    listOfFlags = ['Flag_goodVertices', 'Flag_eeBadScFilter']
    listOfLabels = ['GoodVerticesFilterResult', 'EEBadSCFilterResult']
    process.MiniAODMETFilterProducer = cms.EDProducer('MiniAODTriggerProducer',
        triggers = cms.vstring(*listOfFlags),
        labels = cms.vstring(*listOfLabels),
        bits = cms.InputTag("TriggerResults"),
        #prescales = cms.InputTag("patTrigger"),
        #objects = cms.InputTag("selectedPatTrigger"),
    )
    process.ApplyGoodVerticesFilter = cms.EDFilter('BooleanFlagFilter',
       inputLabel = cms.InputTag('MiniAODMETFilterProducer','GoodVerticesFilterResult'),
       reverseDecision = cms.bool(True)
    )
    process.ApplyEEBadSCFilter = cms.EDFilter('BooleanFlagFilter',
       inputLabel = cms.InputTag('MiniAODMETFilterProducer','EEBadSCFilterResult'),
       reverseDecision = cms.bool(True)
    )
    filters += [process.MiniAODMETFilterProducer, process.ApplyGoodVerticesFilter, process.ApplyEEBadSCFilter]

if abs(options.runFSRFilter)>0:
    process.FSRFilter = cms.EDFilter("MiniAODGenLeptonFSRFilter",
        src = cms.InputTag("prunedGenParticles"),
        drCut = cms.double(0.05),
        #drCut = cms.double(0.4),
    )
    #if options.runFSRFilter<0:
    if options.runFSRFilter>0:
        process.FSRFilter.reverseDecision=cms.bool(True)
    filters += [process.FSRFilter]


#######################################
### MET Uncertainty and Corrections ###
#######################################

if not bool(options.skipMET):
    postfix = 'NewMet'
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    isData = not options.isMC
    runMetCorAndUncFromMiniAOD(process,
                               jetColl=fs_daughter_inputs['jets'],
                               jetCollUnskimmed='slimmedJets',
                               photonColl=fs_daughter_inputs['photons'],
                               electronColl=fs_daughter_inputs['electrons'],
                               muonColl=fs_daughter_inputs['muons'],
                               tauColl=fs_daughter_inputs['taus'],
                               isData=isData,
                               jecUncFile='FinalStateAnalysis/NtupleTools/data/Summer15_25nsV6_{0}_UncertaintySources_AK4PFchs.txt'.format('MC' if options.isMC else 'DATA'),
                               repro74X=True,
                               postfix=postfix,
                               )
    
    collMap = {
        #'jres' : {'Jets'     : 'shiftedPatJetRes{sign}{postfix}'},
        'jres' : {},
        'jes'  : {'Jets'     : 'shiftedPatJetEn{sign}{postfix}'},
        'mes'  : {'Muons'    : 'shiftedPatMuonEn{sign}{postfix}'},
        'ees'  : {'Electrons': 'shiftedPatElectronEn{sign}{postfix}'},
        'tes'  : {'Taus'     : 'shiftedPatTauEn{sign}{postfix}'},
        'ues'  : {},
        'pes'  : {},
    }
    signMap = {
      '+' : 'Up',
      '-' : 'Down',
    }
    metMap = {
      'jres' : 'patPFMetT1JetRes{sign}{postfix}',
      'jes'  : 'patPFMetT1JetEn{sign}{postfix}',
      'mes'  : 'patPFMetT1MuonEn{sign}{postfix}',
      'ees'  : 'patPFMetT1ElectronEn{sign}{postfix}',
      'tes'  : 'patPFMetT1TauEn{sign}{postfix}',
      'ues'  : 'patPFMetT1UnclusteredEn{sign}{postfix}',
      'pes'  : '',
    }
    allowedShifts = ['jres','jes','mes','ees','tes','ues']
    allowedSigns = ['+','-']
    
    # this should be it, but fails
    #process.applyCorrections = cms.Path(getattr(process,'fullPatMetSequence{0}'.format(postfix)))
    # fix things
    getattr(process,'patPFMetT1T2Corr{0}'.format(postfix)).src = cms.InputTag('patJets')
    getattr(process,'patPFMetT2Corr{0}'.format(postfix)).src = cms.InputTag('patJets')
    getattr(process,'patPFMetTxyCorr{0}'.format(postfix)).vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
    process.applyCorrections = cms.Path()
    if options.isMC: process.applyCorrections += process.genMetExtractor
    process.applyCorrections += getattr(process,'patPFMet{0}'.format(postfix))
    process.applyCorrections += process.patJetCorrFactorsReapplyJEC
    process.applyCorrections += process.patJets
    process.applyCorrections += getattr(process,'patPFMetT1Txy{0}'.format(postfix))
    process.applyCorrections += getattr(process,'patPFMetT1{0}'.format(postfix))
    process.applyCorrections += getattr(process,'patPFMetTxy{0}'.format(postfix))
    process.schedule.append(process.applyCorrections)
    fs_daughter_inputs['jets'] = 'patJets'
    
    # embed references to shifts
    process.embedShifts = cms.Path()
    for shift in allowedShifts:
        for sign in allowedSigns:
            # embed shifted objects
            for coll in collMap[shift]:
                modName = '{shift}{sign}{coll}Embedding'.format(shift=shift,sign=signMap[sign],coll=coll)
                pluginName = 'MiniAODShifted{coll}Embedder'.format(coll=coll[:-1])
                dName = coll.lower()
                srcName = fs_daughter_inputs[dName]
                shiftSrcName = collMap[shift][coll].format(sign=signMap[sign],postfix=postfix)
                label = '{shift}{sign}{coll}'.format(shift=shift,sign=signMap[sign],coll=coll)
                module = cms.EDProducer(
                    pluginName,
                    src = cms.InputTag(srcName),
                    shiftSrc = cms.InputTag(shiftSrcName),
                    label = cms.string(label),
                )
                setattr(process,modName,module)
                fs_daughter_inputs[dName] = modName
                #process.embedShifts *= getattr(process,shiftSrcName)
                process.embedShifts *= getattr(process,modName)
            # embed shifted met
            modName = '{shift}{sign}METEmbedding'.format(shift=shift,sign=signMap[sign])
            metName = metMap[shift].format(sign=signMap[sign],postfix=postfix)
            label = '{shift}{sign}MET'.format(shift=shift,sign=signMap[sign])
            module = cms.EDProducer(
                'MiniAODShiftedMETEmbedder',
                src = cms.InputTag(fs_daughter_inputs['pfmet']),
                shiftSrc = cms.InputTag(metName),
                label = cms.string(label),
            )
            setattr(process,modName,module)
            fs_daughter_inputs['pfmet'] = modName
            process.embedShifts *= getattr(process,modName)
    process.schedule.append(process.embedShifts)
    
    
    # switch input to desired one
    if options.metShift: 
        t = options.metShift[:-1]
        d = options.metShift[-1]
        if options.metShift=='all':
            # setup daughters for all
            for shift in allowedShifts:
                if shift!='mes': continue
                for sign in allowedSigns:
                    if sign!='+': continue
                    label = shift + signMap[sign]
                    additional_fs[label] = copy.deepcopy(fs_daughter_inputs)
                    additional_fs[label]['pfmet'] = metMap[shift].format(sign=signMap[sign],postfix=postfix)
                    for coll in collMap[shift]:
                        additional_fs[label][coll.lower()] = collMap[shift][coll].format(sign=signMap[sign],postfix=postfix)
        elif t not in allowedShifts or d not in allowedSigns:
            print 'Warning: {0} is not an allowed MET shift, using unshifted collections'.format(options.metShift)
        else:
            fs_daughter_inputs['pfmet'] = metMap[t].format(sign=signMap[d],postfix=postfix)
            for coll in collMap[t]:
                fs_daughter_inputs[coll.lower()] = collMap[t][coll].format(sign=signMap[d],postfix=postfix)



#process.EventAnalyzer = cms.EDAnalyzer("EventContentAnalyzer")
#process.eventAnalyzerPath = cms.Path(process.EventAnalyzer)
#process.schedule.append(process.eventAnalyzerPath)


#########################################################
### embed some things we need before object selection ###
#########################################################

# HZZ id labels
electronMVANonTrigIDLabel = "BDTIDNonTrig"
electronMVATrigIDLabel = "BDTIDTrig"

#########################################
### calibrate electrons and embed ids ###
#########################################
from FinalStateAnalysis.NtupleTools.customization_electrons import preElectrons
fs_daughter_inputs['electrons'] = preElectrons(process,
                                               options.use25ns,
                                               fs_daughter_inputs['electrons'],
                                               fs_daughter_inputs['vertices'],
                                               electronMVANonTrigIDLabel=electronMVANonTrigIDLabel,
                                               electronMVATrigIDLabel=electronMVATrigIDLabel,
                                               applyEnergyCorrections=bool(options.eCalib),
                                               isMC=bool(options.isMC),
                                               isSync=bool(options.isSync))
for fs in additional_fs:
    additional_fs[fs]['electrons'] = preElectrons(process,
                                                  options.use25ns,
                                                  additional_fs[fs]['electrons'],
                                                  additional_fs[fs]['vertices'],
                                                  electronMVANonTrigIDLabel=electronMVANonTrigIDLabel,
                                                  electronMVATrigIDLabel=electronMVATrigIDLabel,
                                                  applyEnergyCorrections=bool(options.eCalib),
                                                  isMC=bool(options.isMC),
                                                  isSync=bool(options.isSync),
                                                  postfix=fs)

######################
### embed muon IDs ###
######################
from FinalStateAnalysis.NtupleTools.customization_muons import preMuons
fs_daughter_inputs['muons'] = preMuons(process,
                                       options.use25ns,
                                       fs_daughter_inputs['muons'],
                                       fs_daughter_inputs['vertices'],
                                       skipGhost=options.skipGhost)
for fs in additional_fs:
    additional_fs[fs]['muons'] = preMuons(process,
                                          options.use25ns,
                                          additional_fs[fs]['muons'],
                                          additional_fs[fs]['vertices'],
                                          skipGhost=options.skipGhost,
                                          postfix=fs)

#####################
### embed tau IDs ###
#####################
from FinalStateAnalysis.NtupleTools.customization_taus import preTaus
fs_daughter_inputs['taus'] = preTaus(process,
                                     options.use25ns,
                                     fs_daughter_inputs['taus'],
                                     fs_daughter_inputs['vertices'])
for fs in additional_fs:
    additional_fs[fs]['taus'] = preTaus(process,
                                        options.use25ns,
                                        additional_fs[fs]['taus'],
                                        additional_fs[fs]['vertices'],
                                        postfix=fs)

########################
### jet id embedding ###
########################
from FinalStateAnalysis.NtupleTools.customization_jets import preJets
fs_daughter_inputs['jets'] = preJets(process,
                                     options.use25ns,
                                     fs_daughter_inputs['jets'],
                                     fs_daughter_inputs['vertices'],
                                     fs_daughter_inputs['muons'],
                                     fs_daughter_inputs['electrons'],
                                     eCut=parameters['preselection']['j']['e']['selection'],
                                     eDeltaR=parameters['preselection']['j']['e']['deltaR'],
                                     mCut=parameters['preselection']['j']['m']['selection'],
                                     mDeltaR=parameters['preselection']['j']['m']['deltaR'],
                                     jCut=parameters['preselection']['j']['selection'],
                                     jType="AK4PFchs")
for fs in additional_fs:
    additional_fs[fs]['jets'] = preJets(process,
                                        options.use25ns,
                                        additional_fs[fs]['jets'],
                                        additional_fs[fs]['vertices'],
                                        additional_fs[fs]['muons'],
                                        additional_fs[fs]['electrons'],
                                        eCut=parameters['preselection']['j']['e']['selection'],
                                        eDeltaR=parameters['preselection']['j']['e']['deltaR'],
                                        mCut=parameters['preselection']['j']['m']['selection'],
                                        mDeltaR=parameters['preselection']['j']['m']['deltaR'],
                                        jCut=parameters['preselection']['j']['selection'],
                                        jType="AK4PFchs",
                                        postfix=fs)


########################################
### pre selection HZZ customizations ###
########################################
idCheatLabel = "HZZ4lIDPass" # Gets loose ID. For tight ID, append "Tight".
isoCheatLabel = "HZZ4lIsoPass"
if options.hzz:
    from FinalStateAnalysis.NtupleTools.customization_hzz import hzzCustomize
    hzzCustomize(process, fs_daughter_inputs, idCheatLabel, isoCheatLabel, 
                 electronMVANonTrigIDLabel, "dretFSRCand")
    # fs_daughter_inputs entries for electrons, muons, and fsr are automatically 
    # set by hzzCustomize()



###########################
### object preselection ###
###########################
if options.passThru:
    preselections = {}
else:
    preselections = parameters.get('preselection',{})

from FinalStateAnalysis.NtupleTools.object_parameter_selector import setup_selections
process.preselectionSequence = setup_selections(
    process, 
    "Preselection",
    fs_daughter_inputs,
    preselections,
    )
process.FSAPreselection = cms.Path(process.preselectionSequence)
process.schedule.append(process.FSAPreselection)

for fs in additional_fs:
    preSeqName = 'preselectionSequence{0}'.format(fs)
    preSeq = setup_selections(
        process,
        "Preselection",
        additional_fs[fs],
        preselections,
        postfix=fs,
        )
    setattr(process,preSeqName,preSeq)
    for ob in preselections:
        additional_fs[fs][getName(ob)+'s'] = getName(ob)+"Preselection{0}".format(fs)
    setattr(process,'FSAPreselection{0}'.format(fs),cms.Path(getattr(process,preSeqName)))
    process.schedule.append(getattr(process,'FSAPreselection{0}'.format(fs)))



###########################################################################
### The following is embedding that must be done after object selection ###
###########################################################################

###################################
### post electron customization ###
###################################
from FinalStateAnalysis.NtupleTools.customization_electrons import postElectrons
fs_daughter_inputs['electrons'] = postElectrons(process,options.use25ns,fs_daughter_inputs['electrons'],fs_daughter_inputs['jets'])
for fs in additional_fs:
    additional_fs[fs]['electrons'] = postElectrons(process,options.use25ns,additional_fs[fs]['electrons'],additional_fs[fs]['jets'],postfix=fs)


###############################
### post muon customization ###
###############################
from FinalStateAnalysis.NtupleTools.customization_muons import postMuons
fs_daughter_inputs['muons'] = postMuons(process,options.use25ns,fs_daughter_inputs['muons'],fs_daughter_inputs['jets'])
for fs in additional_fs:
    additional_fs[fs]['muons'] = postMuons(process,options.use25ns,additional_fs[fs]['muons'],additional_fs[fs]['jets'],postfix=fs)

##############################
### post tau customization ###
##############################
from FinalStateAnalysis.NtupleTools.customization_taus import postTaus
fs_daughter_inputs['taus'] = postTaus(process,options.use25ns,fs_daughter_inputs['taus'],fs_daughter_inputs['jets'])
for fs in additional_fs:
    additional_fs[fs]['taus'] = postTaus(process,options.use25ns,additional_fs[fs]['taus'],additional_fs[fs]['jets'],postfix=fs)








############################
### Now do the FSA stuff ###
############################

# Make a list of collections to save (in case we're saving 
# collections instead of flat ntuples)
# Note that produce_final_states adds the PATFinalStateEvent and
# the collections it used to build the final states, but not the final
# states themselves, because that fills the file up with unwanted stuff
output_to_keep = []

# Eventually, set buildFSAEvent to False, currently working around bug
# in pat tuples.
produce_final_states(process, 
                     fs_daughter_inputs, 
                     output_to_keep, 
                     process.buildFSASeq,
                     'puTagDoesntMatter', 
                     buildFSAEvent=True,
                     noTracks=True, 
                     runMVAMET=options.runMVAMET,
                     hzz=options.hzz, 
                     rochCor=options.rochCor,
                     eleCor=options.eleCor, 
                     **parameters)
process.buildFSAPath = cms.Path(process.buildFSASeq)
# Don't crash if some products are missing (like tracks)
process.patFinalStateEventProducer.forbidMissing = cms.bool(False)
process.schedule.append(process.buildFSAPath)

for fs in additional_fs:
    setattr(process,'buildFSASeq{0}'.format(fs),cms.Sequence())
    produce_final_states(process, 
                         additional_fs[fs], 
                         output_to_keep, 
                         getattr(process,'buildFSASeq{0}'.format(fs)),
                         'puTagDoesntMatter', 
                         buildFSAEvent=True,
                         noTracks=True, 
                         runMVAMET=options.runMVAMET,
                         hzz=options.hzz, 
                         rochCor=options.rochCor,
                         eleCor=options.eleCor, 
                         postfix=fs, 
                         **parameters)
    setattr(process,'buildFSAPath{0}'.format(fs), cms.Path(getattr(process,'buildFSASeq{0}'.format(fs))))
    getattr(process,'patFinalStateEventProducer{0}'.format(fs)).forbidMissing = cms.bool(False)
    process.schedule.append(getattr(process,'buildFSAPath{0}'.format(fs)))


# Drop the old stuff. (do we still need this?)
process.source.inputCommands = cms.untracked.vstring(
    'keep *',
    'drop PATFinalStatesOwned_finalState*_*_*',
    'drop *_patFinalStateEvent*_*_*'
)

suffix = '' # most analyses don't need to modify the final states

if options.hzz and (options.channels == 'zz' or any(len(c)==4 for c in options.channels.split(','))):
    process.embedHZZSeq = cms.Sequence()
    # Embed matrix elements in relevant final states
    suffix = "HZZ"
    for quadFS in ['ElecElecElecElec', 
                   'ElecElecMuMu',
                   'MuMuMuMu']:
        oldName = "finalState%s"%quadFS
        embedCategoryProducer = cms.EDProducer(
            "MiniAODHZZCategoryEmbedder",
            src = cms.InputTag(oldName),
            tightLepCut = cms.string('userFloat("HZZ4lIDPassTight") > 0.5 && userFloat("HZZ4lIsoPass") > 0.5'),
            bDiscriminator = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            bDiscriminantCut = cms.double(0.89),
            )
        # give the FS collection an intermediate name, with an identifying suffix
        intermediateName = oldName + "HZZCategory"
        setattr(process, intermediateName, embedCategoryProducer)
        process.embedHZZSeq += embedCategoryProducer
        
        embedMEProducer = cms.EDProducer(
            "MiniAODHZZMEEmbedder%s"%quadFS,
            src = cms.InputTag(intermediateName),
            processes = cms.vstring("p0plus_VAJHU",
                                    "p0minus_VAJHU",
                                    "Dgg10_VAMCFM",
                                    "bkg_VAMCFM",
                                    "phjj_VAJHU",
                                    "pvbf_VAJHU",
                                    ),
            fsrLabel = cms.string("dretFSRCand"),
            )
        # give the FS collection the same name as before, but with an identifying suffix
        newName = oldName + suffix
        setattr(process, newName, embedMEProducer)
        process.embedHZZSeq += embedMEProducer
            
    process.embedHZZ = cms.Path(process.embedHZZSeq)
    process.schedule.append(process.embedHZZ)
        


_FINAL_STATE_GROUPS = {
    'zh': 'eeem, eeet, eemt, eett, emmm, emmt, mmmt, mmtt',
    'zz': 'eeee, eemm, mmmm',
    'zgg': 'eegg, mmgg',
    'llt': 'emt, mmt, eet, mmm, emm, mm, ee, em',
    'zg': 'mmg, eeg',
    'zgxtra': 'mgg, emg, egg',
    'dqm': 'e,m,t,g,j',
    '3lep': 'eee, eem, eet, emm, emt, ett, mmm, mmt, mtt, ttt',
    '4lep': 'eeee, eeem, eeet, eemm, eemt, eett, emmm, emmt, emtt, ettt, mmmm, mmmt, mmtt, mttt, tttt',
}

# run dqm
if options.runDQM: options.channels = 'dqm'

# Generate analyzers which build the desired final states.
final_states = [x.strip() for x in options.channels.split(',')]


def order_final_state(state):
    '''
    Sorts final state objects into order expected by FSA.
    
    Sorts string of characters into ordr defined by "order." Invalid 
    characters are ignored, and a warning is pribted to stdout
    
    returns the sorted string
    '''
    order = "emtgj"
    for obj in state:
        if obj not in order:
            print "invalid Final State object "\
                "'%s' ignored" % obj
            state = state.replace(obj, "")
    return ''.join(sorted(state, key=lambda x: order.index(x)))
 
def expanded_final_states(input):
    for fs in input:
        if fs in _FINAL_STATE_GROUPS:
            for subfs in _FINAL_STATE_GROUPS[fs].split(','):
                yield subfs.strip()
        else:
            yield fs

def get_producer_suffix(state):
    '''
    Returns the suffix FSA puts on the end of produer and class names, e.g.
    "ElecElecMuTau" for final state 'eemt'.
    '''
    return ''.join(_producer_translation[obj] for obj in order_final_state(state))


if options.keepPat:
    print "Saving final physics objects instead of making ntuples"

    # Clean unwanted final states (otherwise happens at flat ntuple creation)
    from FinalStateAnalysis.NtupleTools.uniqueness_cut_generator import uniqueness_cuts

    ptCuts = parameters.get('ptCuts', {'e':'0','m':'0','t':'0','g':'0','j':'0'})
    etaCuts = parameters.get('etaCuts', {'e':'10','m':'10','t':'10','g':'10','j':'10'})
    skimCuts = parameters.get('skimCuts', [])

    if options.passThru:
        output_to_keep.append('*')
    else:
        process.finalStateCleaning = cms.Sequence()

        for fs in expanded_final_states(final_states):
            fsCuts = cms.vstring()
            for name, cut in uniqueness_cuts(fs, ptCuts, etaCuts, 
                                             skimCuts=skimCuts,
                                             hzz=options.hzz, 
                                             dblH=options.dblhMode).iteritems():
                fsCuts.append(cut)
                
            fsCleaner = cms.EDProducer(
                "PATFinalStateSelector",
                src = cms.InputTag("finalState%s%s"%(get_producer_suffix(fs), suffix)),
                cuts = fsCuts,
                )
            
            # give the producer a good name so it's easy to find the output
            cleanerName = 'cleanedFinalState%s'%get_producer_suffix(fs)
            setattr(process, cleanerName, fsCleaner)
            process.finalStateCleaning += fsCleaner
            
            # and keep the output in the final file
            output_to_keep.append('*_%s_*_*'%cleanerName)

        process.cleanFinalStates = cms.Path(process.finalStateCleaning)
        process.schedule.append(process.cleanFinalStates)

    # keep final object collections (final states and associated collections already saved)
    for name, label in fs_daughter_inputs.iteritems():
        output_to_keep.append('*_%s_*_*'%label)

    # keep important gen particles
    output_to_keep.append('*_prunedGenParticles_*_*')
    if options.keepPat >= 2:
        print "... Including packedGenParticles"
        output_to_keep.append('*_packedGenParticles_*_*')
        if options.keepPat >= 3:
            print "... Including packedPFCandidates"
            output_to_keep.append('*_packedPFCandidates_*_*')
    output_commands = cms.untracked.vstring('drop *')

    for product in output_to_keep:
        output_commands.append('keep %s'%product)

    process.out = cms.OutputModule(
        "PoolOutputModule",
        fileName=cms.untracked.string(options.outputFile),
        # Drop per-event meta data from dropped objects
        dropMetaData=cms.untracked.string("ALL"),
        outputCommands=output_commands,
        )
    
    process.save = cms.EndPath(process.out)
    process.schedule.append(process.save)
else:
    print "Building ntuple for final states: %s" % ", ".join(final_states)
    for final_state in expanded_final_states(final_states):
        if additional_fs: print 'Adding ntuple {0}'.format(final_state)
        extraJets = options.nExtraJets
        final_state = order_final_state(final_state)
        analyzer = make_ntuple(*final_state, 
                                svFit=options.svFit, 
                                dblhMode=options.dblhMode,
                                runTauSpinner=options.runTauSpinner, 
                                runMVAMET=options.runMVAMET,
                                skimCuts=options.skimCuts, 
                                suffix=suffix,
                                hzz=options.hzz, 
                                nExtraJets=extraJets, 
                                use25ns=options.use25ns,
                                isMC=options.isMC,
                                isShiftedMet=bool(options.metShift),
                                **parameters)
        add_ntuple(final_state, analyzer, process,
                   process.schedule, options.eventView, filters)
        for fs in additional_fs:
            print "Adding additional ntuple with postfix {0}".format(fs)
            analyzer = make_ntuple(*final_state,
                                    svFit=options.svFit, 
                                    dblhMode=options.dblhMode,
                                    runTauSpinner=options.runTauSpinner,
                                    runMVAMET=options.runMVAMET,
                                    skimCuts=options.skimCuts, 
                                    suffix=suffix,
                                    hzz=options.hzz, 
                                    nExtraJets=extraJets,
                                    use25ns=options.use25ns,
                                    isMC=options.isMC,
                                    isShiftedMet=bool(options.metShift),
                                    postfix=fs,
                                    **parameters)
            add_ntuple(final_state+fs, analyzer, process,
                       process.schedule, options.eventView, filters)

            

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.categories.append('FSAEventMissingProduct')
process.MessageLogger.categories.append('UndefinedPreselectionInfo')
process.MessageLogger.categories.append('GsfElectronAlgo')

# Don't go nuts if there are a lot of missing products.
process.MessageLogger.cerr.FSAEventMissingProduct = cms.untracked.PSet(
    limit=cms.untracked.int32(10)
)
# process.MessageLogger.suppresssWarning = cms.untracked.vstring("GsfElectronAlgo","UndefinedPreselectionInfo")
process.MessageLogger.cerr.GsfElectronAlgo = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
process.MessageLogger.cerr.UndefinedPreselectionInfo = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

if options.verbose:
    process.options.wantSummary = cms.untracked.bool(True)
if options.passThru:
    set_passthru(process)

if options.dump:
    print process.dumpPython()

