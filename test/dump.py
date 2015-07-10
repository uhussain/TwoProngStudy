import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeProducerFromMiniAOD")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v3/70000/44B5E6C2-3B08-E511-A88B-0025907253D2.root')
)
process.againstElectronLooseMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('againstElectronLooseMVA5'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.againstElectronMediumMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('againstElectronMediumMVA5'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.againstElectronVLooseMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('againstElectronVLooseMVA5'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.againstMuonLoose3 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('againstMuonLoose3'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.againstMuonTight3 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('againstMuonTight3'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('byLooseCombinedIsolationDeltaBetaCorr3Hits'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('byMediumCombinedIsolationDeltaBetaCorr3Hits'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.byTightCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('byTightCombinedIsolationDeltaBetaCorr3Hits'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.neutralIsoPtSum = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('neutralIsoPtSum'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.puCorrPtSum = cms.EDAnalyzer("MiniAODfakeRate_alt",
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string('puCorrPtSum'),
    taus = cms.InputTag("slimmedTaus"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.p = cms.Path(process.byLooseCombinedIsolationDeltaBetaCorr3Hits+process.byMediumCombinedIsolationDeltaBetaCorr3Hits+process.byTightCombinedIsolationDeltaBetaCorr3Hits+process.neutralIsoPtSum+process.puCorrPtSum+process.againstMuonLoose3+process.againstMuonTight3+process.againstElectronVLooseMVA5+process.againstElectronLooseMVA5+process.againstElectronMediumMVA5)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(100)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MiniAOD_FR_QCD.root')
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

