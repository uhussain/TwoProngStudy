import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.inputFiles = '/store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-250_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/80000/6EBE4F55-4E03-E511-9FE6-0CC47A13D416.root'
options.outputFile = "MiniAOD_effi.root"
options.parseArguments()
process = cms.Process("TreeProducerFromMiniAOD")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.GlobalTag.globaltag = cms.string('PHYS14_25_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)


##################################################
# Main
process.byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits") 
)
process.byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byMediumCombinedIsolationDeltaBetaCorr3Hits")
)
process.byTightCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byTightCombinedIsolationDeltaBetaCorr3Hits")
)
#It tells me that byCombinedIsolationDeltaBetaCorrRaw3Hits is not in the miniAOD
#process.byCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODeffi",
 #   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #  taus = cms.InputTag("slimmedTaus"),
   # jets = cms.InputTag("slimmedJets"),
    #tauID = cms.string("byCombinedIsolationDeltaBetaCorrRaw3Hits")
#)
#process.ChargedIsoPtSum = cms.EDAnalyzer("MiniAODeffi",
 #   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #  taus = cms.InputTag("slimmedTaus"),
   # jets = cms.InputTag("slimmedJets"),
    #tauID = cms.string("chargedIsoPtSum")
#)
process.neutralIsoPtSum= cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("neutralIsoPtSum")
)
process.puCorrPtSum= cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("puCorrPtSum")
)
process.againstMuonLoose3 = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstMuonLoose3")
)
process.againstMuonTight3 = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstMuonTight3")
)
process.againstElectronVLooseMVA5 = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronVLooseMVA5")
)
process.againstElectronLooseMVA5 = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronLooseMVA5")
)
process.againstElectronMediumMVA5 = cms.EDAnalyzer("MiniAODeffi",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronMediumMVA5")
)

###################################################
#Global sequence

process.p = cms.Path(
                     process.byLooseCombinedIsolationDeltaBetaCorr3Hits*
		     process.byMediumCombinedIsolationDeltaBetaCorr3Hits*
		     process.byTightCombinedIsolationDeltaBetaCorr3Hits*
 		     #process.byCombinedIsolationDeltaBetaCorrRaw3Hits*
		     #process.chargedIsoPtSum*
		     process.neutralIsoPtSum*
	 	     process.puCorrPtSum*
		     process.againstMuonLoose3*
	 	     process.againstMuonTight3*
		     process.againstElectronVLooseMVA5*
		     process.againstElectronLooseMVA5*
		     process.againstElectronMediumMVA5
                     )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)

dump_file = open('dump.py','w')
dump_file.write(process.dumpPython())
