import FWCore.ParameterSet.Config as cms
import os


from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')
options.inputFiles = '/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/06293CB4-CEFD-E411-8BE4-0CC47A13D16A.root'
options.outputFile = "MiniAOD_FR_WJets.root"
options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")

#Make the framework shutup
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')

#50 ns global tag for MC replace with 'GR_P_V56' for prompt reco. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Prompt_reconstruction_Global_Tag 
process.GlobalTag.globaltag = 'MCRUN2_74_V9A'

#how many events to run over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)
#output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

##################################################
# Main
process.byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits") 
)
process.byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byMediumCombinedIsolationDeltaBetaCorr3Hits")
)
process.byTightCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
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
process.neutralIsoPtSum= cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("neutralIsoPtSum")
)
process.puCorrPtSum= cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("puCorrPtSum")
)
process.againstMuonLoose3 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstMuonLoose3")
)
process.againstMuonTight3 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstMuonTight3")
)
process.againstElectronVLooseMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronVLooseMVA5")
)
process.againstElectronLooseMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronLooseMVA5")
)
process.againstElectronMediumMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
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

dump_file = open('dump.py','w')
dump_file.write(process.dumpPython())
