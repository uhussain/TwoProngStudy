import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')
#with open('files_SUSYggH.txt') as f:
#    options.inputFiles = f.readlines()
options.inputFiles = 'file:/data/uhussain/signalProduction/steps/CMSSW_8_0_3_patch2/src/Zprime_mchi_5GeV_4.root'  
#options.inputFiles ='/store/mc/RunIISummer16MiniAODv2/ZJetsToNuNu_HT-200To400_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/029F6BD2-5AC9-E611-BD77-0025907DC9D6.root'
#options.inputFiles = '/store/mc/RunIISpring16MiniAODv2/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/065837E2-DA38-E611-9157-008CFA50291C.root'
options.outputFile = "Zprime_GenReco_dR02.root"
options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')


#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8')
#process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#50 ns global tag for MC replace with 'GR_P_V56' for prompt reco. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Prompt_reconstruction_Global_Tag 
from Configuration.AlCa.GlobalTag import GlobalTag
#Make sure Global Tag mathes input file type
#process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_RunIIFall15DR76_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_v6', '')

#how many events to run over
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)


##################################################
# Main
#process.byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODtwoprong",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    taus = cms.InputTag("slimmedTaus"),    
#    PFCandidates = cms.InputTag("packedPFCandidates"),
#    tauID = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits"), 
#    packed = cms.InputTag("packedGenParticles"),
#    pruned = cms.InputTag("prunedGenParticles")
#)
#process.byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODtwoprong",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    taus = cms.InputTag("slimmedTaus"),
#    PFCandidates = cms.InputTag("packedPFCandidates"),
#    tauID = cms.string("byMediumCombinedIsolationDeltaBetaCorr3Hits"),
#    packed = cms.InputTag("packedGenParticles"),
#    pruned = cms.InputTag("prunedGenParticles")
#)
#process.byTightCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODtwoprong",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    taus = cms.InputTag("slimmedTaus"),
#    PFCandidates = cms.InputTag("packedPFCandidates"),
#    tracks = cms.InputTag("tracks"),
#    tauID = cms.string("byTightCombinedIsolationDeltaBetaCorr3Hits"),
#    packed = cms.InputTag("packedGenParticles"),
#    pruned = cms.InputTag("prunedGenParticles")
#)
process.Reco = cms.EDAnalyzer("MiniAODtwoprong",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    jets = cms.InputTag("slimmedJets"),
    tracks = cms.InputTag("tracks"),
    tauID = cms.string("decayModeFindingNewDMs"),
    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles")
)

process.GenLevel = cms.EDAnalyzer("GenTwoProng",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    jets = cms.InputTag("slimmedGenJets"),
    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles")
)
###################################################
#Global sequence

process.p = cms.Path(
             process.Reco*
             process.GenLevel
                     )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)
