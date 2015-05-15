import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeProducerFromMiniAOD")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.GlobalTag.globaltag = cms.string('PHYS14_25_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3C05111C-8C6F-E411-A93E-7845C4F91450.root'
    )
)
 

##################################################
# Main
process.demo = cms.EDAnalyzer("MiniAODtester",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits") 
)

###################################################
#Global sequence

process.p = cms.Path(
                     process.demo
                     )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output.root")
)

dump_file = open('dump.py','w')
dump_file.write(process.dumpPython())




