import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeProducerFromAOD")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')

process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'PHYS14_25_V2',conditions='TrackerAlignmentExtendedError_2011Realistic_v1_mc,TrackerAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonDTAPEObjectsExtended_v0_mc,DTAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonCSCAPEObjectsExtended_v0_mc,CSCAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalSamplesCorrelation_mc,EcalSamplesCorrelationRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseShapes_mc,EcalPulseShapesRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseCovariances_mc,EcalPulseCovariancesRcd,frontier://FrontierProd/CMS_CONDITIONS')
#process.GlobalTag.globaltag = cms.string('PHYS14_25_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# AODSIM        
        #'file:/hdfs/store/mc/Fall13dr/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/4A1577F6-4472-E311-886E-7845C4FC39E3.root'
        'file:/hdfs/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/28C4E0C1-7F6F-E411-AE20-0025905B85EE.root'
        #'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
        #'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/AODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/125A6B71-C56A-E411-9D2B-0025907609BE.root'
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    inputCommands=cms.untracked.vstring(
        'keep *'
        #'drop patTaus_*_*_*',
        #'drop *PFTau*_*_*_*'
    )
)
 
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",#dummy
        fileName = cms.untracked.string('patTuple.root'),
        # save only events passing the full path
        SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
        # save PAT Layer 1 output; you need a '*' to
        # unpack the list of commands 'patEventContent'
        dropMetaData   = cms.untracked.string('DROPPED'),
        outputCommands = cms.untracked.vstring('keep *')
)
#####################################################
  
#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")
# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute')
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
switchJetCollection(
    process,
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    jetCorrections = ( 'AK4PFchs', jetCorrections, "" ),
    outputModules = []
)


#process.patJets.addTagInfos = cms.bool(True)
#process.patJets.addBTagInfo = cms.bool(True)


#--------------------------------------------------------------------------------
# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

# switch on PAT trigger                                                                                                                      
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

process.makePatTrigger = cms.Sequence(process.patTrigger*process.patTriggerEvent)

#--------------------------------------------------------------------------------
# select "good" reconstructed vertices
#process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0"), # CV: cut >= 4 if using 'offlinePrimaryVertices',
                                                                         #         >= 7 if using 'offlinePrimaryVerticesWithBS' as input
    filter = cms.bool(False)                                          
)

process.selectedPrimaryVertexPosition = cms.EDFilter("VertexSelector",
    src = cms.InputTag('selectedPrimaryVertexQuality'),
    cut = cms.string("abs(z) < 24 & abs(position.Rho) < 2."),
    filter = cms.bool(False)                                           
)

process.selectedPrimaryVertexHighestPtTrackSum = cms.EDFilter("PATSingleVertexSelector",
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag('selectedPrimaryVertexPosition'),
    filter = cms.bool(False)                                                    
)

process.selectPrimaryVertex = cms.Sequence(
    process.selectedPrimaryVertexQuality
   * process.selectedPrimaryVertexPosition
   * process.selectedPrimaryVertexHighestPtTrackSum
)
##################################################
# Main

process.tauAnalyzer = cms.EDAnalyzer('tauAnalyzer',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")
)

###################################################
#Global sequence

process.p = cms.Path(process.selectPrimaryVertex *
                     process.makePatMuons*
                     process.makePatElectrons*
                     process.makePatTaus*
                     process.makePatJets*
                     process.makePatMETs*
                     process.makePatTrigger*
                     process.tauAnalyzer
                     )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output.root")
)

# Let it run
process.pathEnd = cms.EndPath(
#        process.out
)

process.out.outputCommands = cms.untracked.vstring('keep *')

process.schedule = cms.Schedule(process.p)






