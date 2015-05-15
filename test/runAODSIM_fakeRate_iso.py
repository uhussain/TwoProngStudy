import FWCore.ParameterSet.Config as cms
import os
#########Var Parsin##########
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.inputFiles = 'file:/hdfs/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/28C4E0C1-7F6F-E411-AE20-0025905B85EE.root'
options.outputFile = "testTau_fakeRate.root"

options.register ('isoDBCone',  0.8, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                  "isoConeSizeForDeltaBeta")

options.parseArguments()
print '========Tau Isolation Cone Configuration======='
print 'isoConeSizeForDeltaBeta =   ',options.isoDBCone,''


process = cms.Process("TreeProducerFromAOD")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')

process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'PHYS14_25_V2',conditions='TrackerAlignmentExtendedError_2011Realistic_v1_mc,TrackerAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonDTAPEObjectsExtended_v0_mc,DTAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonCSCAPEObjectsExtended_v0_mc,CSCAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalSamplesCorrelation_mc,EcalSamplesCorrelationRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseShapes_mc,EcalPulseShapesRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseCovariances_mc,EcalPulseCovariancesRcd,frontier://FrontierProd/CMS_CONDITIONS')
#process.GlobalTag.globaltag = cms.string('PHYS14_25_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    inputCommands=cms.untracked.vstring(
        'keep *'
        #'drop patTaus_*_*_*',
        #'drop *PFTau*_*_*_*'
    )
)

process.TFileService = cms.Service(
   "TFileService",
   fileName = cms.string(options.outputFile)
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
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff") #loading the configuration
# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr.isoConeSizeForDeltaBeta = cms.double(options.isoDBCone)
process.hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr.isoConeSizeForDeltaBeta = cms.double(options.isoDBCone)
process.hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr.isoConeSizeForDeltaBeta = cms.double(options.isoDBCone)
process.hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr.isoConeSizeForDeltaBeta = cms.double(options.isoDBCone)
process.hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits.isoConeSizeForDeltaBeta = cms.double(options.isoDBCone)
process.hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits.isoConeSizeForDeltaBeta = cms.double(options.isoDBCone)
process.hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits.isoConeSizeForDeltaBeta = cms.double(options.isoDBCone)

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

#testing
process.byVLooseIsolation = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation")
)
process.byLooseIsolation = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")
)

#dB
process.byVLooseCombinedIsolationDBSumPtCorr = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr")
)
process.byLooseCombinedIsolationDBSumPtCorr = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr")
)
process.byMediumCombinedIsolationDBSumPtCorr = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr")
)
process.byTightCombinedIsolationDBSumPtCorr = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr")
)

#3Hits
process.byLooseCombinedIsolationDBSumPtCorr3Hits = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits")
)
process.byMediumCombinedIsolationDBSumPtCorr3Hits = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits")
)
process.byTightCombinedIsolationDBSumPtCorr3Hits = cms.EDAnalyzer('fakeRate',
                                     recoTau              = cms.InputTag("hpsPFTauProducer"),
			      	     recoJet              = cms.InputTag("ak4PFJetsCHS"),
                                     recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits")
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
		     process.PFTau*
                     process.byVLooseIsolation*
                     process.byLooseIsolation*
                     process.byVLooseCombinedIsolationDBSumPtCorr*
                     process.byLooseCombinedIsolationDBSumPtCorr*
                     process.byMediumCombinedIsolationDBSumPtCorr*
                     process.byTightCombinedIsolationDBSumPtCorr*
                     process.byLooseCombinedIsolationDBSumPtCorr3Hits*
                     process.byMediumCombinedIsolationDBSumPtCorr3Hits*
                     process.byTightCombinedIsolationDBSumPtCorr3Hits
                     )

# Let it run
process.pathEnd = cms.EndPath(
#        process.out
)

process.out.outputCommands = cms.untracked.vstring('keep *')

process.schedule = cms.Schedule(process.p)

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())




