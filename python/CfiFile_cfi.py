import FWCore.ParameterSet.Config as cms

tauAnalyzer = cms.EDAnalyzer('tauAnalyzer',
                      recoTau              = cms.InputTag("hpsPFTauProducer"),
                      recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")

)
