// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "helpers.h"

//Checks if the reconstructed particle is a neutrino (returns 1 if it is a neutrino)
bool isNeutrino(const reco::Candidate* daughter)
{
  return ( TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}

//Gets visible 4-momentum of a particle's daughter...recursively calls itself until all decay branches are searched
//See GetVisibleP4 function for more detail
reco::Candidate::LorentzVector GetDaughterVisibleP4(const reco::Candidate* daughter){
	reco::Candidate::LorentzVector p4_vis(0,0,0,0);	
 	for(size_t j = 0; j < daughter->numberOfDaughters(); ++j){
 		if (!isNeutrino(daughter->daughter(j)) && daughter->daughter(j)->status() == 1){
 			p4_vis += daughter->daughter(j)->p4();
		}
		if (daughter->daughter(j)->status() == 2){
			p4_vis += GetDaughterVisibleP4(daughter->daughter(j));
		}		
 	}
	return p4_vis;
}

//Gets visible 4-momentum of a particle (tau)
reco::Candidate::LorentzVector GetVisibleP4(const reco::GenParticle* tau){
	reco::Candidate::LorentzVector p4_vis(0,0,0,0);	
 	for(size_t j = 0; j < tau->numberOfDaughters(); ++j){  //looping through first level of decay products 
 		if (!isNeutrino(tau->daughter(j)) && tau->daughter(j)->status() == 1){  //status=1 means no further decay for this daughter
 			p4_vis += tau->daughter(j)->p4();
		} 
		if (tau->daughter(j)->status() == 2){  //status=2 means this daughter decays 
			p4_vis += GetDaughterVisibleP4(tau->daughter(j)); //we check the daughter for visible decay products
		}		
 	}
	return p4_vis;
}

//Checks if the generator-level particle (tau) decays to hadrons (returns 1 if it decays hadronically)
bool isHadronic(const reco::GenParticle* tau){
	bool isHadronic = 1;
	for(size_t j = 0; j < tau->numberOfDaughters(); ++j){ //Loop through daughters of gen. tau
		if (TMath::Abs(tau->daughter(j)->pdgId()) == 11 || TMath::Abs(tau->daughter(j)->pdgId()) == 13){ //Check if the daughter is another lepton
			isHadronic = 0;
		}
	}
	return isHadronic;
}

reco::PFJetRef getJetRef(const reco::PFTau& tau) {
	if (tau.jetRef().isNonnull())
		return tau.jetRef();
	else if (tau.pfTauTagInfoRef()->pfjetRef().isNonnull())
		return tau.pfTauTagInfoRef()->pfjetRef();
	else throw cms::Exception("cant find jet ref");
}

reco::PFJetRef getJetRef(const pat::Tau tau) {
        if (tau.pfJetRef().isNonnull())
                return tau.pfJetRef();
        else throw cms::Exception("cant find jet ref");
}

bool genMatchingMiniAOD(const pat::Tau tau, std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
	bool tau_match=false;
	for (size_t i = 0; i < GenPart.size(); ++i) {
		if (abs(GenPart[i]->pdgId())==15){
			double deltaR = reco::deltaR(tau, *GenPart[i]);
			if (deltaR < maxDR) {
				tau_match=true;
			}
		}
	}
	return tau_match;
}

std::vector<const reco::GenParticle*> getGenParticleCollectionMiniAOD(const edm::Event& evt) {
	std::vector<const reco::GenParticle*> output;
	edm::Handle< std::vector<reco::GenParticle> > handle;
	evt.getByLabel("prunedGenParticles", handle);
	// Loop over objects in current collection
	for (size_t j = 0; j < handle->size(); ++j) {
		const reco::GenParticle& object = handle->at(j);
		if(abs(object.pdgId()) == 15) output.push_back(&object);
	}
	return output;
}

std::vector<const reco::GenParticle*> getGenEleCollectionMiniAOD(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("prunedGenParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                if(abs(object.pdgId()) == 11) output.push_back(&object);
        }
        return output;
}

std::vector<const reco::GenParticle*> getGenMuCollectionMiniAOD(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("prunedGenParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                if(abs(object.pdgId()) == 13) output.push_back(&object);
        }
        return output;
}

// Get collection of generator particles with status 2
std::vector<const reco::GenParticle*> getGenParticleCollection(const edm::Event& evt) {
	std::vector<const reco::GenParticle*> output;
	edm::Handle< std::vector<reco::GenParticle> > handle;
	evt.getByLabel("genParticles", handle);
	// Loop over objects in current collection
	for (size_t j = 0; j < handle->size(); ++j) {
		const reco::GenParticle& object = handle->at(j);
		//if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
		if(abs(object.pdgId()) == 15) output.push_back(&object);
	}
	return output;
}

std::vector<const reco::GenParticle*> getGenEleCollection(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("genParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                //if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
                if(abs(object.pdgId()) == 11) output.push_back(&object);
        }
        return output;
}

std::vector<const reco::GenParticle*> getGenMuCollection(const edm::Event& evt) {
        std::vector<const reco::GenParticle*> output;
        edm::Handle< std::vector<reco::GenParticle> > handle;
        evt.getByLabel("genParticles", handle);
        // Loop over objects in current collection
        for (size_t j = 0; j < handle->size(); ++j) {
                const reco::GenParticle& object = handle->at(j);
                //if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
                if(abs(object.pdgId()) == 13) output.push_back(&object);
        }
        return output;
}


// Method to find the best match between tag tau and gen object. The best matched gen tau object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
//const reco::GenParticle* findBestGenMatch1(const reco::PFTau TagTauObj,
const reco::GenParticle* findBestGenMatch(const reco::PFTau& tauObj,
		std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
	const reco::GenParticle* output = NULL;
	double bestDeltaR = maxDR;
	for (size_t i = 0; i < GenPart.size(); ++i) {
		double deltaR = reco::deltaR(tauObj, *GenPart[i]);
		if (deltaR < maxDR) {
			if (deltaR < bestDeltaR) {
				output = GenPart[i];
				bestDeltaR = deltaR;
			}
		}
	}
	return output;
}

const reco::GenParticle* findBestGenMatch(const pat::Tau& tauObj,
                std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
        const reco::GenParticle* output = NULL;
        double bestDeltaR = maxDR;
        for (size_t i = 0; i < GenPart.size(); ++i) {
                double deltaR = reco::deltaR(tauObj, *GenPart[i]);
                if (deltaR < maxDR) {
                        if (deltaR < bestDeltaR) {
                                output = GenPart[i];
                                bestDeltaR = deltaR;
                        }
                }
        }
        return output;
}

int findBestGenMatchIndex(const pat::Tau& tauObj,
                std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
        double bestDeltaR = maxDR;
	int index = -1;
        for (size_t i = 0; i < GenPart.size(); ++i) {
                double deltaR = reco::deltaR(tauObj, *GenPart[i]);
                if (deltaR < maxDR) {
                        if (deltaR < bestDeltaR) {
                                bestDeltaR = deltaR;
                        }
                }
        }
	return index;
}

const pat::Jet* findBestJetMatch(const pat::Tau& tauObj,
                std::vector<const pat::Jet*>& jet_denom_vec, double maxDR) {
        const pat::Jet* output = NULL;
        double bestDeltaR = maxDR;
        for (size_t i = 0; i < jet_denom_vec.size(); ++i) {
                double deltaR = reco::deltaR(tauObj, *jet_denom_vec[i]);
                if (deltaR < maxDR) {
                        if (deltaR < bestDeltaR) {
                                output = jet_denom_vec[i];
                                bestDeltaR = deltaR;
                        }
                }
        }
        return output;
}

bool isLooseJet(const reco::PFJet jet){
	bool loose = true;
	if (jet.neutralHadronEnergyFraction() >= 0.99) loose = false;
	if (jet.neutralEmEnergyFraction() >= 0.99) loose = false;
	if (jet.numberOfDaughters() <= 1) loose = false; //getPFConstitutents broken in miniAOD
	if (std::abs(jet.eta()) < 2.4) {
		if (jet.chargedHadronEnergyFraction() == 0) loose = false;
		if (jet.chargedHadronMultiplicity() == 0) loose = false;
		if (jet.chargedEmEnergyFraction() >= 0.99) loose = false;
	}
	return loose;
}
bool isMediumJet(const reco::PFJet jet){
	bool medium = true;
	if (jet.neutralHadronEnergyFraction() >= 0.95) medium = false;
	if (jet.neutralEmEnergyFraction() >= 0.95) medium = false;
	if (jet.numberOfDaughters() <= 1) medium = false; //getPFConstitutents broken in miniAOD
	if (std::abs(jet.eta()) < 2.4) {
		if (jet.chargedHadronEnergyFraction() == 0) medium = false;
		if (jet.chargedHadronMultiplicity() == 0) medium = false;
		if (jet.chargedEmEnergyFraction() >= 0.99) medium = false;
	}
	return medium;
}

bool isTightJet(const reco::PFJet jet){
	bool tight = true;
	if (jet.neutralHadronEnergyFraction() >= 0.90) tight = false;
	if (jet.neutralEmEnergyFraction() >= 0.90) tight = false;
	if (jet.numberOfDaughters() <= 1) tight = false; //getPFConstitutents broken in miniAOD
	if (std::abs(jet.eta()) < 2.4) {
		if (jet.chargedHadronEnergyFraction() == 0) tight = false;
		if (jet.chargedHadronMultiplicity() == 0) tight = false;
		if (jet.chargedEmEnergyFraction() >= 0.99) tight = false;
	}
	return tight;
}

bool isLooseJet(const pat::Jet jet){
        bool loose = true;
        if (jet.neutralHadronEnergyFraction() >= 0.99) loose = false;
        if (jet.neutralEmEnergyFraction() >= 0.99) loose = false;
        if (jet.numberOfDaughters() <= 1) loose = false; //getPFConstitutents broken in miniAOD
        if (std::abs(jet.eta()) < 2.4) {
                if (jet.chargedHadronEnergyFraction() == 0) loose = false;
                if (jet.chargedHadronMultiplicity() == 0) loose = false;
                if (jet.chargedEmEnergyFraction() >= 0.99) loose = false;
        }
        return loose;
}
bool isMediumJet(const pat::Jet jet){
        bool medium = true;
        if (jet.neutralHadronEnergyFraction() >= 0.95) medium = false;
        if (jet.neutralEmEnergyFraction() >= 0.95) medium = false;
        if (jet.numberOfDaughters() <= 1) medium = false; //getPFConstitutents broken in miniAOD
        if (std::abs(jet.eta()) < 2.4) {
                if (jet.chargedHadronEnergyFraction() == 0) medium = false;
                if (jet.chargedHadronMultiplicity() == 0) medium = false;
                if (jet.chargedEmEnergyFraction() >= 0.99) medium = false;
        }
        return medium;
}

bool isTightJet(const pat::Jet jet){
        bool tight = true;
        if (jet.neutralHadronEnergyFraction() >= 0.90) tight = false;
        if (jet.neutralEmEnergyFraction() >= 0.90) tight = false;
        if (jet.numberOfDaughters() <= 1) tight = false; //getPFConstitutents broken in miniAOD
        if (std::abs(jet.eta()) < 2.4) {
                if (jet.chargedHadronEnergyFraction() == 0) tight = false;
                if (jet.chargedHadronMultiplicity() == 0) tight = false;
                if (jet.chargedEmEnergyFraction() >= 0.99) tight = false;
        }
        return tight;
}

