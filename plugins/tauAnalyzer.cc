// -*- C++ -*-
//
// Package:    RecoTauTag/tauAnalyzer
// Class:      tauAnalyzer
// 
/**\class tauAnalyzer tauAnalyzer.cc RecoTauTag/tauAnalyzer/plugins/tauAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Sat, 25 Apr 2015 15:51:58 GMT
//
//


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
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TTree.h"
#include "TauTrigMatch.h"

#include "FWCore/Framework/interface/EventSetup.h"
//
// class declaration
//

class tauAnalyzer : public edm::EDAnalyzer {
	public:
		explicit tauAnalyzer(const edm::ParameterSet&);
		~tauAnalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		edm::InputTag tauSrc_;
		edm::InputTag jetSrc_;
		edm::InputTag discriminatorSrc_;
		TTree* tree;
		std::vector<Float_t>* pts_;
		std::vector<Float_t>* etas_;
		std::vector<Int_t>* dmf_;
		std::vector<Int_t>* passDiscr_;
		std::vector<Int_t>* genMatchedJet_;
		std::vector<Int_t>* genMatchedTau_;
		std::vector<Float_t>* genMatchedPt_;
		std::vector<Float_t>* jetRefPts_;
		std::vector<Float_t>* jetPts_;
		std::vector<Float_t>* jetEtas_;
		std::vector<Int_t>* jetIDLoose_;
		std::vector<Int_t>* jetIDMed_;
		std::vector<Int_t>* jetIDTight_;
		double maxDR_;

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
tauAnalyzer::tauAnalyzer(const edm::ParameterSet& cfg)
{
	//now do what ever initialization is needed
	tauSrc_              = cfg.getParameter<edm::InputTag>("recoTau");
	jetSrc_              = cfg.getParameter<edm::InputTag>("recoJet");
	discriminatorSrc_    = cfg.getParameter<edm::InputTag>("recoTauDiscriminator");

	edm::Service<TFileService> fs;
	//ntuple additions
	tree = fs->make<TTree>("Ntuple", "Ntuple");
	pts_ = new std::vector<Float_t>();
	dmf_ = new std::vector<Int_t>();
	passDiscr_ = new std::vector<Int_t>();
	genMatchedTau_ = new std::vector<Int_t>();
	genMatchedJet_ = new std::vector<Int_t>();
	genMatchedPt_ = new std::vector<Float_t>();
	jetRefPts_ = new std::vector<Float_t>();
	jetPts_ = new std::vector<Float_t>();
	jetEtas_ = new std::vector<Float_t>();
	jetIDLoose_ = new std::vector<Int_t>();
	jetIDMed_ = new std::vector<Int_t>();
	jetIDTight_ = new std::vector<Int_t>();
	tree->Branch("pt", "std::vector<float>", &pts_);
	tree->Branch("eta", "std::vector<float>", &etas_);
	tree->Branch("dmf", "std::vector<int>", &dmf_);
	tree->Branch("passDiscr", "std::vector<int>", &passDiscr_);
	tree->Branch("genMatchedJet", "std::vector<int>", &genMatchedJet_);
	tree->Branch("genMatchedTau", "std::vector<int>", &genMatchedTau_);
	tree->Branch("genMatchedPt", "std::vector<float>", &genMatchedPt_);
	tree->Branch("jetRefPt", "std::vector<float>", &jetRefPts_);
	tree->Branch("jetPt", "std::vector<float>", &jetPts_);
	tree->Branch("jetEta", "std::vector<float>", &jetEtas_);
	tree->Branch("jetIDLoose", "std::vector<int>", &jetIDLoose_);
	tree->Branch("jetIDMed", "std::vector<int>", &jetIDMed_);
	tree->Branch("jetIDTight", "std::vector<int>", &jetIDTight_);
	maxDR_ = 0.2;
}


tauAnalyzer::~tauAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	delete pts_;
	delete etas_;
	delete dmf_;
	delete passDiscr_;
	delete genMatchedJet_;
	delete genMatchedTau_;
	delete genMatchedPt_;
	delete jetRefPts_;
	delete jetPts_;
	delete jetEtas_;
	delete jetIDLoose_;
	delete jetIDMed_;
	delete jetIDTight_;
}


//
// member functions
//
// Get collection of generator particles with status 2

std::vector<const reco::GenParticle*> getGenParticleCollection(const edm::Event& evt) {
	std::vector<const reco::GenParticle*> output;
	edm::Handle< std::vector<reco::GenParticle> > handle;
	evt.getByLabel("genParticles", handle);
	// Loop over objects in current collection
	for (size_t j = 0; j < handle->size(); ++j) {
		const reco::GenParticle& object = handle->at(j);
		//if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
		if(object.status() == 2) output.push_back(&object);
	}
	return output;
}


// Get collection of pat::taus
//
std::vector<const reco::PFTau*> getRecoCandCollections(const edm::Event& evt, const edm::InputTag& collection) {
	std::vector<const reco::PFTau*> output;
	edm::Handle<std::vector<reco::PFTau> > handle;
	evt.getByLabel(collection, handle);
	// Loop over objects in current collection
	for (size_t j = 0; j < handle->size(); ++j) {
		const reco::PFTau& object = handle->at(j);
		if(object.pt()>15. && fabs(object.eta())< 3.){
			output.push_back(&object);
		}
	}
	return output;
}
// Method to find the best match between tag tau and gen object. The best matched gen tau object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
const reco::GenParticle* findBestGenMatch(const reco::PFTau TagTauObj,
		std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
	const reco::GenParticle* output = NULL;
	double bestDeltaR = -1;
	for (size_t i = 0; i < GenPart.size(); ++i) {
		double deltaR = reco::deltaR(TagTauObj, *GenPart[i]);
		if (deltaR < maxDR) {
			if (!output || deltaR < bestDeltaR) {
				output = GenPart[i];
				bestDeltaR = deltaR;
			}
		}
	}
	return output;
}

// ------------ method called for each event  ------------
	void
tauAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& iSetup)
{
	using namespace edm;
	Handle<reco::PFTauCollection> tauObjects;
	evt.getByLabel(tauSrc_, tauObjects);

	Handle<reco::PFJetCollection> jetObjects;
	evt.getByLabel(jetSrc_, jetObjects);


	edm::Handle<reco::PFTauDiscriminator> discriminator;
	evt.getByLabel(discriminatorSrc_, discriminator);

	edm::Handle<reco::PFTauDiscriminator> DMF; 
	evt.getByLabel("hpsPFTauDiscriminationByDecayModeFinding",DMF);

	pts_->clear();
	etas_->clear();
	dmf_->clear();
	passDiscr_->clear();
	genMatchedTau_->clear();
	genMatchedJet_->clear();
	genMatchedPt_->clear();
	jetRefPts_->clear();
	jetPts_->clear();
	jetEtas_->clear();
	jetIDLoose_->clear();
	jetIDMed_->clear();
	jetIDTight_->clear();

	for (unsigned int iTau = 0; iTau<tauObjects->size() ; ++iTau){
		reco::PFTauRef tauCandidate(tauObjects, iTau);
		pts_->push_back(tauCandidate->pt());
		etas_->push_back(tauCandidate->eta());
		// check if tau candidate has passed discriminator
		if( (*DMF)[tauCandidate] > 0.5 ) dmf_->push_back(1);
		else dmf_->push_back(0);
		// check if tau candidate has passed discriminator
		if( (*discriminator)[tauCandidate] > 0.5 ){
			// do something with your candidate
			passDiscr_->push_back(1);
		}
		else{
			passDiscr_->push_back(0);
		}
	}//end tau loop 


	std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollection(evt);
	for(unsigned int iTau = 0; iTau<tauObjects->size(); iTau++){
		const reco::PFTau tauCand=tauObjects->at(iTau);
		const reco::Candidate* bestGenMatch = findBestGenMatch(tauCand,GenObjects, maxDR_) ;
		if(bestGenMatch) {
			genMatchedPt_->push_back(bestGenMatch->pt());
			if (abs(bestGenMatch->pdgId())==15){
				genMatchedTau_->push_back(1);
				genMatchedJet_->push_back(0);
			}
			else{
				genMatchedJet_->push_back(1);//Not Necessarily a jet, but certainly not a tau
				genMatchedTau_->push_back(0);
			}
		}
		else {
			genMatchedPt_->push_back(-1);
			genMatchedJet_->push_back(-1);
			genMatchedTau_->push_back(-1);
		}
		if (tauCand.jetRef().isNonnull()){
			jetRefPts_->push_back(tauCand.jetRef()->pt());
		}
		else jetRefPts_->push_back(0);
	}

	for(unsigned int iJet = 0; iJet<jetObjects->size(); iJet++){
		//std::cout<<"getting jetCand"<<std::endl;
		const reco::PFJet jetCand=jetObjects->at(iJet);
		//std::cout<<"getting jetCand Pts"<<std::endl;
		if (fabs(jetCand.eta())<2.3&&jetCand.pt()>20){ 
			jetPts_->push_back(jetCand.pt());
			jetEtas_->push_back(jetCand.pt());
			bool loose = true;
			bool medium = true;
			bool tight = true;
			if (jetCand.neutralHadronEnergyFraction() >= 0.99)
				loose = false;
			if (jetCand.neutralHadronEnergyFraction() >= 0.95)
				medium = false;
			if (jetCand.neutralHadronEnergyFraction() >= 0.90)
				tight = false;

			if (jetCand.neutralEmEnergyFraction() >= 0.99)
				loose = false;
			if (jetCand.neutralEmEnergyFraction() >= 0.95)
				medium = false;
			if (jetCand.neutralEmEnergyFraction() >= 0.90)
				tight = false;

			if (jetCand.numberOfDaughters() <= 1) { //getPFConstitutents broken in miniAOD
				loose = false;
				medium = false;
				tight = false;
			}

			if (std::abs(jetCand.eta()) < 2.4) {
				if (jetCand.chargedHadronEnergyFraction() == 0) {
					loose = false;
					medium = false;
					tight = false;
				}
				if (jetCand.chargedHadronMultiplicity() == 0) {
					loose = false;
					medium = false;
					tight = false;
				}
				if (jetCand.chargedEmEnergyFraction() >= 0.99) {
					loose = false;
					medium = false;
					tight = false;
				}
			}

			jetIDLoose_->push_back(loose);
			jetIDMed_->push_back(medium);
			jetIDTight_->push_back(tight);
		}
	}//endJet
	tree->Fill();  // create TTree

}


// ------------ method called once each job just before starting event loop  ------------
	void 
tauAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
tauAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   tauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   tauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   tauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   tauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	//  edm::ParameterSetDescription desc;
	//  desc.setUnknown();
	//  descriptions.addDefault(desc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(tauAnalyzer);
