// -*- C++ -*-
//
// Package:    RecoTauTag/fakeRate
// Class:      fakeRate
// 
/**\class fakeRate fakeRate.cc RecoTauTag/fakeRate/plugins/fakeRate.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/


// system include files
#include <memory>
#include "TSystem.h"
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

class fakeRate : public edm::EDAnalyzer {
	public:
		explicit fakeRate(const edm::ParameterSet&);
		~fakeRate();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		edm::InputTag tauSrc_;
		edm::InputTag jetSrc_;
		edm::InputTag discriminatorSrc_;

		TTree* tree;
		Float_t tauPt_;
		Float_t tauEta_;
		Int_t tauIndex_;
		Int_t dmf_;
		Int_t passDiscr_;
		Int_t genMatchedJet_;
		Int_t  genMatchedTau_;
		Float_t genMatchedPt_;
		Float_t jetRefPts_;
		Float_t jetRefEtas_;
		Float_t jetPts_;
		Float_t jetEtas_;
		Int_t jetIDLoose_;
		Int_t jetIDMed_;
		Int_t jetIDTight_;
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
fakeRate::fakeRate(const edm::ParameterSet& cfg)
{
	//now do what ever initialization is needed
	tauSrc_              = cfg.getParameter<edm::InputTag>("recoTau");
	jetSrc_              = cfg.getParameter<edm::InputTag>("recoJet");
	discriminatorSrc_    = cfg.getParameter<edm::InputTag>("recoTauDiscriminator");

	edm::Service<TFileService> fs;
	//ntuple additions
	tree = fs->make<TTree>("Ntuple", "Ntuple");

	tree->Branch("tauPt",&tauPt_,"tauPt/F");
	tree->Branch("tauEta",&tauEta_,"tauEta/F");
	tree->Branch("tauIndex",&tauIndex_,"tauIndex/I");
	tree->Branch("dmf",&dmf_,"dmf/I");
	tree->Branch("passDiscr",&passDiscr_,"passDiscr/I");
	tree->Branch("genMatchedJet",&genMatchedJet_,"genMatchedJet/I");
	tree->Branch("genMatchedTau",&genMatchedTau_,"genMatchedTau/I");
	tree->Branch("genMatchedPt",&genMatchedPt_,"genMatchedPt/F");
	tree->Branch("jetRefPt",&jetRefPts_,"jetRefPts/F");
	tree->Branch("jetRefEta",&jetRefEtas_,"jetRefEtas/F");
	tree->Branch("jetPt",&jetPts_,"jetPts/F");
	tree->Branch("jetEta",&jetEtas_,"jetEtas/F");
	tree->Branch("jetIDLoose",jetIDLoose_,"jetIDLoose/I");
	tree->Branch("jetIDMed",jetIDMed_,"jetIDMed/I");
	tree->Branch("jetIDTight",jetIDTight_,"jetIDTight/I");
	maxDR_ = 0.3;
}


fakeRate::~fakeRate()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
// Get collection of generator particles with status 2

reco::PFJetRef getJetRef(const reco::PFTau& tau) {
	if (tau.jetRef().isNonnull())
		return tau.jetRef();
	else if (tau.pfTauTagInfoRef()->pfjetRef().isNonnull())
		return tau.pfTauTagInfoRef()->pfjetRef();
	else throw cms::Exception("cant find jet ref");
}


std::vector<const reco::GenParticle*> getGenParticleCollection2(const edm::Event& evt) {
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


// Method to find the best match between tag tau and gen object. The best matched gen tau object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
//const reco::GenParticle* findBestGenMatch2(const reco::PFTau TagTauObj,
const reco::GenParticle* findBestGenMatch2(const reco::PFTau& TagTauObj,
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



// ------------ method called for each event  ------------
	void
fakeRate::analyze(const edm::Event& evt, const edm::EventSetup& iSetup)
{
	using namespace edm;
	//logPrint ("cout") <<"========ANALYZE======";
	Handle<reco::PFTauCollection> tauObjects;
	evt.getByLabel(tauSrc_, tauObjects);

	Handle<reco::PFJetCollection> jetObjects;
	evt.getByLabel(jetSrc_, jetObjects);

	edm::Handle<reco::PFTauDiscriminator> discriminator;
	evt.getByLabel(discriminatorSrc_, discriminator);

	edm::Handle<reco::PFTauDiscriminator> DMF; 
	evt.getByLabel("hpsPFTauDiscriminationByDecayModeFinding",DMF);


	//Fake rate loop
	//
	
	tauPt_=-999;
	tauEta_=-999;
	tauIndex_=-1;
	dmf_=0;
	passDiscr_=0;
	genMatchedJet_=0;
	genMatchedTau_=0;
	genMatchedPt_=-999;
	jetRefPts_=-999;
	jetRefEtas_=-999;
	jetPts_=-999;
	jetEtas_=-999;
	jetIDLoose_=0;
	jetIDMed_=0;
	jetIDTight_=0;

	std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollection2(evt);
	for(unsigned int iJet = 0; iJet<jetObjects->size(); iJet++){
		const reco::PFJet jetCand=jetObjects->at(iJet);
		jetPts_=-999;
		jetEtas_=-999;
		jetIDLoose_=0;
		jetIDMed_=0;
		jetIDTight_=0;
		if (std::abs(jetCand.eta())<2.3&&jetCand.pt()>20){
			jetPts_=jetCand.pt();
			jetEtas_=jetCand.eta();
			jetIDLoose_=isLooseJet(jetCand);
			jetIDMed_=isMediumJet(jetCand);
			jetIDTight_=isTightJet(jetCand);
			int tau_position=0;
			int bestDR=999;
			tauPt_=-999;
			tauEta_=-999;
			tauIndex_=-1;
			passDiscr_=0;
			jetRefPts_=0;
			jetRefEtas_=0;
			dmf_=0;
			genMatchedTau_=0;
			genMatchedJet_=0;
			genMatchedPt_=0;
			for (unsigned int iTau = 0; iTau<tauObjects->size() ; ++iTau){
				reco::PFTauRef tauCandidate(tauObjects, iTau);
				const reco::PFTau tauCand=tauObjects->at(iTau);
				const reco::Candidate* bestGenMatch = findBestGenMatch2(*tauCandidate,GenObjects, maxDR_) ;
				if(bestGenMatch) {
					genMatchedPt_=bestGenMatch->pt();
					if (abs(bestGenMatch->pdgId())==15){genMatchedTau_=1;}
				}
				else { //this should never happen. Should add a exceptin statement
					genMatchedPt_=-1;//No Gen Match
					genMatchedJet_=-1; //No Gen Match
					genMatchedTau_=-1; //No Gen Match
				}//end found genPArticle match
				if (genMatchedTau_!=1){
					double deltaR = reco::deltaR(*tauCandidate,jetCand);
					//logPrint ("cout") <<"====deltaR: "<<deltaR;
					if (deltaR<maxDR_&&deltaR<bestDR){
						bestDR=deltaR;
						tauPt_=tauCandidate->pt();
						tauEta_=tauCandidate->eta();
						// check if tau candidate has passed discriminator
						if( (*DMF)[tauCandidate] > 0.5 ) dmf_=1;
						else dmf_=0;
						// check if tau candidate has passed discriminator
						if( (*discriminator)[tauCandidate] > 0.5 ){passDiscr_=1;}
						else{passDiscr_=0;}
						reco::PFJetRef jet = getJetRef(*tauCandidate);
						jetRefPts_= jet->pt();
						jetRefEtas_=jet->eta();
						genMatchedJet_=1;
						tauIndex_=tau_position;
					}//end if closest to jet			
				}//end if not gen tau matched
				tau_position++;
			}//end tau loop (efficiency loop)
			//edm::logPrint ("cout") <<"====tau fake: "<<genMatchedJet_;
			//logPrint ("cout") <<"====tau fake pt: "<<tauPt_;
			//logPrint ("cout") <<"====tau fake eta: "<<tauEta_;
			//logPrint ("cout") <<"====tau fake jetrefpt: "<<jetRefPts_;
			//logPrint ("cout") <<"====jetpt: "<<jetPts_;
			//logPrint ("cout") <<"====tau fake jetrefeta: "<<jetRefEtas_;
			//logPrint ("cout") <<"====jetrefeta: "<<jetEtas_;
			//logPrint ("cout") <<"========FILL======";

			std::cout<<"tauPt_: "<<tauPt_<<std::endl;
			std::cout<<"tauEta_: "<<tauEta_<<std::endl;
			std::cout<<"tauIndex_: "<<tauIndex_<<std::endl;
			std::cout<<"dmf_: "<<dmf_<<std::endl;
			std::cout<<"passDiscr_: "<<passDiscr_<<std::endl;
			std::cout<<"genMatchedJet_: "<<genMatchedJet_<<std::endl;
			std::cout<<"genMatchedTau_: "<<genMatchedTau_<<std::endl;
			std::cout<<"genMatchedPt_: "<<genMatchedPt_<<std::endl;
			std::cout<<"jetRefPts_: "<<jetRefPts_<<std::endl;
			std::cout<<"jetRefEtas_: "<<jetRefEtas_<<std::endl;
			std::cout<<"jetRefEtas_: "<<jetRefEtas_<<std::endl;
			std::cout<<"jetPts_: "<<jetPts_<<std::endl;
			std::cout<<"jetEtas_: "<<jetEtas_<<std::endl;
			std::cout<<"jetIDLoose_: "<<jetIDLoose_<<std::endl;
			std::cout<<"jetIDMed_: "<<jetIDMed_<<std::endl;
			std::cout<<"jetIDTight_: "<<jetIDTight_<<std::endl;

			tree->Fill();
		}//end if pt eta
	}//end jet loops

}


// ------------ method called once each job just before starting event loop  ------------
	void 
fakeRate::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
fakeRate::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   fakeRate::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   fakeRate::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   fakeRate::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   fakeRate::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
fakeRate::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	//  edm::ParameterSetDescription desc;
	//  desc.setUnknown();
	//  descriptions.addDefault(desc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(fakeRate);
