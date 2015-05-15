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
#include "helpers.h"

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
		Float_t jetPt_;
		Float_t jetEta_;
		Int_t jetIDLoose_;
		Int_t jetIDMed_;
		Int_t jetIDTight_;
		Int_t genMatchedTau_;
		Int_t isFake_;
		Int_t tauIndex_;
		Int_t dmf_;
		Int_t passDiscr_;
		Float_t tauPt_;
		Float_t tauEta_;
		Float_t jetRefEta_;
		Float_t jetRefPt_;

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
	tree->Branch("jetPt", &jetPt_,"jetPt_/F");
	tree->Branch("tauPt", &tauPt_,"tauPt_/F");
	tree->Branch("jetEta", &jetEta_,"jetEta_/F");
	tree->Branch("tauEta", &tauEta_,"tauEta_/F");
	tree->Branch("jetRefPt",&jetRefPt_,"jetRefPt/F");
	tree->Branch("jetRefEta",&jetRefEta_,"jetRefEta/F");
	tree->Branch("jetIDLoose",&jetIDLoose_,"jetIDLoose/I");
	tree->Branch("jetIDMed",&jetIDMed_,"jetIDMed/I");
	tree->Branch("jetIDTight",&jetIDTight_,"jetIDTight/I");
	tree->Branch("jetIDTight",&jetIDTight_,"jetIDTight/I");
	tree->Branch("genMatchedTau",&genMatchedTau_,"genMatchedTau/I");
	tree->Branch("dmf",&dmf_,"dmf/I");
	tree->Branch("isFake",&isFake_,"isFake/I");
	tree->Branch("tauIndex",&tauIndex_,"tauIndex/I");
	tree->Branch("passDiscr",&passDiscr_,"passDiscr/I");
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

// ------------ method called for each event  ------------
	void
fakeRate::analyze(const edm::Event& evt, const edm::EventSetup& iSetup)
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

	std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollection(evt);
	for(unsigned int iJet = 0; iJet<jetObjects->size(); iJet++){
		const reco::PFJet jetCand=jetObjects->at(iJet);
		if (std::abs(jetCand.eta())<2.3&&jetCand.pt()>20){
			jetPt_=jetCand.pt();
			jetEta_=jetCand.eta();
			jetIDLoose_=isLooseJet(jetCand);
			jetIDMed_=isMediumJet(jetCand);
			jetIDTight_=isTightJet(jetCand);

			//initialize
			int bestDR=999;//right placement?
			int tau_position=-1;
			genMatchedTau_=0;
			isFake_=0;
			dmf_=0;
			passDiscr_=0;
			tauPt_=-999;
			tauEta_=-999;
			tauIndex_=-1;
			jetRefPt_=-999;
			jetRefEta_=-999;

			for (unsigned int iTau = 0; iTau<tauObjects->size() ; ++iTau){
				reco::PFTauRef tauCandidate(tauObjects, iTau);
				const reco::Candidate* bestGenMatch = findBestGenMatch(*tauCandidate,GenObjects, maxDR_) ;
				genMatchedTau_=0;
				tau_position++;
				if(bestGenMatch) {
					if (abs(bestGenMatch->pdgId())==15){genMatchedTau_=1;}
					else {
						genMatchedTau_=0;
						double deltaR = reco::deltaR(*tauCandidate,jetCand);
						if (deltaR<maxDR_&&deltaR<bestDR){
							bestDR=deltaR;
							if( (*DMF)[tauCandidate] > 0.5 ) dmf_=1;
							else dmf_=0;
							if( (*discriminator)[tauCandidate] > 0.5 ){passDiscr_=1;}
							else{passDiscr_=0;}
							tauPt_=tauCandidate->pt();
							tauEta_=tauCandidate->eta();
							reco::PFJetRef jet = getJetRef(*tauCandidate);
							jetRefPt_= jet->pt();
							jetRefEta_=jet->eta();
							tauIndex_=tau_position;
							isFake_=1;
						}//end deltaR
					}//end not genTau
				}//end if best gen match
			}//end tau
			tree->Fill();
		}//end if
	}//end Jet
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
