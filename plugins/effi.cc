// -*- C++ -*-
//
// Package:    RecoTauTag/effi
// Class:      effi
// 
/**\class effi effi.cc RecoTauTag/effi/plugins/effi.cc

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
#include "TauTrigMatch.h"
#include "helpers.h"

#include "FWCore/Framework/interface/EventSetup.h"
//
// class declaration
//

class effi : public edm::EDAnalyzer {
	public:
		explicit effi(const edm::ParameterSet&);
		~effi();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		edm::InputTag tauSrc_;
		edm::InputTag jetSrc_;
		edm::InputTag discriminatorSrc_;
		TTree* treeEF;
		Float_t tauPt_;
		Float_t tauEta_;
		Int_t tauIndex_;
		Int_t dmf_;
		Int_t passDiscr_;
		Int_t genMatchedTau_;
		Float_t genMatchedPt_;
		Float_t jetRefPts_;
		Float_t jetRefEtas_;
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
effi::effi(const edm::ParameterSet& cfg)
{
	//now do what ever initialization is needed
	tauSrc_              = cfg.getParameter<edm::InputTag>("recoTau");
	jetSrc_              = cfg.getParameter<edm::InputTag>("recoJet");
	discriminatorSrc_    = cfg.getParameter<edm::InputTag>("recoTauDiscriminator");

	edm::Service<TFileService> fs;
	//ntuple additions
	treeEF = fs->make<TTree>("Ntuple", "Ntuple");
	treeEF->Branch("tauPt", &tauPt_,"tauPt_/F");
	treeEF->Branch("tauEta", &tauEta_,"tauEta_/F");
	treeEF->Branch("tauIndex", &tauIndex_,"tauIndex_/i");
	treeEF->Branch("dmf", &dmf_,"dmf_/i");
	treeEF->Branch("passDiscr", &passDiscr_,"passDiscr_/i");
	treeEF->Branch("genMatchedTau", &genMatchedTau_,"genMatchedTau_/i");
	treeEF->Branch("genMatchedPt", &genMatchedPt_,"genMatchedPt_/F");
	treeEF->Branch("jetRefPt", &jetRefPts_,"jetRefPts_/F");
	treeEF->Branch("jetRefEta", &jetRefEtas_,"jetRefEtas_/F");
	maxDR_ = 0.3;
}


effi::~effi()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
// ------------ method called for each event  ------------
	void
effi::analyze(const edm::Event& evt, const edm::EventSetup& iSetup)
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

	int tau_position=0;
	std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollection(evt);
	for (unsigned int iTau = 0; iTau<tauObjects->size() ; ++iTau){
		reco::PFTauRef tauCandidate(tauObjects, iTau);
		const reco::PFTau tauCand=tauObjects->at(iTau);
		tauPt_=tauCandidate->pt();
		tauEta_=tauCandidate->eta();

		// check if tau candidate has passed discriminator
		if( (*DMF)[tauCandidate] > 0.5 ) dmf_=1;
		else dmf_=0;
		// check if tau candidate has passed discriminator
		if( (*discriminator)[tauCandidate] > 0.5 ){passDiscr_=1;}
		else{passDiscr_=0;}

		const reco::Candidate* bestGenMatch = findBestGenMatch(*tauCandidate,GenObjects, maxDR_) ;
		if(bestGenMatch) {
			genMatchedPt_=bestGenMatch->pt();
			if (abs(bestGenMatch->pdgId())==15){
				genMatchedTau_=1;
			}
			else{
				genMatchedTau_=0;//Not Necessarily a jet, but certainly not a tau
			}
		}
		else {
			genMatchedPt_=-1;//No Gen Match
			genMatchedTau_=-1; //No Gen Match
		}
		reco::PFJetRef jet = getJetRef(*tauCandidate);
		jetRefPts_= jet->pt();
		jetRefEtas_=jet->eta();

		tauIndex_=tau_position;
		treeEF->Fill();
		tau_position++;
	}//end tau loop (efficiency loop)

}


// ------------ method called once each job just before starting event loop  ------------
	void 
effi::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
effi::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   effi::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   effi::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   effi::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   effi::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
effi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	//  edm::ParameterSetDescription desc;
	//  desc.setUnknown();
	//  descriptions.addDefault(desc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(effi);
