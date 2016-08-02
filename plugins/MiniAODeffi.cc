// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "helpers.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "iostream"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

// function declarations

// class declaration
class MiniAODeffi : public edm::EDAnalyzer {
	public:
		explicit MiniAODeffi(const edm::ParameterSet&);
		~MiniAODeffi();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;
		edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
		std::string tauID_;
		edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;
                edm::EDGetTokenT<std::vector < pat::PackedGenParticle> >packedGenToken_;

		TTree* tree;
		Float_t tauPt_;
		Float_t tauEta_;
		Float_t tauMass_;
		Int_t tauIndex_;
		Int_t nvtx_;
		Int_t dmf_;
		Int_t goodReco_;
		Int_t genTauMatch_;
		double maxDR_;
		bool good_dz_;
		bool good_dr_;
};

MiniAODeffi::MiniAODeffi(const edm::ParameterSet& iConfig):
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
        electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
	prunedGenToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
	packedGenToken_(consumes<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
	tauID_    = iConfig.getParameter<std::string>("tauID");
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>("Ntuple", "Ntuple");
	tree->Branch("tauPt", &tauPt_,"tauPt_/F");
	tree->Branch("tauEta", &tauEta_,"tauEta_/F");
	tree->Branch("tauIndex", &tauIndex_,"tauIndex_/I");
	tree->Branch("nvtx",&nvtx_,"nvtx_/I");
	tree->Branch("dmf",&dmf_,"dmf_/I");
	tree->Branch("goodReco",&goodReco_,"goodReco_/I");
	tree->Branch("tauMass",&tauMass_,"tauMass_/I");
}

MiniAODeffi::~MiniAODeffi()
{
}

	void

MiniAODeffi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	nvtx_=vertices->size();
//	const reco::Vertex &PV = vertices->front();
	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);
 	edm::Handle<std::vector<reco::GenParticle> > genParticles;
 	iEvent.getByToken(prunedGenToken_, genParticles);

 	std::vector<const reco::GenParticle*> GenTaus;
 	std::vector<const reco::GenParticle*> GenEles;
 	std::vector<const reco::GenParticle*> GenMus;
 	//Place generated leptons into separate lists
	for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
		if(TMath::Abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
		if(TMath::Abs(genParticle->pdgId()) == 11) GenEles.push_back(&(*genParticle));
		if(TMath::Abs(genParticle->pdgId()) == 13) GenMus.push_back(&(*genParticle));
	}
	
	if (GenEles.size() > 0 || GenMus.size() > 0 || GenTaus.size() == 0) return; //only want gen taus
	tauIndex_ = 0;
	goodReco_ = -1;
	for(const pat::Tau &tau : *taus){	//Loop through all reconstructed taus
		genTauMatch_ = 0;	//Assume this tau does not match a generated tau
		dmf_ = tau.tauID("decayModeFinding");
		if (!(tau.pt() > 20.0 && TMath::Abs(tau.eta())<2.3 && dmf_>0.5 && tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"))) continue; 
		tauPt_ = tau.pt();
		tauEta_ = tau.eta();
		tauMass_ = tau.mass();
		for (size_t i=0; i < GenTaus.size(); i++){	//Loop through all generated taus to check for match
			reco::Candidate::LorentzVector p4_vis = GetVisibleP4(GenTaus[i]);
			if (reco::deltaR(tau.eta(),tau.phi(),p4_vis.eta(),p4_vis.phi()) < 0.3 && p4_vis.pt() > 20.0 && TMath::Abs(p4_vis.eta())<2.3 && isHadronic(GenTaus[i])){
				genTauMatch_ = 1;
				break;	
			}
		}
		if (genTauMatch_ == 1) { //Tau must meet denominator requirements
			dmf_ = tau.decayMode();  //switch value from "decayModeFinding" so we can look at the actual decay mode 
			goodReco_ = tau.tauID(tauID_) >0.5; //Lepton discriminant for numerator
			tree->Fill(); 
		}
		++tauIndex_;
	}		
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODeffi);
