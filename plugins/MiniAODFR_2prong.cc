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

// function declarations

// class declaration
class MiniAODFR_2prong : public edm::EDAnalyzer {
	public:
		explicit MiniAODFR_2prong(const edm::ParameterSet&);
		~MiniAODFR_2prong();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;
		std::string tauID_;
		edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;
                edm::EDGetTokenT<std::vector < pat::PackedGenParticle> >packedGenToken_;

		TTree* tree;
		Float_t tauPt_;
		Float_t tauEta_;
		Float_t jetPt_;
		Float_t jetEta_;
		Int_t tauIndex_;
		Int_t nvtx_;
		bool genTauMatch_;
		Int_t decayMode_;
		double maxDR_;
		bool good_dz_;
		bool good_dr_;
		Int_t fakeJet_;
};

MiniAODFR_2prong::MiniAODFR_2prong(const edm::ParameterSet& iConfig):
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
	prunedGenToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
	packedGenToken_(consumes<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
	tauID_    = iConfig.getParameter<std::string>("tauID");
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>("Ntuple", "Ntuple");
	tree->Branch("jetPt", &jetPt_,"jetPt_/F");
	tree->Branch("jetEta", &jetEta_,"jetEta_/F");
	tree->Branch("tauIndex", &tauIndex_,"tauIndex_/I");
	tree->Branch("nvtx",&nvtx_,"nvtx_/I");
	tree->Branch("decayMode",&decayMode_,"decayMode_/I");
	tree->Branch("fakeJet",&fakeJet_,"fakeJet_/I");
}

MiniAODFR_2prong::~MiniAODFR_2prong()
{
}

	void

MiniAODFR_2prong::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	nvtx_=vertices->size();
//	const reco::Vertex &PV = vertices->front();
	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);
	edm::Handle<pat::JetCollection> jets;
	iEvent.getByToken(jetToken_, jets);
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
		
        //if (GenTaus.size()>0) return; // skip event if real tau! only look for fake taus 
	nvtx_=vertices->size();

	for (const pat::Jet &jet : *jets) {  //looping through jets
		genTauMatch_ = 0;
		fakeJet_ = 0;
		jetPt_ = jet.pt();
		jetEta_ = jet.eta();
		decayMode_ = -1; 
		if (jetPt_ > 20 && TMath::Abs(jetEta_) < 2.3){  //jets must meet these criteria
			for (size_t i=0; i < GenTaus.size(); i++){	//Loop through all generated taus to check for match
				reco::Candidate::LorentzVector p4_vis = GetVisibleP4(GenTaus[i]);
				if (reco::deltaR(jet.eta(),jet.phi(),p4_vis.eta(),p4_vis.phi()) < 0.3 && p4_vis.pt() > 20.0 && TMath::Abs(p4_vis.eta())<2.3 && isHadronic(GenTaus[i])){
					genTauMatch_ = 1;
					break;	
				}
			}
			if (genTauMatch_) continue; //Skip this jet if it matches a generator tau
			tauIndex_ = 0;
			for (const pat::Tau &tau : *taus) {	//loop through all reco taus
				tauPt_ = tau.pt();
				tauEta_ = tau.eta();
				bool fake_discr = tau.tauID(tauID_)>.5; //does tau fake additional discriminator?
				bool istwoprong =(tau.decayMode()== 5 || tau.decayMode()==6);
				if (tau.tauID("decayModeFindingNewDMs")>0.5  && TMath::Abs(tauEta_)<2.3 && fake_discr && !istwoprong && deltaR(jet,tau)<0.4 && !genTauMatch_) { // if the tau passes the criteria  tauPt_>20?
					decayMode_ = tau.decayMode();
					fakeJet_ = 1; //if tau passes set the jet output to "fake". e.g. this jet faked a tau
					break;	//end tau matching search once we know a tau was faked by an electron
				}
				tauIndex_++;//which number tau was faked (good to know if highest pt or second highest pt tau, etc.)
			}//end matching tau to ele
			if (fakeJet_ == 0) tauIndex_ = -1;
			tree->Fill(); //fill gen. electron tree 
		}
	}		
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODFR_2prong);
