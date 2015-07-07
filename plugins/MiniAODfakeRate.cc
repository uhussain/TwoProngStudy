// system include files
#include <memory>
#include <vector>
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "helpers.h"
#include <iostream>
//
// class declaration
//

class MiniAODfakeRate : public edm::EDAnalyzer {
	public:
		explicit MiniAODfakeRate(const edm::ParameterSet&);
		~MiniAODfakeRate();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;

		std::string tauID_;

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
		Int_t nvtx_;

		double maxDR_;
};

MiniAODfakeRate::MiniAODfakeRate(const edm::ParameterSet& iConfig):
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
{

	tauID_    = iConfig.getParameter<std::string>("tauID");

	edm::Service<TFileService> fs;

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
	tree->Branch("genMatchedTau",&genMatchedTau_,"genMatchedTau/I");
	tree->Branch("dmf",&dmf_,"dmf/I");
	tree->Branch("isFake",&isFake_,"isFake/I");
	tree->Branch("tauIndex",&tauIndex_,"tauIndex/I");
	tree->Branch("passDiscr",&passDiscr_,"passDiscr/I");
	tree->Branch("nvtx",&nvtx_,"nvtx/I");
	maxDR_ = 0.4;
}

MiniAODfakeRate::~MiniAODfakeRate()
{
}

	void
MiniAODfakeRate::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	if (vertices->empty()) return; // skip the event if no PV found
	const reco::Vertex &PV = vertices->front();

	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);
	edm::Handle<pat::JetCollection> jets;
	iEvent.getByToken(jetToken_, jets);

	std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollectionMiniAOD(iEvent);
	std::vector<const reco::GenParticle*> GenEles = getGenEleCollectionMiniAOD(iEvent);
	std::vector<const reco::GenParticle*> GenMus = getGenMuCollectionMiniAOD(iEvent);

	int tau_position=-1;
	genMatchedTau_=0;
	isFake_=0;
	dmf_=0;
	passDiscr_=0;
	tauPt_=-999;
	tauEta_=-999;
	tauIndex_=-1;
	jetRefPt_=-999;
	jetPt_=-999;
	jetEta_=-999;
	jetRefEta_=-999;
	jetIDLoose_=0;
	jetIDMed_=0;
	jetIDTight_=0;

	std::vector <int> tau_position_vec;
/*
	if (GenObjects.size()!=0) {
		//std::cout << "Found " << GenObjects.size() << " Gen Taus!\n"; //output to be deleted later!!!
		return;
	}
*/
/*
	if (GenEles.size()!=0) {
                std::cout << "Found " << GenEles.size() << " Gen Eles!\n"; //output to be deleted later!!!
                //return;
        }

        if (GenMus.size()!=0) {
                std::cout << "Found " << GenMus.size() << " Gen Mus!\n"; //output to be deleted later!!!
                //return;
        }
*/	
	int iJet = -1;
	bool matched = false;
	for (const pat::Jet &jet : *jets) {
		jetIDLoose_=0;
		jetIDMed_=0;
		jetIDTight_=0;
		tau_position=-1;
		genMatchedTau_=0;
		isFake_=0;
		dmf_=0;
		passDiscr_=0;
		tauPt_=-999;
		tauEta_=-999;
		tauIndex_=-1;
		jetRefPt_=-999;
		jetPt_=-999;
		jetEta_=-999;
		jetRefEta_=-999;
		iJet++;
		//std::cout << "Analyzing Jet #" << iJet << "\n"; //output to be deleted later!!!

		if (jet.pt() > 20&&jet.eta()<2.3&&(isLooseJet(jet))) {
			//std::cout << "The Jet PT is >20, jetEta<2.3, it is loose \n";

			jetPt_=jet.pt();
			jetEta_=jet.eta();
			jetIDLoose_=isLooseJet(jet);
			jetIDMed_=isMediumJet(jet);
			jetIDTight_=isTightJet(jet);
			
			tau_position = -1;
	
			genMatchedTau_=0;
			isFake_=0;
			dmf_=0;
			passDiscr_=0;
			tauPt_=-999;
			tauEta_=-999;
			tauIndex_=-1;
			jetRefPt_=-999;
			jetRefEta_=-999;
			
			int taus_matched = 0;
			for (const pat::Tau &tau : *taus) {
				if (GenObjects.size()!=0) {
					return;
                	        }
				genMatchedTau_=0;
				isFake_=0;
				dmf_=0;
				passDiscr_=0;
				tauPt_=tau.pt();
				tauEta_=tau.eta();
				jetRefPt_=-999;
				jetRefEta_=-999;
				
				tau_position++;
				matched = false;
				for(int j=0; abs(j)<abs(tau_position_vec.size()); j++)
					if (tau_position_vec[j] == tau_position) matched = true;
				if (matched == true) continue;
				tauIndex_ = tau_position;
				//std::cout << "Analyzing tau #" << tau_position << "\n";
				double deltaR = reco::deltaR(tau,jet);
                                //std::cout << "delta R is " << deltaR << "\n";
				//std::cout << "decayModeFinding is " << tau.tauID("decayModeFinding") << "\n";
				//std::cout << "discriminator value is " << tau.tauID(tauID_) << "\n";
				if (deltaR>maxDR_){
					//std::cout << "dR>maxDR, skipping tau \n";
					continue;}
				else if (/*tau.tauID("decayModeFinding")>.5&&*/deltaR<maxDR_&&tauPt_>20&&tauEta_<2.3&&tau.tauID(tauID_)>=.5){
					//std::cout << "dR<maxDR, tauPT>20, tauEta<2.3, and tau passes Discriminator... faked! \n";
					tau_position_vec.push_back(tau_position);
					tauPt_=tau.pt();
					tauEta_=tau.eta();
					isFake_ = 1;
					//pat::JetRef jetRef = getJetRef(tau);
					dmf_=tau.tauID("decayModeFinding");//dmf always==1 in miniaod
					passDiscr_ = 1;
					//jetRefPt_ = jetRef->pt();
					//jetRefEta_ = jetRef->eta();
					tauIndex_ = tau_position;
					if (genMatchingMiniAOD(tau,GenObjects,maxDR_)) {genMatchedTau_=1;}
					else {genMatchedTau_ = 0;}
					taus_matched++;
					//std::cout << "Jet number " << iJet << " matches " << taus_matched << " taus\n";
					break;
				} //end elseif
				else {continue;}
			} //end tau for loop
			//std::cout << "filling the tree \n";
			tree->Fill();
		} // end if 
	} // end jet

}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODfakeRate);
