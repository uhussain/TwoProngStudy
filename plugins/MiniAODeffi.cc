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

//
// class declaration
//

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
		Int_t tauIndex_;
		Int_t nvtx_;
		Int_t dmf_;
		Int_t goodReco_;
		double maxDR_;
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
	maxDR_ = 0.3;


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
	const reco::Vertex &PV = vertices->front();
	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);
 	edm::Handle<std::vector<reco::GenParticle> > genParticles;
 	iEvent.getByToken(prunedGenToken_, genParticles);

 	std::vector<const reco::GenParticle*> GenTaus;
 	std::vector<const reco::GenParticle*> GenEles;
 	std::vector<const reco::GenParticle*> GenMus;
 	//add code to make this into GenTaus/GenEles/GenMus
	for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++ ){
	  if(abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
	  if(abs(genParticle->pdgId()) == 11) GenEles.push_back(&(*genParticle));
	  if(abs(genParticle->pdgId()) == 13) GenMus.push_back(&(*genParticle));
	}	
	
	for (size_t i = 0; i < GenTaus.size(); i++) { //Loop through all generated taus
		for(size_t j = 0; j < GenTaus[i]->numberOfDaughters(); ++j){ //Loop through daughters of gen. tau
			if (abs(GenTaus[i]->daughter(j)->pdgId()) == 11 || abs(GenTaus[i]->daughter(j)->pdgId()) == 13){ //Check if the daughter is another lepton
				goto nothadronictau;                           //Skip over non-hadronic tau decays
			}
		}
		if (GenTaus[i]->pt() > 20 && abs(GenTaus[i]->eta())<2.3) {
			tauPt_= GenTaus[i]->pt();
			tauEta_= GenTaus[i]->eta();
			tauIndex_ = 0;
			for (const pat::Tau &tau : *taus) {
				goodReco_=0; //Assume each tau is not reconstructed from the generated tau
				dmf_=tau.tauID("decayModeFinding"); // this is the old DMF; strictly tighter than new DMF
				double deltaR = reco::deltaR(tau, *GenTaus[i]);
				if (tau.pt() > 20 && abs(tau.eta())<2.3 && tau.tauID(tauID_)>.5 && dmf_ > 0.5 && abs(tau.vertex().z() - PV.z())<.2 && deltaR<maxDR_) {
					goodReco_=1;
					break; //Break the tau loop once we find a correctly reconstructed tau
				} // end if tau passes criteria
				tauIndex_++;
			} // end tau for loop
			tree->Fill();
			nothadronictau:;
		} //end if gen tau matches critera
	} //end gen tau for loop   
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODeffi);
