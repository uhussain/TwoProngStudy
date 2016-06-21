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
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "helpers.h"
#include <iostream>
#include <algorithm>
//
// class declaration
//

class MiniAODfakeRate_ZToEE : public edm::EDAnalyzer {
	public:
		explicit MiniAODfakeRate_ZToEE(const edm::ParameterSet&);
		~MiniAODfakeRate_ZToEE();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
                edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;
                edm::EDGetTokenT<std::vector < pat::PackedGenParticle> >packedGenToken_;
		edm::EDGetTokenT<pat::ElectronCollection> electronToken_;

		std::string tauID_;

		TTree* tree;
		Int_t tauIndex_=0;
		Int_t dmf_=0;
		Float_t tauPt_=-999;
		Float_t tauEta_=-999;
		Float_t elePt_=-999;
		Float_t eleEta_=-999;
		Int_t nvtx_=-999;
		Int_t passTau_=0;
		Int_t passEle_=0;

		double maxDR_;
};

MiniAODfakeRate_ZToEE::MiniAODfakeRate_ZToEE(const edm::ParameterSet& iConfig):
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	prunedGenToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
	packedGenToken_(consumes<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
	electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
{

	tauID_    = iConfig.getParameter<std::string>("tauID");

	edm::Service<TFileService> fs;

	tree = fs->make<TTree>("Ntuple", "Ntuple");
	tree->Branch("elePt", &elePt_,"elePt_/F");
	tree->Branch("eleEta", &eleEta_,"eleEta_/F");
	tree->Branch("tauPt", &tauPt_,"tauPt_/F");
	tree->Branch("tauEta", &tauEta_,"tauEta_/F");
	tree->Branch("dmf",&dmf_,"dmf/I");
	tree->Branch("tauIndex",&tauIndex_,"tauIndex/I");
	tree->Branch("nvtx",&nvtx_,"nvtx/I");
	tree->Branch("passTau",&passTau_,"passTau/I");
	tree->Branch("passEle",&passEle_,"passEle/I");
	maxDR_ = 0.4;
}

MiniAODfakeRate_ZToEE::~MiniAODfakeRate_ZToEE()
{
}

	void
MiniAODfakeRate_ZToEE::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	if (vertices->empty()) return; // skip the event if no PV found
	//const reco::Vertex &PV = vertices->front();

	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);

	edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByToken(prunedGenToken_, genParticles);

	edm::Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electronToken_, electrons);

	std::vector<const reco::GenParticle*> GenTaus;
	std::vector<const reco::GenParticle*> GenEles;
	std::vector<const reco::GenParticle*> GenMus;

	for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++ ){
	  if(abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
	}
	for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++ ){
	  if(abs(genParticle->pdgId()) == 11) GenEles.push_back(&(*genParticle));
	}
	for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++ ){
	  if(abs(genParticle->pdgId()) == 13) GenMus.push_back(&(*genParticle));
	}

	std::vector <const pat::Electron*> ele_vec;
	std::vector <const reco::GenParticle*> used_GenTau_vec;
        nvtx_=vertices->size();

	for (const pat::Electron &ele : *electrons) {
		passEle_ = 0;
		elePt_ = ele.pt();
		eleEta_ = ele.eta();
		if (ele.pt() > 20 and abs(ele.eta()) < 2.3){  //add more discrimnators?
			passEle_ = 1;
			tree->Fill();
		}
	}
	tauIndex_ = 0;
        for (const pat::Tau &tau : *taus) {
		if (GenTaus.size() == 0) continue;
		passTau_ = 0;
		tauPt_=tau.pt();
		tauEta_=tau.eta();
		dmf_ = tau.tauID("decayModeFinding");
		if (vertices->empty()) continue; // skip the tau if no PV found
		//std::cout << "Analyzing Tau number " << tau_position << "\n";
		std::cout << GenTaus.size();
		bool pass_discr = tau.tauID(tauID_)>.5;
                if (tau.pt() > 20 && abs(tau.eta())<2.3 && pass_discr) { // if the tau passes the critera
			passTau_ = 1;
			const reco::GenParticle* bestGenTau = findBestGenMatch(tau,GenTaus,maxDR_);
			const reco::GenParticle* bestGenEle = findBestGenMatch(tau,GenEles,maxDR_);
			const reco::GenParticle* bestGenMu = findBestGenMatch(tau,GenMus,maxDR_);
			//std::cout << "The best Gen Object is " << bestGenTau << "\n";
			//std::cout << "The size of Gen Obj is " << GenTaus.size() << "\n";
			if (bestGenTau != NULL) {
				GenTaus.erase(std::remove(GenTaus.begin(), GenTaus.end(), bestGenTau), GenTaus.end());
				continue;
			}
                        else if (bestGenEle != NULL) {
                                GenEles.erase(std::remove(GenEles.begin(), GenEles.end(), bestGenEle), GenEles.end());
				continue;
                        } // end if tau is near a gen electron
                        else if (bestGenMu != NULL) {
                                GenMus.erase(std::remove(GenMus.begin(), GenMus.end(), bestGenMu), GenMus.end());
				continue;
			} // end if tau is near a gen object
                	tree->Fill();
		}// end if tau meets pt, eta, discriminator cutoffs
		tauIndex_++;
        } // end numerator for loop
} // end analyze
//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODfakeRate_ZToEE);
