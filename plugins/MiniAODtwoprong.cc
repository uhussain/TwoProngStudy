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
#include <iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include <stdio.h>
#include <fstream>
// function declarations
std::ofstream myfile ("myfile.txt");


// class declaration
class MiniAODtwoprong : public edm::EDAnalyzer {
    public:
        explicit MiniAODtwoprong(const edm::ParameterSet&);
        ~MiniAODtwoprong();

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken;
        edm::EDGetTokenT<pat::TauCollection> tauToken;
        edm::EDGetTokenT<std::vector <pat::PackedCandidate>> PFCandidateToken;
        std::string tauID;
        edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken;
        edm::EDGetTokenT<std::vector < pat::PackedGenParticle> >packedGenToken;

        TTree* tree;
        Int_t nvtx;
        Float_t recoTauPt;
        Float_t recoTauEta;
        Float_t recoTauMass;
        Int_t recoDecayMode;
        Float_t tauPt;
        Float_t tauEta;
        Float_t tauMass;
        Int_t decayMode;
        Int_t leadPDG;
        Float_t pT1;
        Float_t dR1;
        Float_t pT2;
        Float_t dR2;
        Float_t pT3;
        Float_t dR3;
        Float_t dxyErr1;
        Float_t dxyErr2;
        Float_t dxyErr3;
        Float_t dzErr1;
        Float_t dzErr2;
        Float_t dzErr3;
        Float_t numHits1;
        Float_t numHits2;
        Float_t numHits3;
        Float_t numPixHits1;
        Float_t numPixHits2;
        Float_t numPixHits3;

};

MiniAODtwoprong::MiniAODtwoprong(const edm::ParameterSet& iConfig):
    vtxToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    tauToken(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    PFCandidateToken(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
    prunedGenToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken(consumes<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
    tauID = iConfig.getParameter<std::string>("tauID");
    edm::Service<TFileService> fs;

    tree = fs->make<TTree>("Ntuple", "Ntuple");
    tree->Branch("nvtx",&nvtx,"nvtx/I");
    tree->Branch("recoTauPt", &recoTauPt,"recoTauPt/F");
    tree->Branch("recoTauEta", &recoTauEta,"recoTauEta/F");
    tree->Branch("recoTauMass",&recoTauMass,"recoTauMass/F");
    tree->Branch("recoDecayMode",&recoDecayMode,"recoDecayMode/I");
    tree->Branch("tauPt", &tauPt,"tauPt/F");
    tree->Branch("tauEta", &tauEta,"tauEta/F");
    tree->Branch("tauMass",&tauMass,"tauMass/F");
    tree->Branch("decayMode",&decayMode,"decayMode/I");
    tree->Branch("leadPDG",&leadPDG,"leadPDG/I");
    tree->Branch("pT1",&pT1,"pT1/F");
    tree->Branch("dR1",&dR1,"dR1/F");
    tree->Branch("pT2",&pT2,"pT2/F");
    tree->Branch("dR2",&dR2,"dR2/F");
    tree->Branch("pT3",&pT3,"pT3/F");
    tree->Branch("dR3",&dR3,"dR3/F");
    tree->Branch("dxyErr1",&dxyErr1,"dxyErr1/F");
    tree->Branch("dxyErr2",&dxyErr2,"dxyErr2/F");
    tree->Branch("dxyErr3",&dxyErr3,"dxyErr3/F");
    tree->Branch("dzErr1",&dzErr1,"dzErr1/F");
    tree->Branch("dzErr2",&dzErr2,"dzErr2/F");
    tree->Branch("dzErr3",&dzErr3,"dzErr3/F");
    tree->Branch("numHits1",&numHits1,"numHits1/F");
    tree->Branch("numHits2",&numHits2,"numHits2/F");
    tree->Branch("numHits3",&numHits3,"numHits3/F");
    tree->Branch("numPixHits1",&numPixHits1,"numPixHits1/F");
    tree->Branch("numPixHits2",&numPixHits2,"numPixHits2/F");
    tree->Branch("numPixHits3",&numPixHits3,"numPixHits3/F");

}

MiniAODtwoprong::~MiniAODtwoprong()
{
}

    void

MiniAODtwoprong::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken, vertices);
    nvtx=vertices->size();
//  const reco::Vertex &PV = vertices->front();
    edm::Handle<pat::TauCollection> taus;
    iEvent.getByToken(tauToken, taus);
    edm::Handle<std::vector<reco::GenParticle> > genParticles;
    iEvent.getByToken(prunedGenToken, genParticles);
    edm::Handle<std::vector<pat::PackedCandidate>> PFCandidates;
    iEvent.getByToken(PFCandidateToken, PFCandidates);

    std::vector<const reco::GenParticle*> GenTaus;
    std::vector<const reco::GenParticle*> GenEles;
    std::vector<const reco::GenParticle*> GenMus;
    //Place generated leptons into separate lists
    for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
        if(TMath::Abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
        if(TMath::Abs(genParticle->pdgId()) == 11) GenEles.push_back(&(*genParticle));
        if(TMath::Abs(genParticle->pdgId()) == 13) GenMus.push_back(&(*genParticle));
    }
    
    //Make list of possible PF charged pions to be used in tau reconstruction
    std::vector<const pat::PackedCandidate*> hadronCandidates;
    for(std::vector<pat::PackedCandidate>::const_iterator candidate= PFCandidates->begin(); candidate != PFCandidates->end(); candidate++){
        if(TMath::Abs(candidate->pdgId()) == 211) hadronCandidates.push_back(&(*candidate));
    }
    for (size_t i=0;i<GenTaus.size();i++){
        decayMode = GetDecayMode(GenTaus[i]);
	    reco::Candidate::LorentzVector p4_vis = GetVisibleP4(GenTaus[i]);  //only look at visible decay products of gen. tau
        tauPt = p4_vis.pt();
        tauEta = p4_vis.eta();
        tauMass = p4_vis.mass();
        if (decayMode/10==3 && tauPt>18 && TMath::Abs(tauEta)<2.3){  //3-prong requirement; pT and eta cuts   
            for(const pat::Tau & tau : *taus){
                recoDecayMode = tau.decayMode();
                recoTauPt = tau.pt();
                recoTauEta = tau.eta();
                recoTauMass = tau.mass();
                int numSignalHadrons = 0;
                reco::CandidatePtr leadChargedHadr = tau.leadChargedHadrCand();
                reco::CandidatePtrVector signalCands = tau.signalCands();  //vector of the PF objects used in tau reconstruction
                for(size_t ii=0; ii<signalCands.size(); ii++){
                    if (TMath::Abs(signalCands[ii]->pdgId())==211) numSignalHadrons++;
                }
                if (recoTauPt>18 && TMath::Abs(recoTauEta)<2.3 && tau.tauID("decayModeFindingNewDMs")>0.5 && tau.tauID(tauID)>0.5 && (recoDecayMode==5 || recoDecayMode==6) && reco::deltaR(tau.eta(),tau.phi(),p4_vis.eta(),p4_vis.phi())<0.3 && numSignalHadrons == 2){  //2-prong requirement (decay mode 5,6); dR gen. tau matching; lead track tagged as charged pion  
                    int n=1;
                    const pat::PackedCandidate* firstHadron = hadronCandidates[0];  //signal Candidates do not contain hit info, so we will need to find the equivalent packed candidate
                    const pat::PackedCandidate* secondHadron = hadronCandidates[0];
                    for(size_t j=0; j<signalCands.size(); j++){
                        if (TMath::Abs(signalCands[j]->pdgId())==211){  //There should be two pions used in reconstruction of a 2-prong tau
                            if (n==1){
                                pT1 = signalCands[j]->pt();
                                dR1 = reco::deltaR(signalCands[j]->phi(),signalCands[j]->eta(),tau.phi(),tau.eta());
                                dxyErr1 = signalCands[j]->dxyError();
                                dzErr1 = signalCands[j]->dzError();
                                for(size_t l=0; l<hadronCandidates.size(); l++){
                                    if (pT1 == hadronCandidates[l]->pt()) firstHadron = hadronCandidates[l];                               
                                }
                                numHits1 = firstHadron->numberOfHits();
                                numPixHits1 = firstHadron->numberOfPixelHits();
                            }
                            if (n==2){
                                pT2 = signalCands[j]->pt();
                                dR2 = reco::deltaR(signalCands[j]->phi(),signalCands[j]->eta(),tau.phi(),tau.eta());
                                dxyErr2 = signalCands[j]->dxyError();
                                dzErr2 = signalCands[j]->dzError();
                                for(size_t l=0; l<hadronCandidates.size(); l++){
                                    if (pT1 == hadronCandidates[l]->pt()) secondHadron = hadronCandidates[l];
                                }
                                numHits2 = secondHadron->numberOfHits();
                                numPixHits2 = secondHadron->numberOfPixelHits();

                            }
                            if (n==3) std::cout << "problem";
                            n++;
                        }
                    }
                    const pat::PackedCandidate* thirdHadron = findThirdHadron(hadronCandidates,signalCands,tau);  //find nearest PF charged pion not already used in reconstruction
                    pT3 = thirdHadron->pt();
                    dR3 = reco::deltaR(thirdHadron->phi(),thirdHadron->eta(),tau.phi(),tau.eta());
                    dxyErr3 = thirdHadron->dxyError();
                    dzErr3 = thirdHadron->dzError();
                    numHits3 = thirdHadron->numberOfHits();
                    numPixHits3 = thirdHadron->numberOfPixelHits();
                    leadPDG = leadChargedHadr->pdgId();
                    myfile << tauID << " " << pT1 << " " << pT2 << " " <<pT3 << " " << dR1 << " " << dR2 << " " << dR3 << " " << dxyErr1 << " "<< dxyErr2 << " "<< dxyErr3 <<  " " << dzErr1 << " "<< dzErr2 << " "<< dzErr3 << " " << numHits1 << " " << numHits2 << " " << numHits3 << " " << numPixHits1 << " " << numPixHits2 << " " << numPixHits3 <<  "\n";
                    tree->Fill();
                }
            }
        }
    }
}
//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODtwoprong);
