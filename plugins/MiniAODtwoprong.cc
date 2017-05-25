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
#include "DataFormats/Math/interface/LorentzVector.h"
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


// class declaration
class MiniAODtwoprong : public edm::EDAnalyzer {
    public:
        explicit MiniAODtwoprong(const edm::ParameterSet&);
        ~MiniAODtwoprong();

       	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters);
        const pat::PackedCandidate* findThirdHadron(std::vector<const pat::PackedCandidate*> hadronCands, reco::CandidatePtrVector signalCands, std::vector<const reco::GenParticle*> daughters, int *recoTrack, const reco::GenParticle* *genTrack3, double *trackDR);
        int GetDecayMode(std::vector<const reco::GenParticle*>& daughters);
        reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
        void eraseHadronCands(std::vector<const pat::PackedCandidate*>& hadronCands, reco::CandidatePtrVector signalCands);
        bool isNeutrino(const reco::Candidate* daughter);

    private:
		    virtual void beginJob() override;
		    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		    virtual void endJob() override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken;
        edm::EDGetTokenT<std::vector<pat::Tau> >tauToken;
        edm::EDGetTokenT<std::vector <pat::PackedCandidate>> PFCandidateToken;
        std::string tauID;
        edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken;
        edm::EDGetTokenT<std::vector <pat::PackedGenParticle> >packedGenToken;

        TTree* tree;
        int nvtx;
        double recoTauPt;
        double recoTauEta;
        double recoTauMass;
        int recoDecayMode;
        double tauPt;
        double tauEta;
        double tauMass;
        int decayMode;
        int leadPDG;
        int hadrPDG;
        double pT1;
        double dR1;
        double pT2;
        double dR2;
        double pT3;
        double dR3;
        double dxy1;
        double dxy2;
        double dxy3;
        double dz1;
        double dz2;
        double dz3;
        double dxyErr1;
        double dxyErr2;
        double dxyErr3;
        double dzErr1;
        double dzErr2;
        double dzErr3;
        double numHits1;
        double numHits2;
        double numHits3;
        double numPixHits1;
        double numPixHits2;
        double numPixHits3;
        int   recoTrack;
        double dxy;
        double dxy_Sig;
        double dxy_error;
        double flightLengthSig;
        int hasSecondaryVertex;
        double gendxy3;
        double gendz3;
        double trackDR;
        double trackDpT;
        double gendxyErr3;
        double gendzErr3;
        const reco::GenParticle* genTrack3;

};

MiniAODtwoprong::MiniAODtwoprong(const edm::ParameterSet& iConfig):
    vtxToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    tauToken(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus"))),
    PFCandidateToken(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
    prunedGenToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken(consumes<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
    tauID = iConfig.getParameter<std::string>("tauID");
    edm::Service<TFileService> fs;

    tree = fs->make<TTree>("Ntuple", "Ntuple");
    tree->Branch("nvtx",&nvtx,"nvtx/I");
    tree->Branch("recoTauPt", &recoTauPt,"recoTauPt/D");
    tree->Branch("recoTauEta", &recoTauEta,"recoTauEta/D");
    tree->Branch("recoTauMass",&recoTauMass,"recoTauMass/D");
    tree->Branch("recoDecayMode",&recoDecayMode,"recoDecayMode/I");
    tree->Branch("tauPt", &tauPt,"tauPt/D");
    tree->Branch("tauEta", &tauEta,"tauEta/D");
    tree->Branch("tauMass",&tauMass,"tauMass/D");
    tree->Branch("decayMode",&decayMode,"decayMode/I");
    tree->Branch("leadPDG",&leadPDG,"leadPDG/I");
    tree->Branch("pT1",&pT1,"pT1/D");
    tree->Branch("dR1",&dR1,"dR1/D");
    tree->Branch("pT2",&pT2,"pT2/D");
    tree->Branch("dR2",&dR2,"dR2/D");
    tree->Branch("pT3",&pT3,"pT3/D");
    tree->Branch("dR3",&dR3,"dR3/D");
    tree->Branch("dxy1",&dxy1,"dxy1/D");
    tree->Branch("dxy2",&dxy2,"dxy2/D");
    tree->Branch("dxy3",&dxy3,"dxy3/D");
    tree->Branch("dz1",&dz1,"dz1/D");
    tree->Branch("dz2",&dz2,"dz2/D");
    tree->Branch("dz3",&dz3,"dz3/D");
    tree->Branch("dxyErr1",&dxyErr1,"dxyErr1/D");
    tree->Branch("dxyErr2",&dxyErr2,"dxyErr2/D");
    tree->Branch("dxyErr3",&dxyErr3,"dxyErr3/D");
    tree->Branch("dzErr1",&dzErr1,"dzErr1/D");
    tree->Branch("dzErr2",&dzErr2,"dzErr2/D");
    tree->Branch("dzErr3",&dzErr3,"dzErr3/D");
    tree->Branch("numHits1",&numHits1,"numHits1/D");
    tree->Branch("numHits2",&numHits2,"numHits2/D");
    tree->Branch("numHits3",&numHits3,"numHits3/D");
    tree->Branch("numPixHits1",&numPixHits1,"numPixHits1/D");
    tree->Branch("numPixHits2",&numPixHits2,"numPixHits2/D");
    tree->Branch("numPixHits3",&numPixHits3,"numPixHits3/D");
    tree->Branch("recoTrack",&recoTrack,"recoTrack/I");
    tree->Branch("dxy",&dxy,"dxy/D");
    tree->Branch("dxy_Sig",&dxy_Sig,"dxy_Sig/D");
    tree->Branch("dxy_error",&dxy_error,"dxy_error/D");
    tree->Branch("flightLengthSig",&flightLengthSig,"flightLengthSig/D");
    tree->Branch("hasSecondaryVertex",&hasSecondaryVertex,"hasSecondaryVertex/I");
    tree->Branch("gendxy3",&gendxy3,"gendxy3/D");
    tree->Branch("gendz3",&gendz3,"gendz3/D");
    tree->Branch("trackDR",&trackDR,"trackDR/D");
    tree->Branch("trackDpT",&trackDpT,"trackDpT/D");
    tree->Branch("gendxyErr3",&gendxyErr3,"gendxyErr3/D");
    tree->Branch("gendzErr3",&gendzErr3,"gendzErr3/D");
}

MiniAODtwoprong::~MiniAODtwoprong()
{
}

void MiniAODtwoprong::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken, vertices);
    nvtx=vertices->size();
    //const reco::Vertex &PV = vertices->front();
    
    edm::Handle<std::vector<pat::Tau> > taus;
    iEvent.getByToken(tauToken, taus);
    
    edm::Handle<std::vector<reco::GenParticle> > genParticles;
    iEvent.getByToken(prunedGenToken, genParticles);
    
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticles;
    iEvent.getByToken(packedGenToken, packedGenParticles);
    
    edm::Handle<std::vector<pat::PackedCandidate>> PFCandidates;
    iEvent.getByToken(PFCandidateToken, PFCandidates);

    //Pruned generator leptons
    std::vector<const reco::GenParticle*> GenTaus;
    std::vector<const reco::GenParticle*> GenEles;
    std::vector<const reco::GenParticle*> GenMus;
    //Place generated leptons into separate lists
    int genindex = 0;
    for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
      genindex++;
        if(TMath::Abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
    }
    //std::cout<<"genIndex: "<<genindex<<std::endl;
    //Packed generator hadrons 
    std::vector<const pat::PackedGenParticle*> packedHadrons;
    for(std::vector<pat::PackedGenParticle>::const_iterator genParticle = packedGenParticles->begin(); genParticle != packedGenParticles->end(); genParticle++){
        if(TMath::Abs(genParticle->pdgId()) == 211) packedHadrons.push_back(&(*genParticle));
    }

    //Make list of possible PF charged pions to be used in tau reconstruction
    std::vector<const pat::PackedCandidate*> hadronCandidates;
    for(std::vector<pat::PackedCandidate>::const_iterator candidate= PFCandidates->begin(); candidate != PFCandidates->end(); candidate++){
        if(TMath::Abs(candidate->pdgId()) == 211) hadronCandidates.push_back(&(*candidate));
    }
    //std::cout<<"GenTausSize(): "<<GenTaus.size()<<std::endl; 
    for (auto genTau : GenTaus){
        std::vector<const reco::GenParticle*> genTauDaughters;
        findDaughters(genTau, genTauDaughters);
        decayMode = GetDecayMode(genTauDaughters);
	      reco::Candidate::LorentzVector p4_vis = GetVisibleP4(genTauDaughters);  //only look at visible decay products of gen. tau
        tauPt = (float) p4_vis.pt();
        //std::cout<<"tauPt: "<<tauPt<<std::endl;
        tauEta = (float) p4_vis.eta();
        tauMass = (float) p4_vis.mass();
        if (decayMode/10==3 && tauPt>18 && TMath::Abs(tauEta)<2.3){  //3-prong requirement; pT and eta cuts   
            for(uint32_t j=0; j < taus->size(); j++){
                const pat::Tau &tau = (*taus)[j];
                recoDecayMode = tau.decayMode();
                recoTauPt = tau.pt();
                recoTauEta = tau.eta();
                recoTauMass = tau.mass();
                dxy = tau.dxy();
                dxy_error = tau.dxy_error();
                dxy_Sig = tau.dxy_Sig();
                flightLengthSig = tau.flightLengthSig();
                hasSecondaryVertex = tau.hasSecondaryVertex();
                int numSignalHadrons = 0;
                const reco::CandidatePtr leadChargedHadr = tau.leadChargedHadrCand();
                hadrPDG = leadChargedHadr->pdgId();
                std::cout<<"leadChargedHadrId: "<<hadrPDG<<std::endl;
                reco::CandidatePtrVector signalCands = tau.signalCands();  //vector of the PF objects used in tau reconstruction
                for(size_t i=0; i<signalCands.size(); i++){
                    if (TMath::Abs(signalCands[i]->pdgId())==211) numSignalHadrons++;
                }
                std::cout<<"numSignalHadrons: "<<numSignalHadrons<<std::endl;
                if (recoTauPt>18 && TMath::Abs(recoTauEta)<2.3 && tau.tauID("decayModeFindingNewDMs")>0.5 && tau.tauID(tauID)>0.5 && (recoDecayMode==5 || recoDecayMode==6) && reco::deltaR(tau.eta(),tau.phi(),p4_vis.eta(),p4_vis.phi())<0.3){  //2-prong requirement (decay mode 5,6); dR gen. tau matching; lead track tagged as charged pion  

                    int n=1;//Used to distinguish the two signal pions when filling the tree
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
                                dxy1 = firstHadron->dxy();
                                dz1 = firstHadron->dz();
                                numHits1 = firstHadron->numberOfHits();
                                numPixHits1 = firstHadron->numberOfPixelHits();
                            }
                            if (n==2){
                                pT2 = signalCands[j]->pt();
                                dR2 = reco::deltaR(signalCands[j]->phi(),signalCands[j]->eta(),tau.phi(),tau.eta());
                                dxyErr2 = signalCands[j]->dxyError();
                                dzErr2 = signalCands[j]->dzError();
                                for(size_t l=0; l<hadronCandidates.size(); l++){
                                    if (pT2 == hadronCandidates[l]->pt()) secondHadron = hadronCandidates[l];
                                }
                                dxy2 = secondHadron->dxy();
                                dz2 = secondHadron->dz();
                                numHits2 = secondHadron->numberOfHits();
                                numPixHits2 = secondHadron->numberOfPixelHits();
                            }
                            if (n==3) std::cout << "problem";
                            n++;
                        }
                    }
                    const pat::PackedCandidate* thirdHadron = findThirdHadron(hadronCandidates,signalCands,genTauDaughters, &recoTrack, &genTrack3, &trackDR);  //find nearest PF charged pion not already used in reconstruction
                    std::cout<<"thirdHadronID: "<<thirdHadron->pdgId();
                    std::cout << "\nReco 3rd pT: "<< genTrack3->pt();
                    const pat::PackedGenParticle* thirdGenHadron = packedHadrons[0];
                    float minDiffpT = 1.0; 
                    for (size_t i=0; i<packedHadrons.size(); ++i){
                        if (TMath::Abs(genTrack3->pt()-packedHadrons[i]->pt()) < minDiffpT){
                            minDiffpT = TMath::Abs(genTrack3->pt()-packedHadrons[i]->pt());
                            thirdGenHadron = packedHadrons[i];
                        } 
                    }
                    trackDpT = TMath::Abs(thirdHadron->pt() - thirdGenHadron->pt());
                    gendxy3 = thirdGenHadron->dxy();
                    gendz3 = TMath::Abs(tau.vertex().z()-genTrack3->vz());
                    gendxyErr3 = genTrack3->dxyError();
                    gendzErr3 = genTrack3->dzError();
                    std::cout <<" with gen pT of: " << thirdGenHadron->pt();
                    if (recoTrack == 1){
                        pT3 = thirdHadron->pt();
                        dR3 = reco::deltaR(thirdHadron->phi(),thirdHadron->eta(),tau.phi(),tau.eta());
                        dxy3 = thirdHadron->dxy();
                        dz3 = thirdHadron->dz();
                        dxyErr3 = thirdHadron->dxyError();
                        dzErr3 = thirdHadron->dzError();
                        numHits3 = thirdHadron->numberOfHits();
                        numPixHits3 = thirdHadron->numberOfPixelHits();
                    }
                    else{
                        pT3 = genTrack3->pt();
    
                        dR3 = reco::deltaR(genTrack3->phi(),genTrack3->eta(),tau.phi(),tau.eta());
                        dxy3 = -1;//genTrack3->dxy();
                        dz3 = -1;//;genTrack3->dz();
                        //pat::genParticle genTrack3= pat::genParticle(genTrack3);
                        dxyErr3 = genTrack3->dxyError();
                        dzErr3 = genTrack3->dzError();
                        numHits3 = -1;
                        numPixHits3 = -1;
                    }
                }

                    leadPDG = leadChargedHadr->pdgId();
            }
        }
    }

                    tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAODtwoprong::beginJob()
{
}


// ------------ method called once each job just after ending the event loop  ------------
	void 
MiniAODtwoprong::endJob() 
{
}

void MiniAODtwoprong::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
{
    unsigned numDaughters = mother->numberOfDaughters();
    if (numDaughters == 0) std::cout << " none ";
    for (unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
        const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();
        if (daughter->status() == 1){  //status = 1 is a final state daughter
            daughters.push_back(daughter); 
        }
        if (daughter->status() == 2){  //status = 2 is an intermediate daughter; will decay further
            daughters.push_back(daughter); 
            findDaughters(daughter, daughters);
        }
    }
}

//Returns a decay mode for a generator tau (input its vector of daughters)
int MiniAODtwoprong::GetDecayMode(std::vector<const reco::GenParticle*>& daughters){
	int decayMode = -1;
	std::vector<int> counts(4,0);
    for(size_t i=0; i<daughters.size(); ++i){
		int pdg = daughters[i]->pdgId();
		if (TMath::Abs(pdg)==11) ++counts[0];  //electrons
		if (TMath::Abs(pdg)==13) ++counts[1];  //muons
		if (TMath::Abs(pdg)==111) ++counts[2];  //neutral pions
        if (TMath::Abs(pdg)==211) ++counts[3];  //charged pions
 	}
	if (counts[0] > 0) decayMode = 3;
	if (counts[1] > 0) decayMode = 4;
	if (counts[3] > 0) decayMode = 10*counts[3]+counts[2];  //First digit is # of prongs (1 or 3), second is # of neutral pions
	return decayMode;
}

bool MiniAODtwoprong::isNeutrino(const reco::Candidate* daughter)
{
  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}
//Gets visible 4-momentum of a particle from list of daughters
reco::Candidate::LorentzVector MiniAODtwoprong::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
	reco::Candidate::LorentzVector p4_vis(0,0,0,0);	
 	for(size_t i = 0; i < daughters.size(); ++i){
 		if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
 			p4_vis += daughters[i]->p4();
		}
 	}
	return p4_vis;
}

void MiniAODtwoprong::eraseHadronCands(std::vector<const pat::PackedCandidate*>& hadronCands, reco::CandidatePtrVector signalCands){
    for(size_t i=0; i<signalCands.size(); ++i){
        for(size_t j=0; j<hadronCands.size(); ++j){
            if(signalCands[i]->pt()==hadronCands[j]->pt()) hadronCands.erase(hadronCands.begin()+j);
        }
    }
}
//Finds closest charged pion not already used in 2-prong tau reconstruction
const pat::PackedCandidate* MiniAODtwoprong::findThirdHadron(std::vector<const pat::PackedCandidate*> hadronCands, reco::CandidatePtrVector signalCands, std::vector<const reco::GenParticle*> daughters,int *recoTrack, const reco::GenParticle* *genTrack3, double *trackDR){
    const pat::PackedCandidate* thirdHadron = hadronCands[0];
    std::vector<const reco::GenParticle*> genTracks;  //make a list of charged pion daughters
    for(size_t i=0; i<daughters.size(); ++i){
        if (TMath::Abs(daughters[i]->pdgId())==211) genTracks.push_back(daughters[i]);
    }
    for(size_t i=0; i<signalCands.size(); ++i){  //loop through the generator charged pions
        if (TMath::Abs(signalCands[i]->pdgId())==211){
            float minDR = 1.0;
            int index = -1;
            for(size_t j=0; j<genTracks.size(); ++j){  //look for a signal candidate match
                float dR = reco::deltaR(genTracks[j]->eta(),genTracks[j]->phi(),signalCands[i]->eta(),signalCands[i]->phi());
                if (dR < minDR){
                    minDR = dR;
                    index = j;
                }
            }
            genTracks.erase(genTracks.begin()+index);  //remove the generator charged pion
        }
    }
    *genTrack3 = genTracks[0];
    eraseHadronCands(hadronCands,signalCands);
    float minDR = 1.0; //reusing minDR to match third track with hadron candidate
    for(size_t i=0; i<hadronCands.size(); ++i){
        if (hadronCands[i]->pdgId()==genTracks[0]->pdgId()){
            float dR = reco::deltaR(genTracks[0]->eta(),genTracks[0]->phi(),hadronCands[i]->eta(),hadronCands[i]->phi());
            if (dR < minDR && TMath::Abs(hadronCands[i]->pt()-genTracks[0]->pt()) < 5){
                minDR = dR;
                thirdHadron = hadronCands[i];
            }   
        }
    }
    *trackDR = minDR;
    if (minDR > 0.02) *recoTrack = 0;
    else *recoTrack = 1;
    std::cout <<"\n\nThird reco track by dR: dPT = " << TMath::Abs(thirdHadron->pt()-genTracks[0]->pt()) << " dR = " << reco::deltaR(thirdHadron->eta(),thirdHadron->phi(),genTracks[0]->eta(),genTracks[0]->phi());     
    std::cout <<"\nthird pt: " <<thirdHadron->pt() << " gen pt: " <<genTracks[0]->pt();
    
    return thirdHadron;
}
void
MiniAODtwoprong::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	  edm::ParameterSetDescription desc;
	  desc.setUnknown();
	  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODtwoprong);
