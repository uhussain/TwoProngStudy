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
#include <TROOT.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include <stdio.h>
#include <fstream>
#include "util.h"

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
        edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken;
        edm::EDGetTokenT<std::vector<pat::Tau> >tauToken;
        edm::EDGetTokenT<std::vector <pat::PackedCandidate>> PFCandidateToken;
        std::string tauID;
        edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken;
        edm::EDGetTokenT<std::vector <pat::PackedGenParticle> >packedGenToken;

        TTree* tree;
        int nvtx;
        int nChHadrj1;
        double recoTrkj1;
        double Trk1Pt;
        double Trk1PtFrac;
        double Trk1Eta;
        double Trk1Phi;
        double ChHadr1_pt;
        double ChHadr1_PtDiff;
        double ChHadr1_ptfrac;
        double ChHadr1_eta;
        double ChHadr1_phi;
        double ChHadr1_phiAtVtx;//this is identical to phi() for the vast majority
        //of the particles, but the two might differ for some of them if the calorimeters had contributed significantly in defining the
        //4-vector of the particle.
        float ChHadr1_dxy;//longitudinal and transverse impact parameters with respect to the PV: dxy(), dzAssociatedPVAssociatedPV().
        float ChHadr1_dzAssociatedPV;
        double ChHadr1_vx;//returns the position of the point of closest approach to the PV 
        double ChHadr1_vy;
        double ChHadr1_vz;
        double ChHadr1_vertexRef;//reference to the PV itself
        double ChHadr1_fromPV;//returns a number between 3 and 0 to define how tight the association with the first PV is
        double ChHadr1_numHits;
        double ChHadr1_numPixelHits;


        double Trk2Pt;
        double Trk2PtFrac;
        double Trk2Eta;
        double Trk2Phi;
        
        double dRTrk12; //deltaR between two highest PtTracks
        
        double ChHadr2_pt;
        double ChHadr2_PtDiff;
        double ChHadr2_ptfrac;
        double ChHadr2_eta;
        double ChHadr2_phi;
        double ChHadr2_phiAtVtx;//this is identical to phi() for the vast majority
        //of the particles, but the two might differ for some of them if the calorimeters had contributed significantly in defining the
        //4-vector of the particle.
        float ChHadr2_dxy;//longitudinal and transverse impact parameters with respect to the PV: dxy(), dzAssociatedPVAssociatedPV().
        float ChHadr2_dzAssociatedPV;
        double ChHadr2_vx;//returns the x-coordinate of vertex position
        double ChHadr2_vy;
        double ChHadr2_vz;
        double ChHadr2_vertexRef;//reference to the PV itself
        double ChHadr2_fromPV;//returns a number between 3 and 0 to define how tight the association with the first PV is
        double ChHadr2_numHits;
        double ChHadr2_numPixelHits;
};

MiniAODtwoprong::MiniAODtwoprong(const edm::ParameterSet& iConfig):
    vtxToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    jetsToken(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    tauToken(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus"))),
    PFCandidateToken(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
    prunedGenToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken(consumes<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
    tauID = iConfig.getParameter<std::string>("tauID");
    edm::Service<TFileService> fs;

    tree = fs->make<TTree>("Ntuple", "Ntuple");
    tree->Branch("nvtx",&nvtx,"nvtx/I");
    tree->Branch("nChHadrj1",&nChHadrj1,"nChHadrj1");
    tree->Branch("recoTrkj1",&recoTrkj1,"recoTrkj1");
    
    tree->Branch("Trk1Pt",&Trk1Pt,"Trk1Pt/D");
    tree->Branch("Trk1Eta",&Trk1Eta,"Trk1Eta/D");
    tree->Branch("Trk1Phi",&Trk1Phi,"Trk1Phi/D");
    tree->Branch("Trk1PtFrac",&Trk1PtFrac,"Trk1PtFrac/D");
    
    tree->Branch("ChHadr1_pt",&ChHadr1_pt,"ChHadr1_pt/D");
    tree->Branch("ChHadr1_PtDiff",&ChHadr1_PtDiff,"ChHadr1_PtDiff/D");
    tree->Branch("ChHadr1_ptfrac",&ChHadr1_ptfrac,"ChHadr1_ptfrac/D");
    tree->Branch("ChHadr1_eta",&ChHadr1_eta,"ChHadr1_eta/D");
    tree->Branch("ChHadr1_phi",&ChHadr1_phi,"ChHadr1_phi/D");
    tree->Branch("ChHadr1_phiAtVtx",&ChHadr1_phiAtVtx,"ChHadr1_phiAtVtx/D");
    tree->Branch("ChHadr1_dxy",&ChHadr1_dxy,"ChHadr1_dxy/F");
    tree->Branch("ChHadr1_dzAssociatedPV",&ChHadr1_dzAssociatedPV,"ChHadr1_dzAssociatedPV/F");
    tree->Branch("ChHadr1_vx",&ChHadr1_vx,"ChHadr1_vx/D");
    tree->Branch("ChHadr1_vy",&ChHadr1_vy,"ChHadr1_vy/D");
    tree->Branch("ChHadr1_vz",&ChHadr1_vz,"ChHadr1_vz/D");
    tree->Branch("ChHadr1_fromPV",&ChHadr1_fromPV,"ChHadr1_fromPV/D");
    tree->Branch("ChHadr1_numHits",&ChHadr1_numHits,"ChHadr1_numHits/D");
    tree->Branch("ChHadr1_numPixelHits",&ChHadr1_numPixelHits,"ChHadr1_numPixelHits/D");

    tree->Branch("Trk2Pt",&Trk2Pt,"Trk2Pt/D");
    tree->Branch("Trk2Eta",&Trk2Eta,"Trk2Eta/D");
    tree->Branch("Trk2Phi",&Trk2Phi,"Trk2Phi/D");
    tree->Branch("Trk2PtFrac",&Trk2PtFrac,"Trk2PtFrac/D");

    tree->Branch("dRTrk12",&dRTrk12,"dRTrk12/D");
    
    tree->Branch("ChHadr2_pt",&ChHadr2_pt,"ChHadr2_pt/D");
    tree->Branch("ChHadr2_PtDiff",&ChHadr2_PtDiff,"ChHadr2_PtDiff/D");
    tree->Branch("ChHadr2_ptfrac",&ChHadr2_ptfrac,"ChHadr2_ptfrac/D");
    tree->Branch("ChHadr2_eta",&ChHadr2_eta,"ChHadr2_eta/D");
    tree->Branch("ChHadr2_phi",&ChHadr2_phi,"ChHadr2_phi/D");
    tree->Branch("ChHadr2_phiAtVtx",&ChHadr2_phiAtVtx,"ChHadr2_phiAtVtx/D");
    tree->Branch("ChHadr2_dxy",&ChHadr2_dxy,"ChHadr2_dxy/F");
    tree->Branch("ChHadr2_dzAssociatedPV",&ChHadr2_dzAssociatedPV,"ChHadr2_dzAssociatedPV/F");
    tree->Branch("ChHadr2_vx",&ChHadr2_vx,"ChHadr2_vx/D");
    tree->Branch("ChHadr2_vy",&ChHadr2_vy,"ChHadr2_vy/D");
    tree->Branch("ChHadr2_vz",&ChHadr2_vz,"ChHadr2_vz/D");
    tree->Branch("ChHadr2_fromPV",&ChHadr2_fromPV,"ChHadr2_fromPV/D");
    tree->Branch("ChHadr2_numHits",&ChHadr2_numHits,"ChHadr2_numHits/D");
    tree->Branch("ChHadr2_numPixelHits",&ChHadr2_numPixelHits,"ChHadr2_numPixelHits/D");
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
    
    edm::Handle<std::vector<pat::Jet> > ak4jets;
    iEvent.getByToken(jetsToken, ak4jets);

    edm::Handle<std::vector<pat::Tau> > taus;
    iEvent.getByToken(tauToken, taus);
    
    edm::Handle<std::vector<reco::GenParticle> > genParticles;
    iEvent.getByToken(prunedGenToken, genParticles);
    
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticles;
    iEvent.getByToken(packedGenToken, packedGenParticles);
    
    edm::Handle<std::vector<pat::PackedCandidate>> PFCandidates;
    iEvent.getByToken(PFCandidateToken, PFCandidates);

    std::vector<const pat::PackedGenParticle*> packedHadrons;
    for(std::vector<pat::PackedGenParticle>::const_iterator genParticle = packedGenParticles->begin(); genParticle != packedGenParticles->end(); genParticle++){
        if(TMath::Abs(genParticle->pdgId()) == 211) packedHadrons.push_back(&(*genParticle));
    }

    //Make list of possible PF charged pions to be used in tau reconstruction
    std::vector<const pat::PackedCandidate*> hadronCandidates;
    for(std::vector<pat::PackedCandidate>::const_iterator candidate= PFCandidates->begin(); candidate != PFCandidates->end(); candidate++){
        if(TMath::Abs(candidate->pdgId()) == 211) hadronCandidates.push_back(&(*candidate));
    }
    for(uint32_t i=0; i < ak4jets->size(); i++){
        recoTrkj1 = -99; //count the number of reco Tracks within dR of leading Jet
        const pat::Jet &jet = (*ak4jets)[i];
        if (i==0){   
        double difference = 0;
        //These vectors will store the PtDifference between a packed PFCandidate and a nearbyJet
        std::vector<std::pair<double,const pat::PackedCandidate*>> PtDiffChHadr;//pdgId = abs(211) 
        std::vector<std::pair<double,const pat::PackedCandidate*>> PtDiffNeutrHadr;//pdgId = abs(130)
        std::vector<std::pair<double,const pat::PackedCandidate*>> PtDiffPhotons;//pdgId = abs(22)
        //Loop over all PFCandidates
        for(uint32_t j = 0; j < PFCandidates->size(); j++) {
            const pat::PackedCandidate &pfCand = (*PFCandidates)[j]; 
            if(reco::deltaR(jet.eta(),jet.phi(),pfCand.eta(),pfCand.phi())< 0.1){
                if(abs(pfCand.pdgId())== 211) {
                  difference = abs(jet.pt()-pfCand.pt());
                  PtDiffChHadr.push_back({difference,&pfCand});}

                if(abs(pfCand.pdgId())== 130) {
                  difference = abs(jet.pt()-pfCand.pt());
                  PtDiffNeutrHadr.push_back(std::make_pair(difference,&pfCand));}

                if(abs(pfCand.pdgId())== 22) {
                  difference = abs(jet.pt()-pfCand.pt());
                  PtDiffPhotons.push_back(std::make_pair(difference,&pfCand));}
      }
    }
        //Sort these PtDifference vectors by ascending order in PtDifference from the Jet
        std::sort(PtDiffChHadr.begin(),PtDiffChHadr.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
        //std::sort(PtDiffChHadr.begin(),PtDiffChHadr.end(),pairCompare);
        std::sort(PtDiffNeutrHadr.begin(),PtDiffNeutrHadr.end(),pairCompare); //or you can use the specific pairCompare function. Sort has some issues with template that I cannot fix right now?
        std::sort(PtDiffPhotons.begin(),PtDiffPhotons.end(),pairCompare);//pairCompare is defined in util.h

        if(PtDiffChHadr.size()>0){
          nChHadrj1 = PtDiffChHadr.size();
          for (int i=0;i<nChHadrj1;i++){
            const pat::PackedCandidate &ChHadr = *(PtDiffChHadr.at(i).second);
            if(ChHadr.bestTrack()!=nullptr){//method points to null if there is no track
              recoTrkj1++;
            }

            const pat::PackedCandidate &ChHadr1 = *(PtDiffChHadr.at(0).second);
            if(ChHadr1.bestTrack()!=nullptr){//lets store some info about Trk1 associated with highestPt Hadron
              const reco::Track &Trk = *(ChHadr1.bestTrack());
              Trk1Pt = Trk.pt();
              Trk1PtFrac = (Trk1Pt/jet.pt());
              Trk1Eta = Trk.eta();
              Trk1Phi = Trk.phi();
              }
            ChHadr1_pt = ChHadr1.pt();
            ChHadr1_PtDiff = PtDiffChHadr.at(0).first;
            ChHadr1_ptfrac = (ChHadr1_pt/jet.pt());
            ChHadr1_eta = ChHadr1.eta();
            ChHadr1_phi = ChHadr1.phi();
            ChHadr1_phiAtVtx = ChHadr1.phiAtVtx();
            ChHadr1_dxy = ChHadr1.dxy();
            ChHadr1_dzAssociatedPV = ChHadr1.dzAssociatedPV();
            ChHadr1_vx = ChHadr1.vx();
            ChHadr1_vy = ChHadr1.vy();
            ChHadr1_vz = ChHadr1.vz();
            ChHadr1_fromPV = ChHadr1.fromPV();
            ChHadr1_numHits = ChHadr1.numberOfHits();
            ChHadr1_numPixelHits = ChHadr1.numberOfPixelHits();

            const pat::PackedCandidate &ChHadr2 = *(PtDiffChHadr.at(1).second);
            if(ChHadr2.bestTrack()!=nullptr){//lets store some info about Trk2 associated with second highestPt charged Hadron
              const reco::Track &Trk = *(ChHadr2.bestTrack());
              Trk2Pt = Trk.pt();
              Trk2PtFrac = (Trk2Pt/jet.pt());
              Trk2Eta = Trk.eta();
              Trk2Phi = Trk.phi();
              }
            ChHadr2_pt = ChHadr2.pt();
            ChHadr2_PtDiff = PtDiffChHadr.at(1).first;
            ChHadr2_ptfrac = (ChHadr2_pt/jet.pt());
            ChHadr2_eta = ChHadr2.eta();
            ChHadr2_phi = ChHadr2.phi();
            ChHadr2_phiAtVtx = ChHadr2.phiAtVtx();
            ChHadr2_dxy = ChHadr2.dxy();
            ChHadr2_dzAssociatedPV = ChHadr2.dzAssociatedPV();
            ChHadr2_vx = ChHadr2.vx();
            ChHadr2_vy = ChHadr2.vy();
            ChHadr2_vz = ChHadr2.vz();
            ChHadr2_fromPV = ChHadr2.fromPV();
            ChHadr2_numHits = ChHadr2.numberOfHits();
            ChHadr2_numPixelHits = ChHadr2.numberOfPixelHits();
          }
            dRTrk12 = reco::deltaR(Trk1Eta,Trk1Phi,Trk2Eta,Trk2Phi);
          }

  } 
}

 //   for(uint32_t j=0; j < taus->size(); j++){
 //       const pat::Tau &tau = (*taus)[j];
 //       recoDecayMode = tau.decayMode();
 //       recoTauPt = tau.pt();
 //       recoTauEta = tau.eta();
 //       recoTauMass = tau.mass();
 //       dxy = tau.dxy();
 //       dxy_error = tau.dxy_error();
 //       dxy_Sig = tau.dxy_Sig();
 //       flightLengthSig = tau.flightLengthSig();
 //       hasSecondaryVertex = tau.hasSecondaryVertex();
 //       int numSignalHadrons = 0;
 //       const reco::CandidatePtr leadChargedHadr = tau.leadChargedHadrCand();
 //       hadrPDG = leadChargedHadr->pdgId();
 //       std::cout<<"leadChargedHadrId: "<<hadrPDG<<std::endl;
 //       reco::CandidatePtrVector signalCands = tau.signalCands();  //vector of the PF objects used in tau reconstruction
 //       for(size_t i=0; i<signalCands.size(); i++){
 //           if (TMath::Abs(signalCands[i]->pdgId())==211) numSignalHadrons++;
 //       }
 //       std::cout<<"numSignalHadrons: "<<numSignalHadrons<<std::endl;
 //       if (recoTauPt>18 && TMath::Abs(recoTauEta)<2.3 && tau.tauID("decayModeFindingNewDMs")>0.5 && tau.tauID(tauID)>0.5 && (recoDecayMode==5 || recoDecayMode==6) && reco::deltaR(tau.eta(),tau.phi(),p4_vis.eta(),p4_vis.phi())<0.3){  //2-prong requirement (decay mode 5,6); dR gen. tau matching; lead track tagged as charged pion  

 //           int n=1;//Used to distinguish the two signal pions when filling the tree
 //           const pat::PackedCandidate* firstHadron = hadronCandidates[0];  //signal Candidates do not contain hit info, so we will need to find the equivalent packed candidate
 //           const pat::PackedCandidate* secondHadron = hadronCandidates[0];
 //           for(size_t j=0; j<signalCands.size(); j++){
 //               if (TMath::Abs(signalCands[j]->pdgId())==211){  //There should be two pions used in reconstruction of a 2-prong tau
 //                   if (n==1){
 //                       pT1 = signalCands[j]->pt();
 //                       dR1 = reco::deltaR(signalCands[j]->phi(),signalCands[j]->eta(),tau.phi(),tau.eta());
 //                       dxyErr1 = signalCands[j]->dxyError();
 //                       dzErr1 = signalCands[j]->dzError();
 //                       for(size_t l=0; l<hadronCandidates.size(); l++){
 //                           if (pT1 == hadronCandidates[l]->pt()) firstHadron = hadronCandidates[l];                               
 //                       }
 //                       dxy1 = firstHadron->dxy();
 //                       dz1 = firstHadron->dz();
 //                       numHits1 = firstHadron->numberOfHits();
 //                       numPixHits1 = firstHadron->numberOfPixelHits();
 //                   }
 //                   if (n==2){
 //                       pT2 = signalCands[j]->pt();
 //                       dR2 = reco::deltaR(signalCands[j]->phi(),signalCands[j]->eta(),tau.phi(),tau.eta());
 //                       dxyErr2 = signalCands[j]->dxyError();
 //                       dzErr2 = signalCands[j]->dzError();
 //                       for(size_t l=0; l<hadronCandidates.size(); l++){
 //                           if (pT2 == hadronCandidates[l]->pt()) secondHadron = hadronCandidates[l];
 //                       }
 //                       dxy2 = secondHadron->dxy();
 //                       dz2 = secondHadron->dz();
 //                       numHits2 = secondHadron->numberOfHits();
 //                       numPixHits2 = secondHadron->numberOfPixelHits();
 //                   }
 //                   if (n==3) std::cout << "problem";
 //                   n++;
 //               }
 //           }
 //           const pat::PackedCandidate* thirdHadron = findThirdHadron(hadronCandidates,signalCands,genTauDaughters, &recoTrack, &genTrack3, &trackDR);  //find nearest PF charged pion not already used in reconstruction
 //           std::cout<<"thirdHadronID: "<<thirdHadron->pdgId();
 //           std::cout << "\nReco 3rd pT: "<< genTrack3->pt();
 //           const pat::PackedGenParticle* thirdGenHadron = packedHadrons[0];
 //           float minDiffpT = 1.0; 
 //           for (size_t i=0; i<packedHadrons.size(); ++i){
 //               if (TMath::Abs(genTrack3->pt()-packedHadrons[i]->pt()) < minDiffpT){
 //                   minDiffpT = TMath::Abs(genTrack3->pt()-packedHadrons[i]->pt());
 //                   thirdGenHadron = packedHadrons[i];
 //               } 
 //           }
 //           trackDpT = TMath::Abs(thirdHadron->pt() - thirdGenHadron->pt());
 //           gendxy3 = thirdGenHadron->dxy();
 //           gendz3 = TMath::Abs(tau.vertex().z()-genTrack3->vz());
 //           gendxyErr3 = genTrack3->dxyError();
 //           gendzErr3 = genTrack3->dzError();
 //           std::cout <<" with gen pT of: " << thirdGenHadron->pt();
 //           if (recoTrack == 1){
 //               pT3 = thirdHadron->pt();
 //               dR3 = reco::deltaR(thirdHadron->phi(),thirdHadron->eta(),tau.phi(),tau.eta());
 //               dxy3 = thirdHadron->dxy();
 //               dz3 = thirdHadron->dz();
 //               dxyErr3 = thirdHadron->dxyError();
 //               dzErr3 = thirdHadron->dzError();
 //               numHits3 = thirdHadron->numberOfHits();
 //               numPixHits3 = thirdHadron->numberOfPixelHits();
 //           }
 //           else{
 //               pT3 = genTrack3->pt();
 //   
 //               dR3 = reco::deltaR(genTrack3->phi(),genTrack3->eta(),tau.phi(),tau.eta());
 //               dxy3 = -1;//genTrack3->dxy();
 //               dz3 = -1;//;genTrack3->dz();
 //               //pat::genParticle genTrack3= pat::genParticle(genTrack3);
 //               dxyErr3 = genTrack3->dxyError();
 //               dzErr3 = genTrack3->dzError();
 //               numHits3 = -1;
 //               numPixHits3 = -1;
 //           }
 //       }

 //           leadPDG = leadChargedHadr->pdgId();
 //           }
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

//void MiniAODtwoprong::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
//{
//    unsigned numDaughters = mother->numberOfDaughters();
//    if (numDaughters == 0) std::cout << " none ";
//    for (unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
//        const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();
//        if (daughter->status() == 1){  //status = 1 is a final state daughter
//            daughters.push_back(daughter); 
//        }
//        if (daughter->status() == 2){  //status = 2 is an intermediate daughter; will decay further
//            daughters.push_back(daughter); 
//            findDaughters(daughter, daughters);
//        }
//    }
//}
//
////Returns a decay mode for a generator tau (input its vector of daughters)
//int MiniAODtwoprong::GetDecayMode(std::vector<const reco::GenParticle*>& daughters){
//	int decayMode = -1;
//	std::vector<int> counts(4,0);
//    for(size_t i=0; i<daughters.size(); ++i){
//		int pdg = daughters[i]->pdgId();
//		if (TMath::Abs(pdg)==11) ++counts[0];  //electrons
//		if (TMath::Abs(pdg)==13) ++counts[1];  //muons
//		if (TMath::Abs(pdg)==111) ++counts[2];  //neutral pions
//        if (TMath::Abs(pdg)==211) ++counts[3];  //charged pions
// 	}
//	if (counts[0] > 0) decayMode = 3;
//	if (counts[1] > 0) decayMode = 4;
//	if (counts[3] > 0) decayMode = 10*counts[3]+counts[2];  //First digit is # of prongs (1 or 3), second is # of neutral pions
//	return decayMode;
//}
//
//bool MiniAODtwoprong::isNeutrino(const reco::Candidate* daughter)
//{
//  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
//}
////Gets visible 4-momentum of a particle from list of daughters
//reco::Candidate::LorentzVector MiniAODtwoprong::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
//	reco::Candidate::LorentzVector p4_vis(0,0,0,0);	
// 	for(size_t i = 0; i < daughters.size(); ++i){
// 		if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
// 			p4_vis += daughters[i]->p4();
//		}
// 	}
//	return p4_vis;
//}
//
//void MiniAODtwoprong::eraseHadronCands(std::vector<const pat::PackedCandidate*>& hadronCands, reco::CandidatePtrVector signalCands){
//    for(size_t i=0; i<signalCands.size(); ++i){
//        for(size_t j=0; j<hadronCands.size(); ++j){
//            if(signalCands[i]->pt()==hadronCands[j]->pt()) hadronCands.erase(hadronCands.begin()+j);
//        }
//    }
//}
////Finds closest charged pion not already used in 2-prong tau reconstruction
//const pat::PackedCandidate* MiniAODtwoprong::findThirdHadron(std::vector<const pat::PackedCandidate*> hadronCands, reco::CandidatePtrVector signalCands, std::vector<const reco::GenParticle*> daughters,int *recoTrack, const reco::GenParticle* *genTrack3, double *trackDR){
//    const pat::PackedCandidate* thirdHadron = hadronCands[0];
//    std::vector<const reco::GenParticle*> genTracks;  //make a list of charged pion daughters
//    for(size_t i=0; i<daughters.size(); ++i){
//        if (TMath::Abs(daughters[i]->pdgId())==211) genTracks.push_back(daughters[i]);
//    }
//    for(size_t i=0; i<signalCands.size(); ++i){  //loop through the generator charged pions
//        if (TMath::Abs(signalCands[i]->pdgId())==211){
//            float minDR = 1.0;
//            int index = -1;
//            for(size_t j=0; j<genTracks.size(); ++j){  //look for a signal candidate match
//                float dR = reco::deltaR(genTracks[j]->eta(),genTracks[j]->phi(),signalCands[i]->eta(),signalCands[i]->phi());
//                if (dR < minDR){
//                    minDR = dR;
//                    index = j;
//                }
//            }
//            genTracks.erase(genTracks.begin()+index);  //remove the generator charged pion
//        }
//    }
//    *genTrack3 = genTracks[0];
//    eraseHadronCands(hadronCands,signalCands);
//    float minDR = 1.0; //reusing minDR to match third track with hadron candidate
//    for(size_t i=0; i<hadronCands.size(); ++i){
//        if (hadronCands[i]->pdgId()==genTracks[0]->pdgId()){
//            float dR = reco::deltaR(genTracks[0]->eta(),genTracks[0]->phi(),hadronCands[i]->eta(),hadronCands[i]->phi());
//            if (dR < minDR && TMath::Abs(hadronCands[i]->pt()-genTracks[0]->pt()) < 5){
//                minDR = dR;
//                thirdHadron = hadronCands[i];
//            }   
//        }
//    }
//    *trackDR = minDR;
//    if (minDR > 0.02) *recoTrack = 0;
//    else *recoTrack = 1;
//    std::cout <<"\n\nThird reco track by dR: dPT = " << TMath::Abs(thirdHadron->pt()-genTracks[0]->pt()) << " dR = " << reco::deltaR(thirdHadron->eta(),thirdHadron->phi(),genTracks[0]->eta(),genTracks[0]->phi());     
//    std::cout <<"\nthird pt: " <<thirdHadron->pt() << " gen pt: " <<genTracks[0]->pt();
//    
//    return thirdHadron;
//}
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
