// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
//#include "helpers.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <TROOT.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include <stdio.h>
#include <fstream>

// class declaration
class AODtwoprong : public edm::EDAnalyzer {
    public:
        explicit AODtwoprong(const edm::ParameterSet&);
        ~AODtwoprong();

       	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters);
        const reco::PFCandidate* findThirdHadron(std::vector<const reco::PFCandidate*> hadronCands, reco::CandidatePtrVector signalCands, std::vector<const reco::GenParticle*> daughters, int *recoTrack, const reco::GenParticle* *genTrack3, double *trackDR);
        int GetDecayMode(std::vector<const reco::GenParticle*>& daughters);
        reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
        void eraseHadronCands(std::vector<const reco::PFCandidate*>& hadronCands, reco::CandidatePtrVector signalCands);
        bool isNeutrino(const reco::Candidate* daughter);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);           
           
    private:
		    virtual void beginJob() override;
		    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		    virtual void endJob() override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken;
        edm::EDGetTokenT<std::vector<reco::PFJet> > jetsToken;
        edm::EDGetTokenT<std::vector <reco::PFCandidate>> PFCandidateToken;
        edm::EDGetTokenT<std::vector<reco::Track>> TrackToken;
        edm::EDGetTokenT<std::vector <reco::GenParticle>> GenToken;
        edm::EDGetTokenT<std::vector<reco::PFMET> > pfMETsToken;
        edm::EDGetTokenT<std::vector<reco::CaloMET> > caloMETToken;

        TTree* RecoTree;
        int     run_;
        int  event_;
        int     lumis_;
        int nvtx;
        
        float pfMET;
        float caloMET;
        double j1Pt;
        double j1Eta;
        double j1Phi;
        double j1CHF;
        double j1NHF;
        float j1ConsEtaPhiSpread;
       
        //recoTrk variables
        int nPosTrkj1;
        int nNegTrkj1;
        int nTotTrkj1;
        double posTrk1Pt;
        double negTrk1Pt;
        double posTrk1Charge;
        double negTrk1Charge;

        double posTrk2Pt;
        double negTrk2Pt;
        double posTrk2Charge;
        double negTrk2Charge;

        double posTrk3Pt;
        double negTrk3Pt;
        double posTrk3Charge;
        double negTrk3Charge;

        double dRPosTrk1negTrk1;
        
        //PosCharged Hadron Variables
        int nChHadrPosj1;
        int ChHadrPosrecoTrkj1;

        double ChHadrPosTrk1Pt;
        double ChHadrPosTrk1_Charge;
        double ChHadrPosTrk1PtFrac;
        double ChHadrPosTrk1Eta;
        double ChHadrPosTrk1Phi;
        double ChHadrPos1_pt;
        double ChHadrPos1_Charge;
        double ChHadrPos1_PtDiff;
        double ChHadrPos1_ptfrac;
        double ChHadrPos1_eta;
        double ChHadrPos1_phi;
        double ChHadrPos1_vx;//returns the position of the point of closest approach to the PV 
        double ChHadrPos1_vy;
        double ChHadrPos1_vz;


        double ChHadrPosTrk2Pt; 
        double ChHadrPosTrk2_Charge;
        double ChHadrPosTrk2PtFrac;
        double ChHadrPosTrk2Eta;
        double ChHadrPosTrk2Phi;
        
        double dRChHadrPosTrk12; //deltaR between two highest PtTracks
        
        double ChHadrPos2_pt;
        double ChHadrPos2_Charge;
        double ChHadrPos2_PtDiff;
        double ChHadrPos2_ptfrac;
        double ChHadrPos2_eta;
        double ChHadrPos2_phi;
        double ChHadrPos2_vx;//returns the x-coordinate of vertex position
        double ChHadrPos2_vy;
        double ChHadrPos2_vz;

        double dRChHadrPos12; //deltaR between two highest PtTracks
        double ChHadrPos12_Ptfrac;

        double ChHadrPosTrk12_Ptfrac;

        //Neg ChHadr Hadron Variables
        int nChHadrNegj1;

        int ChHadrNegrecoTrkj1;

        double ChHadrNegTrk1Pt;
        double ChHadrNegTrk1_Charge;
        double ChHadrNegTrk1PtFrac;
        double ChHadrNegTrk1Eta;
        double ChHadrNegTrk1Phi;
        double ChHadrNeg1_pt;
        double ChHadrNeg1_Charge;
        double ChHadrNeg1_PtDiff;
        double ChHadrNeg1_ptfrac;
        double ChHadrNeg1_eta;
        double ChHadrNeg1_phi;
        double ChHadrNeg1_vx;//returns the position of the point of closest approach to the PV 
        double ChHadrNeg1_vy;
        double ChHadrNeg1_vz;


        double ChHadrNegTrk2Pt; 
        double ChHadrNegTrk2_Charge;
        double ChHadrNegTrk2PtFrac;
        double ChHadrNegTrk2Eta;
        double ChHadrNegTrk2Phi;
        
        double dRChHadrNegTrk12; //deltaR between two highest PtTracks
        
        double ChHadrNeg2_pt;
        double ChHadrNeg2_Charge;
        double ChHadrNeg2_PtDiff;
        double ChHadrNeg2_ptfrac;
        double ChHadrNeg2_eta;
        double ChHadrNeg2_phi;
        double ChHadrNeg2_vx;//returns the x-coordinate of vertex position
        double ChHadrNeg2_vy;
        double ChHadrNeg2_vz;

        double dRChHadrNeg12; //deltaR between two highest PtTracks
        double ChHadrNeg12_Ptfrac;

        double ChHadrNegTrk12_Ptfrac;

        double dRPosHadr1NegHadr1;

        //Photons in GenParticles
        int nPhoj1;
        double PhoTotalCharge;

        double Pho1_pt;
        double Pho1_Charge;
        double Pho1_PtDiff;
        double Pho1_ptfrac;
        double Pho1_eta;
        double Pho1_phi;
        double Pho1_vx;//returns the position of the point of closest approach to the PV 
        double Pho1_vy;
        double Pho1_vz;
        
        double Pho2_pt;
        double Pho2_Charge;
        double Pho2_PtDiff;
        double Pho2_ptfrac;
        double Pho2_eta;
        double Pho2_phi;
        double Pho2_vx;//returns the x-coordinate of vertex position
        double Pho2_vy;
        double Pho2_vz;

        double Pho12_Ptfrac; 
        double dRPho12; //deltaR between two highest PtTracks

        int nNeutralj1;

        double Neutral1_pt;
        double Neutral1_Charge;
        double Neutral1_PtDiff;
        double Neutral1_ptfrac;
        double Neutral1_eta;
        double Neutral1_phi;
        double Neutral1_vx;//returns the position of the point of closest approach to the PV 
        double Neutral1_vy;
        double Neutral1_vz;
        
        double Neutral2_pt;
        double Neutral2_Charge;
        double Neutral2_PtDiff;
        double Neutral2_ptfrac;
        double Neutral2_eta;
        double Neutral2_phi;
        double Neutral2_vx;//returns the x-coordinate of vertex position
        double Neutral2_vy;
        double Neutral2_vz;

};

AODtwoprong::AODtwoprong(const edm::ParameterSet& iConfig):
    vtxToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    jetsToken(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    PFCandidateToken(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
    TrackToken(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
    GenToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    pfMETsToken(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("MET"))),
    caloMETToken(consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("CaloMET")))
{
    edm::Service<TFileService> fs;

    RecoTree = fs->make<TTree>("RecoTree", "RecoTree");
    RecoTree->Branch("nvtx",&nvtx,"nvtx/I");

    //recoTrack Variables

    RecoTree->Branch("nPosTrkj1",&nPosTrkj1,"nPosTrkj1/I");
    RecoTree->Branch("nNegTrkj1",&nNegTrkj1,"nNegTrkj1/I"); 
    RecoTree->Branch("nTotTrkj1",&nTotTrkj1,"nTotTrkj1/I");
    RecoTree->Branch("posTrk1Pt",&posTrk1Pt,"posTrk1Pt/D");
    RecoTree->Branch("posTrk2Pt",&posTrk2Pt,"posTrk2Pt/D"); 
    RecoTree->Branch("posTrk3Pt",&posTrk3Pt,"posTrk3Pt/D");
    RecoTree->Branch("posTrk1Charge",&posTrk1Charge,"posTrk1Charge/D");
    RecoTree->Branch("posTrk2Charge",&posTrk2Charge,"posTrk2Charge/D");
    RecoTree->Branch("posTrk3Charge",&posTrk3Charge,"posTrk3Charge/D");
    
    RecoTree->Branch("negTrk1Pt",&negTrk1Pt,"negTrk1Pt/D");
    RecoTree->Branch("negTrk2Pt",&negTrk2Pt,"negTrk2Pt/D"); 
    RecoTree->Branch("negTrk3Pt",&negTrk3Pt,"negTrk3Pt/D");
    RecoTree->Branch("negTrk1Charge",&negTrk1Charge,"negTrk1Charge/D");
    RecoTree->Branch("negTrk2Charge",&negTrk2Charge,"negTrk2Charge/D");
    RecoTree->Branch("negTrk3Charge",&negTrk3Charge,"negTrk3Charge/D");
    RecoTree->Branch("dRPosTrk1negTrk1",&dRPosTrk1negTrk1,"dRPosTrk1negTrk1/D");

    RecoTree->Branch("pfMET", &pfMET);
    RecoTree->Branch("caloMET", &caloMET);
    RecoTree->Branch("j1Pt",&j1Pt,"j1Pt/D");
    RecoTree->Branch("j1Eta",&j1Eta,"j1Eta/D");
    RecoTree->Branch("j1Phi",&j1Phi,"j1Phi/D");
    RecoTree->Branch("j1CHF",       &j1CHF,"j1CHF/D");
    RecoTree->Branch("j1NHF",       &j1NHF,"j1NHF/D");
    RecoTree->Branch("j1ConsEtaPhiSpread",&j1ConsEtaPhiSpread,"j1ConsEtaPhiSpread/F");
    
    RecoTree->Branch("nChHadrPosj1",&nChHadrPosj1,"nChHadrPosj1/I");
    RecoTree->Branch("ChHadrPosrecoTrkj1",&ChHadrPosrecoTrkj1,"ChHadrPosrecoTrkj1/I");
    
    RecoTree->Branch("ChHadrPosTrk1Pt",&ChHadrPosTrk1Pt,"ChHadrPosTrk1Pt/D"); 
    RecoTree->Branch("ChHadrPosTrk1Eta",&ChHadrPosTrk1Eta,"ChHadrPosTrk1Eta/D");
    RecoTree->Branch("ChHadrPosTrk1Phi",&ChHadrPosTrk1Phi,"ChHadrPosTrk1Phi/D");
    RecoTree->Branch("ChHadrPosTrk1PtFrac",&ChHadrPosTrk1PtFrac,"ChHadrPosTrk1PtFrac/D");
    
    RecoTree->Branch("ChHadrPos1_pt",&ChHadrPos1_pt,"ChHadrPos1_pt/D");
    RecoTree->Branch("ChHadrPos1_Charge",&ChHadrPos1_Charge,"ChHadrPos1_Charge/D");
    RecoTree->Branch("ChHadrPos1_PtDiff",&ChHadrPos1_PtDiff,"ChHadrPos1_PtDiff/D");
    RecoTree->Branch("ChHadrPos1_ptfrac",&ChHadrPos1_ptfrac,"ChHadrPos1_ptfrac/D");
    RecoTree->Branch("ChHadrPos1_eta",&ChHadrPos1_eta,"ChHadrPos1_eta/D");
    RecoTree->Branch("ChHadrPos1_phi",&ChHadrPos1_phi,"ChHadrPos1_phi/D");
    RecoTree->Branch("ChHadrPos1_vx",&ChHadrPos1_vx,"ChHadrPos1_vx/D");
    RecoTree->Branch("ChHadrPos1_vy",&ChHadrPos1_vy,"ChHadrPos1_vy/D");
    RecoTree->Branch("ChHadrPos1_vz",&ChHadrPos1_vz,"ChHadrPos1_vz/D");

    RecoTree->Branch("ChHadrPosTrk2Pt",&ChHadrPosTrk2Pt,"ChHadrPosTrk2Pt/D"); 
    RecoTree->Branch("ChHadrPosTrk2Eta",&ChHadrPosTrk2Eta,"ChHadrPosTrk2Eta/D");
    RecoTree->Branch("ChHadrPosTrk2Phi",&ChHadrPosTrk2Phi,"ChHadrPosTrk2Phi/D");
    RecoTree->Branch("ChHadrPosTrk2PtFrac",&ChHadrPosTrk2PtFrac,"ChHadrPosTrk2PtFrac/D");
 
    RecoTree->Branch("ChHadrPos2_pt",&ChHadrPos2_pt,"ChHadrPos2_pt/D");
    RecoTree->Branch("ChHadrPos2_Charge",&ChHadrPos2_Charge,"ChHadrPos2_Charge/D");
    RecoTree->Branch("ChHadrPos2_PtDiff",&ChHadrPos2_PtDiff,"ChHadrPos2_PtDiff/D");
    RecoTree->Branch("ChHadrPos2_ptfrac",&ChHadrPos2_ptfrac,"ChHadrPos2_ptfrac/D");
    RecoTree->Branch("ChHadrPos2_eta",&ChHadrPos2_eta,"ChHadrPos2_eta/D");
    RecoTree->Branch("ChHadrPos2_phi",&ChHadrPos2_phi,"ChHadrPos2_phi/D");
    RecoTree->Branch("ChHadrPos2_vx",&ChHadrPos2_vx,"ChHadrPos2_vx/D");
    RecoTree->Branch("ChHadrPos2_vy",&ChHadrPos2_vy,"ChHadrPos2_vy/D");
    RecoTree->Branch("ChHadrPos2_vz",&ChHadrPos2_vz,"ChHadrPos2_vz/D");

    RecoTree->Branch("dRChHadrPosTrk12",&dRChHadrPosTrk12,"dRChHadrPosTrk12/D");

    RecoTree->Branch("ChHadrPosTrk12_Ptfrac",&ChHadrPosTrk12_Ptfrac,"ChHadrPosTrk12_Ptfrac/D"); 
    
    RecoTree->Branch("dRChHadrPos12",&dRChHadrPos12,"dRChHadrPos12/D");
    RecoTree->Branch("ChHadrPos12_Ptfrac",&ChHadrPos12_Ptfrac,"ChHadrPos12_Ptfrac/D"); 

    //Neg Charged Hadrons
    RecoTree->Branch("nChHadrNegj1",&nChHadrNegj1,"nChHadrNegj1/I");

    RecoTree->Branch("ChHadrNegrecoTrkj1",&ChHadrNegrecoTrkj1,"ChHadrNegrecoTrkj1/I");
    
    RecoTree->Branch("ChHadrNegTrk1Pt",&ChHadrNegTrk1Pt,"ChHadrNegTrk1Pt/D"); 
    RecoTree->Branch("ChHadrNegTrk1Eta",&ChHadrNegTrk1Eta,"ChHadrNegTrk1Eta/D");
    RecoTree->Branch("ChHadrNegTrk1Phi",&ChHadrNegTrk1Phi,"ChHadrNegTrk1Phi/D");
    RecoTree->Branch("ChHadrNegTrk1PtFrac",&ChHadrNegTrk1PtFrac,"ChHadrNegTrk1PtFrac/D");
    
    RecoTree->Branch("ChHadrNeg1_pt",&ChHadrNeg1_pt,"ChHadrNeg1_pt/D");
    RecoTree->Branch("ChHadrNeg1_Charge",&ChHadrNeg1_Charge,"ChHadrNeg1_Charge/D");
    RecoTree->Branch("ChHadrNeg1_PtDiff",&ChHadrNeg1_PtDiff,"ChHadrNeg1_PtDiff/D");
    RecoTree->Branch("ChHadrNeg1_ptfrac",&ChHadrNeg1_ptfrac,"ChHadrNeg1_ptfrac/D");
    RecoTree->Branch("ChHadrNeg1_eta",&ChHadrNeg1_eta,"ChHadrNeg1_eta/D");
    RecoTree->Branch("ChHadrNeg1_phi",&ChHadrNeg1_phi,"ChHadrNeg1_phi/D");
    RecoTree->Branch("ChHadrNeg1_vx",&ChHadrNeg1_vx,"ChHadrNeg1_vx/D");
    RecoTree->Branch("ChHadrNeg1_vy",&ChHadrNeg1_vy,"ChHadrNeg1_vy/D");
    RecoTree->Branch("ChHadrNeg1_vz",&ChHadrNeg1_vz,"ChHadrNeg1_vz/D");

    RecoTree->Branch("ChHadrNegTrk2Pt",&ChHadrNegTrk2Pt,"ChHadrNegTrk2Pt/D"); 
    RecoTree->Branch("ChHadrNegTrk2Eta",&ChHadrNegTrk2Eta,"ChHadrNegTrk2Eta/D");
    RecoTree->Branch("ChHadrNegTrk2Phi",&ChHadrNegTrk2Phi,"ChHadrNegTrk2Phi/D");
    RecoTree->Branch("ChHadrNegTrk2PtFrac",&ChHadrNegTrk2PtFrac,"ChHadrNegTrk2PtFrac/D");
 
    RecoTree->Branch("ChHadrNeg2_pt",&ChHadrNeg2_pt,"ChHadrNeg2_pt/D");
    RecoTree->Branch("ChHadrNeg2_Charge",&ChHadrNeg2_Charge,"ChHadrNeg2_Charge/D");
    RecoTree->Branch("ChHadrNeg2_PtDiff",&ChHadrNeg2_PtDiff,"ChHadrNeg2_PtDiff/D");
    RecoTree->Branch("ChHadrNeg2_ptfrac",&ChHadrNeg2_ptfrac,"ChHadrNeg2_ptfrac/D");
    RecoTree->Branch("ChHadrNeg2_eta",&ChHadrNeg2_eta,"ChHadrNeg2_eta/D");
    RecoTree->Branch("ChHadrNeg2_phi",&ChHadrNeg2_phi,"ChHadrNeg2_phi/D");
    RecoTree->Branch("ChHadrNeg2_vx",&ChHadrNeg2_vx,"ChHadrNeg2_vx/D");
    RecoTree->Branch("ChHadrNeg2_vy",&ChHadrNeg2_vy,"ChHadrNeg2_vy/D");
    RecoTree->Branch("ChHadrNeg2_vz",&ChHadrNeg2_vz,"ChHadrNeg2_vz/D");


    RecoTree->Branch("dRPosHadr1NegHadr1",&dRPosHadr1NegHadr1,"dRPosHadr1NegHadr1/D");

    RecoTree->Branch("dRChHadrNegTrk12",&dRChHadrNegTrk12,"dRChHadrNegTrk12/D");

    RecoTree->Branch("ChHadrNegTrk12_Ptfrac",&ChHadrNegTrk12_Ptfrac,"ChHadrNegTrk12_Ptfrac/D"); 
    
    RecoTree->Branch("dRChHadrNeg12",&dRChHadrNeg12,"dRChHadrNeg12/D");
    RecoTree->Branch("ChHadrNeg12_Ptfrac",&ChHadrNeg12_Ptfrac,"ChHadrNeg12_Ptfrac/D"); 
    //Photons from pi0

    RecoTree->Branch("nPhoj1",&nPhoj1,"nPhoj1/I");
    RecoTree->Branch("PhoTotalCharge",&PhoTotalCharge,"PhoTotalCharge/I");

    RecoTree->Branch("Pho1_pt",&Pho1_pt,"Pho1_pt/D"); 
    RecoTree->Branch("Pho1_Charge",&Pho1_Charge,"Pho1_Charge/D");
    RecoTree->Branch("Pho1_PtDiff",&Pho1_PtDiff,"Pho1_PtDiff/D");
    RecoTree->Branch("Pho1_ptfrac",&Pho1_ptfrac,"Pho1_ptfrac/D");
    RecoTree->Branch("Pho1_eta",&Pho1_eta,"Pho1_eta/D");
    RecoTree->Branch("Pho1_phi",&Pho1_phi,"Pho1_phi/D");
    RecoTree->Branch("Pho1_vx",&Pho1_vx,"Pho1_vx/D");
    RecoTree->Branch("Pho1_vy",&Pho1_vy,"Pho1_vy/D");
    RecoTree->Branch("Pho1_vz",&Pho1_vz,"Pho1_vz/D");
 
    RecoTree->Branch("Pho2_pt",&Pho2_pt,"Pho2_pt/D");
    RecoTree->Branch("Pho2_Charge",&Pho2_Charge,"Pho2_Charge/D");
    RecoTree->Branch("Pho2_PtDiff",&Pho2_PtDiff,"Pho2_PtDiff/D");
    RecoTree->Branch("Pho2_ptfrac",&Pho2_ptfrac,"Pho2_ptfrac/D");
    RecoTree->Branch("Pho2_eta",&Pho2_eta,"Pho2_eta/D");
    RecoTree->Branch("Pho2_phi",&Pho2_phi,"Pho2_phi/D");
    RecoTree->Branch("Pho2_vx",&Pho2_vx,"Pho2_vx/D");
    RecoTree->Branch("Pho2_vy",&Pho2_vy,"Pho2_vy/D");
    RecoTree->Branch("Pho2_vz",&Pho2_vz,"Pho2_vz/D");

    RecoTree->Branch("Pho12_Ptfrac",&Pho12_Ptfrac,"Pho12_Ptfrac/D");
    RecoTree->Branch("dRPho12",&dRPho12,"dRPho12/D");

    //Neutrals in PencilJet that are not Photons

    RecoTree->Branch("nNeutralj1",&nNeutralj1,"nNeutralj1/I");

    RecoTree->Branch("Neutral1_pt",&Neutral1_pt,"Neutral1_pt/D"); 
    RecoTree->Branch("Neutral1_Charge",&Neutral1_Charge,"Neutral1_Charge/D");
    RecoTree->Branch("Neutral1_PtDiff",&Neutral1_PtDiff,"Neutral1_PtDiff/D");
    RecoTree->Branch("Neutral1_ptfrac",&Neutral1_ptfrac,"Neutral1_ptfrac/D");
    RecoTree->Branch("Neutral1_eta",&Neutral1_eta,"Neutral1_eta/D");
    RecoTree->Branch("Neutral1_phi",&Neutral1_phi,"Neutral1_phi/D");
    RecoTree->Branch("Neutral1_vx",&Neutral1_vx,"Neutral1_vx/D");
    RecoTree->Branch("Neutral1_vy",&Neutral1_vy,"Neutral1_vy/D");
    RecoTree->Branch("Neutral1_vz",&Neutral1_vz,"Neutral1_vz/D");
 
    RecoTree->Branch("Neutral2_pt",&Neutral2_pt,"Neutral2_pt/D");
    RecoTree->Branch("Neutral2_Charge",&Neutral2_Charge,"Neutral2_Charge/D");
    RecoTree->Branch("Neutral2_PtDiff",&Neutral2_PtDiff,"Neutral2_PtDiff/D");
    RecoTree->Branch("Neutral2_ptfrac",&Neutral2_ptfrac,"Neutral2_ptfrac/D");
    RecoTree->Branch("Neutral2_eta",&Neutral2_eta,"Neutral2_eta/D");
    RecoTree->Branch("Neutral2_phi",&Neutral2_phi,"Neutral2_phi/D");
    RecoTree->Branch("Neutral2_vx",&Neutral2_vx,"Neutral2_vx/D");
    RecoTree->Branch("Neutral2_vy",&Neutral2_vy,"Neutral2_vy/D");
    RecoTree->Branch("Neutral2_vz",&Neutral2_vz,"Neutral2_vz/D");
}

AODtwoprong::~AODtwoprong()
{
}

void AODtwoprong::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken, vertices);
    nvtx=vertices->size();
    //const reco::Vertex &PV = vertices->front();
    
    edm::Handle<std::vector<reco::PFJet> > ak4jets;
    iEvent.getByToken(jetsToken, ak4jets);

    
    edm::Handle<std::vector<reco::GenParticle> > GenParticles;
    iEvent.getByToken(GenToken, GenParticles);
     
    edm::Handle<std::vector<reco::PFCandidate>> PFCandidates;
    iEvent.getByToken(PFCandidateToken, PFCandidates);

    edm::Handle<std::vector<reco::Track> > recoTracks;
    iEvent.getByToken(TrackToken,recoTracks);

    edm::Handle<std::vector<reco::PFMET> > pfMETs;
    iEvent.getByToken(pfMETsToken, pfMETs);


    edm::Handle<std::vector<reco::CaloMET> > caloMETs;
    iEvent.getByToken(caloMETToken, caloMETs);
    
    nPosTrkj1=0;
    nNegTrkj1=0;
    nTotTrkj1=0;
    posTrk1Pt=0;
    negTrk1Pt=0;
    posTrk1Charge=0;
    negTrk1Charge=0;

    posTrk2Pt=0;
    negTrk2Pt=0;
    posTrk2Charge=0;
    negTrk2Charge=0;

    posTrk3Pt=0;
    negTrk3Pt=0;
    posTrk3Charge=0;
    negTrk3Charge=0;

    dRPosTrk1negTrk1=0;
    
    //Charged Hadron Variables
    nChHadrPosj1=0;
    ChHadrPosrecoTrkj1=0;

    ChHadrPosTrk1Pt=0;
    ChHadrPosTrk1_Charge=0;
    ChHadrPosTrk1PtFrac=0;
    ChHadrPosTrk1Eta=0;
    ChHadrPosTrk1Phi=0;
    ChHadrPos1_pt=0;
    ChHadrPos1_Charge=0;
    ChHadrPos1_PtDiff=0;
    ChHadrPos1_ptfrac=0;
    ChHadrPos1_eta=0;
    ChHadrPos1_phi=0;
    ChHadrPos1_vx=0;//returns the position of the point of closest approach to the PV 
    ChHadrPos1_vy=0;
    ChHadrPos1_vz=0;


    ChHadrPosTrk2Pt=0; 
    ChHadrPosTrk2_Charge=0;
    ChHadrPosTrk2PtFrac=0;
    ChHadrPosTrk2Eta=0;
    ChHadrPosTrk2Phi=0;
    
    dRChHadrPosTrk12=0; //deltaR between two highest PtTracks
    
    ChHadrPos2_pt=0;
    ChHadrPos2_Charge=0;
    ChHadrPos2_PtDiff=0;
    ChHadrPos2_ptfrac=0;
    ChHadrPos2_eta=0;
    ChHadrPos2_phi=0;
    ChHadrPos2_vx=0;
    ChHadrPos2_vy=0;
    ChHadrPos2_vz=0;
    dRChHadrPos12=0; 
    ChHadrPos12_Ptfrac=0;
    ChHadrPosTrk12_Ptfrac=0;
    
    //NegChHadr Hadron Variables
    nChHadrNegj1=0;

    ChHadrNegrecoTrkj1=0;

    ChHadrNegTrk1Pt=0;
    ChHadrNegTrk1_Charge=0;
    ChHadrNegTrk1PtFrac=0;
    ChHadrNegTrk1Eta=0;
    ChHadrNegTrk1Phi=0;
    ChHadrNeg1_pt=0;
    ChHadrNeg1_Charge=0;
    ChHadrNeg1_PtDiff=0;
    ChHadrNeg1_ptfrac=0;
    ChHadrNeg1_eta=0;
    ChHadrNeg1_phi=0;
    ChHadrNeg1_vx=0;//returns the position of the point of closest approach to the PV 
    ChHadrNeg1_vy=0;
    ChHadrNeg1_vz=0;


    ChHadrNegTrk2Pt=0; 
    ChHadrNegTrk2_Charge=0;
    ChHadrNegTrk2PtFrac=0;
    ChHadrNegTrk2Eta=0;
    ChHadrNegTrk2Phi=0;
    
    dRChHadrNegTrk12=0; //deltaR between two highest PtTracks
    
    ChHadrNeg2_pt=0;
    ChHadrNeg2_Charge=0;
    ChHadrNeg2_PtDiff=0;
    ChHadrNeg2_ptfrac=0;
    ChHadrNeg2_eta=0;
    ChHadrNeg2_phi=0;
    ChHadrNeg2_vx=0;
    ChHadrNeg2_vy=0;
    ChHadrNeg2_vz=0;
    dRChHadrNeg12=0;
    ChHadrNeg12_Ptfrac=0;
    ChHadrNegTrk12_Ptfrac=0;



    //Photons in GenParticles
    nPhoj1=0;
    PhoTotalCharge=0;

    Pho1_pt=0;
    Pho1_Charge=0;
    Pho1_PtDiff=0;
    Pho1_ptfrac=0;
    Pho1_eta=0;
    Pho1_phi=0;
    Pho1_vx=0;
    Pho1_vy=0;
    Pho1_vz=0;
    
    Pho2_pt=0;
    Pho2_Charge=0;
    Pho2_PtDiff=0;
    Pho2_ptfrac=0;
    Pho2_eta=0;
    Pho2_phi=0;
    Pho2_vx=0;
    Pho2_vy=0;
    Pho2_vz=0;

    Pho12_Ptfrac=0; 
    dRPho12=0; 
    nNeutralj1=0;

    Neutral1_pt=0;
    Neutral1_Charge=0;
    Neutral1_PtDiff=0;
    Neutral1_ptfrac=0;
    Neutral1_eta=0;
    Neutral1_phi=0;
    Neutral1_vx=0;//returns the position of the point of closest approach to the PV 
    Neutral1_vy=0;
    Neutral1_vz=0;
    
    Neutral2_pt=0;
    Neutral2_Charge=0;
    Neutral2_PtDiff=0;
    Neutral2_ptfrac=0;
    Neutral2_eta=0;
    Neutral2_phi=0;
    Neutral2_vx=0;//returns the x-coordinate of vertex position
    Neutral2_vy=0;
    Neutral2_vz=0;
    //PhoTotalCharge = 0;
    int nZprime = 0;
    //dRPosTrk1negTrk1= 0=0;

    run_    = iEvent.id().run();
    event_  = iEvent.id().event();
    lumis_  = iEvent.luminosityBlock();
    
    std::cout<<"NewEvent"<<std::endl;
    std::cout<<run_<<":"<<event_<<":"<<lumis_<<std::endl;

    std::vector<const reco::GenParticle*> packedHadrons;
    for(std::vector<reco::GenParticle>::const_iterator genParticle = GenParticles->begin(); genParticle != GenParticles->end(); genParticle++){
        //std::cout<<"genParticle: "<<genParticle->pdgId()<<std::endl;
        if(TMath::Abs(genParticle->pdgId()) == 211) packedHadrons.push_back(&(*genParticle));
    }

    //Make list of possible PF charged pions to be used in tau reconstruction
    std::vector<const reco::PFCandidate*> hadronCandidates;
    for(std::vector<reco::PFCandidate>::const_iterator candidate= PFCandidates->begin(); candidate != PFCandidates->end(); candidate++){
        if(TMath::Abs(candidate->pdgId()) == 211) hadronCandidates.push_back(&(*candidate));
    }
    ChHadrPosrecoTrkj1 = 0; //count the number of reco Tracks within dR of leading Jet
    
    ChHadrNegrecoTrkj1 = 0; //count the number of reco Tracks within dR of leading Jet
    for(uint32_t i=0; i < ak4jets->size(); i++){
        const reco::PFJet &jet = (*ak4jets)[i];
        if (i==0){
        j1Pt = jet.pt();
        j1Eta = jet.eta();
        j1Phi = jet.phi();
        j1CHF = jet.chargedHadronEnergyFraction();
        j1NHF = jet.neutralHadronEnergyFraction();
        j1ConsEtaPhiSpread = jet.constituentEtaPhiSpread();
        double difference = 0;
        double PosTrkdifference = 0;
        double NegTrkdifference = 0;
        //These vectors will store the PtDifference between a packed PFCandidate and a nearbyJet
        std::vector<std::pair<double,const reco::PFCandidate*>> PtDiffChHadrPos;//pdgId = 211 
        std::vector<std::pair<double,const reco::PFCandidate*>> PtDiffChHadrNeg;//pdgId = -211
        std::vector<std::pair<double,const reco::PFCandidate*>> PtDiffPhotons;//pdgId = abs(22)
        std::vector<std::pair<double,const reco::PFCandidate*>> PtDiffNeutral; //miscellaneous neutral particles that are not photons
        std::vector<std::pair<double,const reco::Track*>>PtDiffPosTrk;//match the Trks to genlevel positive hadrons 
        std::vector<std::pair<double,const reco::Track*>>PtDiffNegTrk;//match the Trks to genlevel negative hadrons

        for(uint32_t j=0;j<GenParticles->size();j++){
          if(abs((*GenParticles)[j].pdgId())==600001){
            const reco::Candidate * Zprime = &(*GenParticles)[j];
            nZprime++;
            //Loop over all GenParticles
            if(nZprime<2){
            for(uint32_t k = 0; k < GenParticles->size(); k++) {
              const reco::GenParticle &genCand = (*GenParticles)[k]; 
              //Find ancestor of this genCand
              const reco::Candidate * genCand_Ancestor = &(*(genCand.mother(0)));
              //Now we confirm that if this genCand has Zprime ancestor, lets try matching it to appropriate PFCandidate
              if(genCand_Ancestor != nullptr && isAncestor(Zprime,genCand_Ancestor)){
                //match Tracks with appropriate GenCandidates
                float Trkdiff = 0; //diff between TrkCandPt and GenPt
                std::vector<std::pair<float,unsigned>> TrkPtdifference; //save index of the TrkCandidate that matches the genParticle from Zprime
                for(unsigned l = 0; l < recoTracks->size(); l++) {
                  const reco::Track &recoTrk = (*recoTracks)[l];
                  if(reco::deltaR(genCand.eta(),genCand.phi(),recoTrk.eta(),recoTrk.phi())<0.1){
                    Trkdiff=abs(genCand.pt()-recoTrk.pt());
                    TrkPtdifference.push_back(std::make_pair(Trkdiff,l));}}
                //Positive Tracks pdgId = +211 
                std::sort(TrkPtdifference.begin(),TrkPtdifference.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
                if(TrkPtdifference.size()>0 && (TrkPtdifference.at(0).first<20)){
                  unsigned n = TrkPtdifference.at(0).second;
                  const reco::Track &TrkCand = (*recoTracks)[n];
                  if(reco::deltaR(jet.eta(),jet.phi(),TrkCand.eta(),TrkCand.phi())< 0.2){
                     if(TrkCand.charge()==1 && genCand.charge()==1 && genCand.status()==1) { 
                      // std::cout<<"genCand matching with PosTrk: "<<genCand.pdgId()<<std::endl;
                      // std::cout<<"genCandPt: "<<genCand.pt()<<std::endl;
                      // std::cout<<"TrkPt: " <<TrkCand.pt()<<std::endl;
                      // std::cout<<"TrkPtdifference: "<<TrkPtdifference.at(0).first<<std::endl;
                      // std::cout<<"TrkCharge: "<<TrkCand.charge()<<std::endl;
                       PosTrkdifference = abs(jet.pt()-TrkCand.pt());
                       PtDiffPosTrk.push_back({PosTrkdifference,&TrkCand});}
                
                    //Negative Tracks with pdgId = -211
                     if(TrkCand.charge()==-1 && genCand.charge()== -1 && genCand.status()==1) { 
                      // std::cout<<"genCand matching with NegTrk: "<<genCand.pdgId()<<std::endl;
                      // std::cout<<"genCandPt: "<<genCand.pt()<<std::endl;
                      // std::cout<<"TrkPt: " <<TrkCand.pt()<<std::endl;
                      // std::cout<<"TrkPtdifference: "<<TrkPtdifference.at(0).first<<std::endl; 
                      // std::cout<<"TrkCharge: "<<TrkCand.charge()<<std::endl;
                       NegTrkdifference = abs(jet.pt()-TrkCand.pt());
                       PtDiffNegTrk.push_back({NegTrkdifference,&TrkCand});}
                  }
                }
                //We have to see which PFCandidate matches this GenCandidate so loop over all PFCandidates
                float diff = 0; //diff between PFCandPt and GenPt
                std::vector<std::pair<float,unsigned>> PtDifference; //save index of the PFCandidate that matches the genParticle from Zprime
                for(unsigned l = 0; l < PFCandidates->size(); l++) {
                  const reco::PFCandidate &pfCand = (*PFCandidates)[l];
                  if(reco::deltaR(genCand.eta(),genCand.phi(),pfCand.eta(),pfCand.phi())<0.1){
                    diff=abs(genCand.pt()-pfCand.pt());
                    PtDifference.push_back(std::make_pair(diff,l));}}
                //std::cout<<"No.of PFCands matched to GenParticle: "<<PtDifference.size()<<std::endl; 
                std::sort(PtDifference.begin(),PtDifference.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
                for(unsigned m=0;m<PtDifference.size();m++){
                  if(PtDifference.at(m).first<20){
                     unsigned q = PtDifference.at(m).second;
                     const reco::PFCandidate &PFCand = (*PFCandidates)[q];
                     if(reco::deltaR(jet.eta(),jet.phi(),PFCand.eta(),PFCand.phi())< 0.2){
                        if((PFCand.pdgId()== 211 || PFCand.pdgId()==321) && genCand.charge()==1 && genCand.status()==1) { 
                          //std::cout<<"genCand matching with ChHadr: "<<abs(genCand.pdgId())<<std::endl;
                          //std::cout<<"PtDifference: "<<PtDifference.at(0).first<<std::endl;
                          difference = abs(jet.pt()-PFCand.pt());
                          PtDiffChHadrPos.push_back({difference,&PFCand});}
                          
                        if((PFCand.pdgId()== -211 || PFCand.pdgId()==-321) && genCand.charge()==-1 && genCand.status()==1) {
                          //std::cout<<"genCand matching with ChHadrNeg: "<<abs(genCand.pdgId())<<std::endl;
                          //std::cout<<"PtDifference: "<<PtDifference.at(0).first<<std::endl;
                          difference = abs(jet.pt()-PFCand.pt());
                          PtDiffChHadrNeg.push_back(std::make_pair(difference,&PFCand));}

                        if(abs(PFCand.pdgId())== 22 && genCand.charge()==0 && genCand.status()==1) {
                          //std::cout<<"genCand matching with recoPhotons: "<<abs(genCand.pdgId())<<std::endl;
                          //std::cout<<"PtDifference: "<<PtDifference.at(0).first<<std::endl;
                          difference = abs(jet.pt()-PFCand.pt());
                          PtDiffPhotons.push_back(std::make_pair(difference,&PFCand));}

                        //miscellaneous neutral particles that are mostly K-Long mesons and since K-long decays to pions, leptons, neutrinos, even pi0s
                        //This genMatching is tricky because it could match to a whole bunch of genLevel particles
                        if((PFCand.charge()==0 && PFCand.pdgId()==130) && genCand.status()==1){
                          std::cout<<"Neutral PFCand: "<<PFCand.pdgId()<<std::endl;
                          std::cout<<"genCand matching with neutral particles in PencilJet: "<<abs(genCand.pdgId())<<std::endl;
                          std::cout<<"PtDifference: "<<PtDifference.at(0).first<<std::endl;
                          difference = abs(jet.pt()-PFCand.pt());
                          PtDiffNeutral.push_back(std::make_pair(difference,&PFCand));
                        }
      }

                }
              }//closing loop for PtDifference which PFCandidates with PtDiff from genCandidate in question
                  }

              }}}}   
    
        //Sort these PtDifference vectors by ascending order in PtDifference from the Jet
          std::sort(PtDiffChHadrPos.begin(),PtDiffChHadrPos.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
          //std::sort(PtDiffChHadr.begin(),PtDiffChHadr.end(),pairCompare);
          std::sort(PtDiffChHadrNeg.begin(),PtDiffChHadrNeg.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;}); 
          std::sort(PtDiffPhotons.begin(),PtDiffPhotons.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
          
          std::sort(PtDiffNeutral.begin(),PtDiffNeutral.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
          //sort the tracks

          std::sort(PtDiffPosTrk.begin(),PtDiffPosTrk.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;}); 
          std::sort(PtDiffNegTrk.begin(),PtDiffNegTrk.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
          
          //Get track story
          nPosTrkj1=PtDiffPosTrk.size();
          nNegTrkj1=PtDiffNegTrk.size();
          nTotTrkj1=nPosTrkj1+nNegTrkj1;
          //Positive Tracks
          for(int i=0;i<nPosTrkj1;i++){
            const reco::Track &PosTrk = *(PtDiffPosTrk.at(i).second);
            if(i==0){
              posTrk1Pt = PosTrk.pt();
              posTrk1Charge = PosTrk.charge();
            }

            if(i==1){
              posTrk2Pt = PosTrk.pt();
              posTrk2Charge = PosTrk.charge();
            }
            if(i==2){
              posTrk3Pt = PosTrk.pt();
              posTrk3Charge = PosTrk.charge();
            }
            if(i==2){ 
              break;
            }
          }
          //Negative Tracks
          for(int i=0;i<nNegTrkj1;i++){
            const reco::Track &NegTrk = *(PtDiffNegTrk.at(i).second);
            if(i==0){
              negTrk1Pt = NegTrk.pt();
              negTrk1Charge = NegTrk.charge();
            }

            if(i==1){
              negTrk2Pt = NegTrk.pt();
              negTrk2Charge = NegTrk.charge();
            }
            if(i==2){
              negTrk3Pt = NegTrk.pt();
              negTrk3Charge = NegTrk.charge();
            }
            if(i==2){ 
              break;
            }
          }
          if(nPosTrkj1 >0 && nNegTrkj1>0){
            const reco::Track &PosTrk1 = *(PtDiffPosTrk.at(0).second);
            const reco::Track &NegTrk1 = *(PtDiffNegTrk.at(0).second);
            dRPosTrk1negTrk1=deltaR(PosTrk1.eta(),PosTrk1.phi(),NegTrk1.eta(),NegTrk1.phi());
          }
          //PosCharged Hadron begin
          nChHadrPosj1 = PtDiffChHadrPos.size();
          for (int i=0;i<nChHadrPosj1;i++){
            const reco::PFCandidate &ChHadrPos = *(PtDiffChHadrPos.at(i).second);
            //std::cout<<"ChHadr: "<<ChHadr.charge();
            if(ChHadrPos.bestTrack()!=nullptr){//method points to null if there is no track 
              ChHadrPosrecoTrkj1++;
              //std::cout<<"TrkCharge: "<<Trk.charge()<<std::endl;
            }}
          if(nChHadrPosj1>0){
            const reco::PFCandidate &ChHadrPos1 = *(PtDiffChHadrPos.at(0).second);
            if(ChHadrPos1.bestTrack()!=nullptr){//lets store some info about Trk1 associated with highestPt Hadron
              const reco::Track &Trk = *(ChHadrPos1.bestTrack());
              ChHadrPosTrk1Pt = Trk.pt();
              ChHadrPosTrk1_Charge = Trk.charge();
              ChHadrPosTrk1PtFrac = (ChHadrPosTrk1Pt/jet.pt());
              ChHadrPosTrk1Eta = Trk.eta();
              ChHadrPosTrk1Phi = Trk.phi();
              }
            ChHadrPos1_pt = ChHadrPos1.pt();
            ChHadrPos1_Charge = ChHadrPos1.charge();
            ChHadrPos1_PtDiff = PtDiffChHadrPos.at(0).first;
            ChHadrPos1_ptfrac = (ChHadrPos1_pt/jet.pt());
            ChHadrPos1_eta = ChHadrPos1.eta();
            ChHadrPos1_phi = ChHadrPos1.phi();
            ChHadrPos1_vx = ChHadrPos1.vx();
            ChHadrPos1_vy = ChHadrPos1.vy();
            ChHadrPos1_vz = ChHadrPos1.vz();
          }//Close the loop for first Charged Hadron 
          //Start to look at second charged Hadron in the event  
          if(PtDiffChHadrPos.size()>=2){
            const reco::PFCandidate &ChHadrPos2 = *(PtDiffChHadrPos.at(1).second);
            if(ChHadrPos2.bestTrack()!=nullptr){//lets store some info about Trk2 associated with second highestPt charged Hadron
              const reco::Track &Trk = *(ChHadrPos2.bestTrack());
              ChHadrPosTrk2Pt = Trk.pt();
              ChHadrPosTrk2_Charge = Trk.charge();
              ChHadrPosTrk2PtFrac = (ChHadrPosTrk2Pt/jet.pt());
              ChHadrPosTrk2Eta = Trk.eta();
              ChHadrPosTrk2Phi = Trk.phi();
              }
            ChHadrPos2_pt = ChHadrPos2.pt();
            ChHadrPos2_Charge = ChHadrPos2.charge();
            ChHadrPos2_PtDiff = PtDiffChHadrPos.at(1).first;
            ChHadrPos2_ptfrac = (ChHadrPos2_pt/jet.pt());
            ChHadrPos2_eta = ChHadrPos2.eta();
            ChHadrPos2_phi = ChHadrPos2.phi();
            ChHadrPos2_vx = ChHadrPos2.vx();
            ChHadrPos2_vy = ChHadrPos2.vy();
            ChHadrPos2_vz = ChHadrPos2.vz();
            } 
            dRChHadrPosTrk12 = reco::deltaR(ChHadrPosTrk1Eta,ChHadrPosTrk1Phi,ChHadrPosTrk2Eta,ChHadrPosTrk2Phi);//Reco ChargeHadron analysis ends
            
            dRChHadrPos12 = reco::deltaR(ChHadrPos1_eta,ChHadrPos1_phi,ChHadrPos2_eta,ChHadrPos2_phi);//Reco ChargeHadron analysis ends
            ChHadrPos12_Ptfrac = (ChHadrPos1_pt+ChHadrPos2_pt)/(j1Pt);

            ChHadrPosTrk12_Ptfrac = (ChHadrPosTrk1Pt+ChHadrPosTrk2Pt)/(j1Pt);
      
            //Neg Charged Hadrons
            nChHadrNegj1 = PtDiffChHadrNeg.size();
            
            for (int i=0;i<nChHadrNegj1;i++){
              const reco::PFCandidate &ChHadrNeg = *(PtDiffChHadrNeg.at(i).second);
              //std::cout<<"ChHadr: "<<ChHadr.charge();
              if(ChHadrNeg.bestTrack()!=nullptr){//method points to null if there is no track 
                ChHadrNegrecoTrkj1++;
                //std::cout<<"TrkCharge: "<<Trk.charge()<<std::endl;
              }}
            if(nChHadrNegj1>0){
              const reco::PFCandidate &ChHadrNeg1 = *(PtDiffChHadrNeg.at(0).second);
              if(ChHadrNeg1.bestTrack()!=nullptr){//lets store some info about Trk1 associated with highestPt Hadron
                const reco::Track &Trk = *(ChHadrNeg1.bestTrack());
                ChHadrNegTrk1Pt = Trk.pt();
                ChHadrNegTrk1_Charge = Trk.charge();
                ChHadrNegTrk1PtFrac = (ChHadrNegTrk1Pt/jet.pt());
                ChHadrNegTrk1Eta = Trk.eta();
                ChHadrNegTrk1Phi = Trk.phi();
                }
              ChHadrNeg1_pt = ChHadrNeg1.pt();
              ChHadrNeg1_Charge = ChHadrNeg1.charge();
              ChHadrNeg1_PtDiff = PtDiffChHadrNeg.at(0).first;
              ChHadrNeg1_ptfrac = (ChHadrNeg1_pt/jet.pt());
              ChHadrNeg1_eta = ChHadrNeg1.eta();
              ChHadrNeg1_phi = ChHadrNeg1.phi();
              ChHadrNeg1_vx = ChHadrNeg1.vx();
              ChHadrNeg1_vy = ChHadrNeg1.vy();
              ChHadrNeg1_vz = ChHadrNeg1.vz();
            }//Close the loop for first Charged Hadron 
            //Start to look at second charged Hadron in the event  
            if(PtDiffChHadrNeg.size()>=2){
              const reco::PFCandidate &ChHadrNeg2 = *(PtDiffChHadrNeg.at(1).second);
              if(ChHadrNeg2.bestTrack()!=nullptr){//lets store some info about Trk2 associated with second highestPt charged Hadron
                const reco::Track &Trk = *(ChHadrNeg2.bestTrack());
                ChHadrNegTrk2Pt = Trk.pt();
                ChHadrNegTrk2_Charge = Trk.charge();
                ChHadrNegTrk2PtFrac = (ChHadrNegTrk2Pt/jet.pt());
                ChHadrNegTrk2Eta = Trk.eta();
                ChHadrNegTrk2Phi = Trk.phi();
                }
              ChHadrNeg2_pt = ChHadrNeg2.pt();
              ChHadrNeg2_Charge = ChHadrNeg2.charge();
              ChHadrNeg2_PtDiff = PtDiffChHadrNeg.at(1).first;
              ChHadrNeg2_ptfrac = (ChHadrNeg2_pt/jet.pt());
              ChHadrNeg2_eta = ChHadrNeg2.eta();
              ChHadrNeg2_phi = ChHadrNeg2.phi();
              ChHadrNeg2_vx = ChHadrNeg2.vx();
              ChHadrNeg2_vy = ChHadrNeg2.vy();
              ChHadrNeg2_vz = ChHadrNeg2.vz();
              } 
              dRChHadrNegTrk12 = reco::deltaR(ChHadrNegTrk1Eta,ChHadrNegTrk1Phi,ChHadrNegTrk2Eta,ChHadrNegTrk2Phi);//Reco ChargeHadron analysis ends
              
              dRChHadrNeg12 = reco::deltaR(ChHadrNeg1_eta,ChHadrNeg1_phi,ChHadrNeg2_eta,ChHadrNeg2_phi);//Reco ChargeHadron analysis ends
              ChHadrNeg12_Ptfrac = (ChHadrNeg1_pt+ChHadrNeg2_pt)/(j1Pt);

              ChHadrNegTrk12_Ptfrac = (ChHadrNegTrk1Pt+ChHadrNegTrk2Pt)/(j1Pt);

            //detaR between leading positive and negative hadron 
            dRPosHadr1NegHadr1=deltaR(ChHadrPos1_eta,ChHadrPos1_phi,ChHadrNeg1_eta,ChHadrNeg1_phi);

            //RecoPhotons
            nPhoj1 = PtDiffPhotons.size();
            for (int i=0;i<nPhoj1;i++){
              const reco::PFCandidate &Pho = *(PtDiffPhotons.at(i).second);
              PhoTotalCharge = PhoTotalCharge + Pho.charge();
            }
            if(nPhoj1>0){
              const reco::PFCandidate &Pho1 = *(PtDiffPhotons.at(0).second);
              Pho1_pt = Pho1.pt();
              Pho1_Charge = Pho1.charge();
              Pho1_PtDiff = PtDiffPhotons.at(0).first;
              Pho1_ptfrac = (Pho1_pt/jet.pt());
              Pho1_eta = Pho1.eta();
              Pho1_phi = Pho1.phi();
              Pho1_vx = Pho1.vx();
              Pho1_vy = Pho1.vy();
              Pho1_vz = Pho1.vz();
            }
              //std::cout<<"Where is this error?"<<std::endl;
              if(PtDiffPhotons.size()>=2){
              const reco::PFCandidate &Pho2 = *(PtDiffPhotons.at(1).second);
              Pho2_pt = Pho2.pt();
              Pho2_Charge = Pho2.charge();
              Pho2_PtDiff = PtDiffPhotons.at(1).first;
              Pho2_ptfrac = (Pho2_pt/jet.pt());
              Pho2_eta = Pho2.eta();
              Pho2_phi = Pho2.phi();
              Pho2_vx = Pho2.vx();
              Pho2_vy = Pho2.vy();
              Pho2_vz = Pho2.vz();
              } 
              dRPho12 = reco::deltaR(Pho1_eta,Pho1_phi,Pho2_eta,Pho2_phi); //Reco Neutral Hadron Analysis ends 
              Pho12_Ptfrac = (Pho1_pt+Pho2_pt)/(j1Pt);

            //Neutral Particles inside the PencilJet

            nNeutralj1 = PtDiffNeutral.size();
            std::cout<<"No.of Neutral particles in PencilJet that are not Photons: "<<nNeutralj1<<std::endl;
            if(nNeutralj1>0){
              const reco::PFCandidate &Neutral1 = *(PtDiffNeutral.at(0).second);
              Neutral1_pt = Neutral1.pt();
              Neutral1_Charge = Neutral1.charge();
              Neutral1_PtDiff = PtDiffNeutral.at(0).first;
              Neutral1_ptfrac = (Neutral1_pt/jet.pt());
              Neutral1_eta = Neutral1.eta();
              Neutral1_phi = Neutral1.phi();
              Neutral1_vx = Neutral1.vx();
              Neutral1_vy = Neutral1.vy();
              Neutral1_vz = Neutral1.vz();
            }
              //std::cout<<"Where is this error?"<<std::endl;
              if(PtDiffNeutral.size()>=2){
              const reco::PFCandidate &Neutral2 = *(PtDiffNeutral.at(1).second);
              Neutral2_pt = Neutral2.pt();
              Neutral2_Charge = Neutral2.charge();
              Neutral2_PtDiff = PtDiffNeutral.at(1).first;
              Neutral2_ptfrac = (Neutral2_pt/jet.pt());
              Neutral2_eta = Neutral2.eta();
              Neutral2_phi = Neutral2.phi();
              Neutral2_vx = Neutral2.vx();
              Neutral2_vy = Neutral2.vy();
              Neutral2_vz = Neutral2.vz();
              } 
          } //closing the leading jet
        }//closing the ak4jets loop
           
            pfMET = -99,caloMET = -99;
            pfMET = (*pfMETs)[0].pt();
            caloMET = (*caloMETs)[0].pt();

            RecoTree->Fill();
}
bool AODtwoprong::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
        if(ancestor == particle ) return true;
 
        //otherwise loop on mothers, if any and return true if the ancestor is found
        for(size_t i=0;i< particle->numberOfMothers();i++)
          {
              if(isAncestor(ancestor,particle->mother(i))) return true;
          }
          //if we did not return yet, then particle and ancestor are not relatives
        return false;
}


// ------------ method called once each job just before starting event loop  ------------
void 
AODtwoprong::beginJob()
{
}


// ------------ method called once each job just after ending the event loop  ------------
	void 
AODtwoprong::endJob() 
{
}

//void AODtwoprong::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
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
//int AODtwoprong::GetDecayMode(std::vector<const reco::GenParticle*>& daughters){
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
//bool AODtwoprong::isNeutrino(const reco::Candidate* daughter)
//{
//  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
//}
////Gets visible 4-momentum of a particle from list of daughters
//reco::Candidate::LorentzVector AODtwoprong::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
//	reco::Candidate::LorentzVector p4_vis(0,0,0,0);	
// 	for(size_t i = 0; i < daughters.size(); ++i){
// 		if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
// 			p4_vis += daughters[i]->p4();
//		}
// 	}
//	return p4_vis;
//}
//
//void AODtwoprong::eraseHadronCands(std::vector<const reco::PFCandidate*>& hadronCands, reco::CandidatePtrVector signalCands){
//    for(size_t i=0; i<signalCands.size(); ++i){
//        for(size_t j=0; j<hadronCands.size(); ++j){
//            if(signalCands[i]->pt()==hadronCands[j]->pt()) hadronCands.erase(hadronCands.begin()+j);
//        }
//    }
//}
////Finds closest charged pion not already used in 2-prong tau reconstruction
//const reco::PFCandidate* AODtwoprong::findThirdHadron(std::vector<const reco::PFCandidate*> hadronCands, reco::CandidatePtrVector signalCands, std::vector<const reco::GenParticle*> daughters,int *recoTrack, const reco::GenParticle* *genTrack3, double *trackDR){
//    const reco::PFCandidate* thirdHadron = hadronCands[0];
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
AODtwoprong::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	  edm::ParameterSetDescription desc;
	  desc.setUnknown();
	  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(AODtwoprong);
