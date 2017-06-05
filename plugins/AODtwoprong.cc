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

        TTree* RecoTree;
        int     run_;
        int  event_;
        int     lumis_;
        int nvtx;
        
        double j1Pt;
        double j1Eta;
        double j1Phi;
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
        //Charged Hadron Variables
        int nChHadrj1;
        double ChHadrTotalCharge; 
        int recoTrkj1;
        double TrkTotalCharge;

        double Trk1Pt;
        double Trk1_Charge;
        double Trk1PosCharge; 
        double Trk1NegCharge;
        int nTrk1Pos;
        int nTrk1Neg;
        double Trk1PtFrac;
        double Trk1Eta;
        double Trk1Phi;
        double ChHadr1_pt;
        double ChHadr1_Charge;
        double ChHadr1_PtDiff;
        double ChHadr1_ptfrac;
        double ChHadr1_eta;
        double ChHadr1_phi;
        double ChHadr1_vx;//returns the position of the point of closest approach to the PV 
        double ChHadr1_vy;
        double ChHadr1_vz;


        double Trk2Pt; 
        double Trk2_Charge;
        double Trk2PosCharge; 
        double Trk2NegCharge;
        int nTrk2Pos;
        int nTrk2Neg;
        double Trk2PtFrac;
        double Trk2Eta;
        double Trk2Phi;
        
        double dRTrk12; //deltaR between two highest PtTracks
        
        double ChHadr2_pt;
        double ChHadr2_Charge;
        double ChHadr2_PtDiff;
        double ChHadr2_ptfrac;
        double ChHadr2_eta;
        double ChHadr2_phi;
        double ChHadr2_vx;//returns the x-coordinate of vertex position
        double ChHadr2_vy;
        double ChHadr2_vz;

        double dRChHadr12; //deltaR between two highest PtTracks
        double ChHadr12_Ptfrac;
        int ChHadr12_Charge;

        double Trk12_Ptfrac;
        int Trk12_Charge;

        //Neutral Hadron Variables
        int nNeutrHadrj1;

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
};

AODtwoprong::AODtwoprong(const edm::ParameterSet& iConfig):
    vtxToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    jetsToken(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    PFCandidateToken(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
    TrackToken(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
    GenToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned")))
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


    RecoTree->Branch("j1Pt",&j1Pt,"j1Pt/D");
    RecoTree->Branch("j1Eta",&j1Eta,"j1Eta/D");
    RecoTree->Branch("j1Phi",&j1Phi,"j1Phi/D");
    RecoTree->Branch("j1ConsEtaPhiSpread",&j1ConsEtaPhiSpread,"j1ConsEtaPhiSpread/F");
    
    RecoTree->Branch("nChHadrj1",&nChHadrj1,"nChHadrj1/I");
    RecoTree->Branch("TrkTotalCharge",&TrkTotalCharge,"TrkTotalCharge/D");
    RecoTree->Branch("recoTrkj1",&recoTrkj1,"recoTrkj1/I");
    RecoTree->Branch("ChHadrTotalCharge",&ChHadrTotalCharge,"ChHadrTotalCharge/D");
    
    RecoTree->Branch("Trk1Pt",&Trk1Pt,"Trk1Pt/D"); 
    RecoTree->Branch("Trk1PosCharge",&Trk1PosCharge,"Trk1PosCharge/D"); 
    RecoTree->Branch("Trk1NegCharge",&Trk1NegCharge,"Trk1NegCharge/D");
    RecoTree->Branch("nTrk1Pos",&nTrk1Pos,"nTrk1Pos/I"); 
    RecoTree->Branch("nTrk1Neg",&nTrk1Neg,"nTrk1Neg/I");
    RecoTree->Branch("Trk1Eta",&Trk1Eta,"Trk1Eta/D");
    RecoTree->Branch("Trk1Phi",&Trk1Phi,"Trk1Phi/D");
    RecoTree->Branch("Trk1PtFrac",&Trk1PtFrac,"Trk1PtFrac/D");
    
    RecoTree->Branch("ChHadr1_pt",&ChHadr1_pt,"ChHadr1_pt/D");
    RecoTree->Branch("ChHadr1_Charge",&ChHadr1_Charge,"ChHadr1_Charge/D");
    RecoTree->Branch("ChHadr1_PtDiff",&ChHadr1_PtDiff,"ChHadr1_PtDiff/D");
    RecoTree->Branch("ChHadr1_ptfrac",&ChHadr1_ptfrac,"ChHadr1_ptfrac/D");
    RecoTree->Branch("ChHadr1_eta",&ChHadr1_eta,"ChHadr1_eta/D");
    RecoTree->Branch("ChHadr1_phi",&ChHadr1_phi,"ChHadr1_phi/D");
    RecoTree->Branch("ChHadr1_vx",&ChHadr1_vx,"ChHadr1_vx/D");
    RecoTree->Branch("ChHadr1_vy",&ChHadr1_vy,"ChHadr1_vy/D");
    RecoTree->Branch("ChHadr1_vz",&ChHadr1_vz,"ChHadr1_vz/D");

    RecoTree->Branch("Trk2Pt",&Trk2Pt,"Trk2Pt/D"); 
    RecoTree->Branch("Trk2PosCharge",&Trk2PosCharge,"Trk2PosCharge/D"); 
    RecoTree->Branch("Trk2NegCharge",&Trk2NegCharge,"Trk2NegCharge/D");
    RecoTree->Branch("nTrk2Pos",&nTrk2Pos,"nTrk2Pos/I"); 
    RecoTree->Branch("nTrk2Neg",&nTrk2Neg,"nTrk2Neg/I");
    RecoTree->Branch("Trk2Eta",&Trk2Eta,"Trk2Eta/D");
    RecoTree->Branch("Trk2Phi",&Trk2Phi,"Trk2Phi/D");
    RecoTree->Branch("Trk2PtFrac",&Trk2PtFrac,"Trk2PtFrac/D");
 
    RecoTree->Branch("ChHadr2_pt",&ChHadr2_pt,"ChHadr2_pt/D");
    RecoTree->Branch("ChHadr2_Charge",&ChHadr2_Charge,"ChHadr2_Charge/D");
    RecoTree->Branch("ChHadr2_PtDiff",&ChHadr2_PtDiff,"ChHadr2_PtDiff/D");
    RecoTree->Branch("ChHadr2_ptfrac",&ChHadr2_ptfrac,"ChHadr2_ptfrac/D");
    RecoTree->Branch("ChHadr2_eta",&ChHadr2_eta,"ChHadr2_eta/D");
    RecoTree->Branch("ChHadr2_phi",&ChHadr2_phi,"ChHadr2_phi/D");
    RecoTree->Branch("ChHadr2_vx",&ChHadr2_vx,"ChHadr2_vx/D");
    RecoTree->Branch("ChHadr2_vy",&ChHadr2_vy,"ChHadr2_vy/D");
    RecoTree->Branch("ChHadr2_vz",&ChHadr2_vz,"ChHadr2_vz/D");

    RecoTree->Branch("dRTrk12",&dRTrk12,"dRTrk12/D");

    RecoTree->Branch("Trk12_Ptfrac",&Trk12_Ptfrac,"Trk12_Ptfrac/D"); 
    RecoTree->Branch("Trk12_Charge",&Trk12_Charge,"Trk12_Charge/I");
    
    RecoTree->Branch("dRChHadr12",&dRChHadr12,"dRChHadr12/D");
    RecoTree->Branch("ChHadr12_Ptfrac",&ChHadr12_Ptfrac,"ChHadr12_Ptfrac/D"); 
    RecoTree->Branch("ChHadr12_Charge",&ChHadr12_Charge,"ChHadr12_Charge/I");

    //Neutral Hadrons
    RecoTree->Branch("nNeutrHadrj1",&nNeutrHadrj1,"nNeutrHadrj1/I");

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
    
    ChHadrTotalCharge = -50;
    TrkTotalCharge = -50;
    Trk1PosCharge = -50;
    Trk1NegCharge = -50; 
    Trk2PosCharge = -50;
    Trk2NegCharge = -50;
    nTrk1Pos = nTrk1Neg = 0;
    nTrk2Pos = nTrk2Neg = 0;
    ChHadr1_pt = ChHadr2_pt = Trk1Pt = Trk2Pt = 0;
    PhoTotalCharge = 0;
    int nZprime = 0;

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
    recoTrkj1 = 0; //count the number of reco Tracks within dR of leading Jet
    for(uint32_t i=0; i < ak4jets->size(); i++){
        const reco::PFJet &jet = (*ak4jets)[i];
        if (i==0){
        j1Pt = jet.pt();
        j1Eta = jet.eta();
        j1Phi = jet.phi();
        j1ConsEtaPhiSpread = jet.constituentEtaPhiSpread();
        double difference = 0;
        double PosTrkdifference = 0;
        double NegTrkdifference = 0;
        //These vectors will store the PtDifference between a packed PFCandidate and a nearbyJet
        std::vector<std::pair<double,const reco::PFCandidate*>> PtDiffChHadr;//pdgId = abs(211) 
        std::vector<std::pair<double,const reco::PFCandidate*>> PtDiffNeutrHadr;//pdgId = abs(130)
        std::vector<std::pair<double,const reco::PFCandidate*>> PtDiffPhotons;//pdgId = abs(22)
        
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
                if(TrkPtdifference.size()>0 && (TrkPtdifference.at(0).first<10)){
                  unsigned n = TrkPtdifference.at(0).second;
                  const reco::Track &TrkCand = (*recoTracks)[n];
                  if(reco::deltaR(jet.eta(),jet.phi(),TrkCand.eta(),TrkCand.phi())< 0.2){
                     if(TrkCand.charge()==1 && genCand.charge()==1 && genCand.status()==1) { 
                       std::cout<<"genCand matching with PosTrk: "<<genCand.pdgId()<<std::endl;
                       std::cout<<"genCandPt: "<<genCand.pt()<<std::endl;
                       std::cout<<"TrkPt: " <<TrkCand.pt()<<std::endl;
                       std::cout<<"TrkPtdifference: "<<TrkPtdifference.at(0).first<<std::endl;
                       std::cout<<"TrkCharge: "<<TrkCand.charge()<<std::endl;
                       PosTrkdifference = abs(jet.pt()-TrkCand.pt());
                       PtDiffPosTrk.push_back({PosTrkdifference,&TrkCand});}
                
                    //Negative Tracks with pdgId = -211
                     if(TrkCand.charge()==-1 && genCand.charge()== -1 && genCand.status()==1) { 
                       std::cout<<"genCand matching with NegTrk: "<<genCand.pdgId()<<std::endl;
                       std::cout<<"genCandPt: "<<genCand.pt()<<std::endl;
                       std::cout<<"TrkPt: " <<TrkCand.pt()<<std::endl;
                       std::cout<<"TrkPtdifference: "<<TrkPtdifference.at(0).first<<std::endl; 
                       std::cout<<"TrkCharge: "<<TrkCand.charge()<<std::endl;
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
                
                std::sort(PtDifference.begin(),PtDifference.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
                if(PtDifference.size()>0 && (PtDifference.at(0).first<20)){
                  unsigned m = PtDifference.at(0).second;
                  const reco::PFCandidate &PFCand = (*PFCandidates)[m];
                  if(reco::deltaR(jet.eta(),jet.phi(),PFCand.eta(),PFCand.phi())< 0.2){
                     if(abs(PFCand.pdgId())== 211) { 
                       //std::cout<<"genCand matching with ChHadr: "<<abs(genCand.pdgId())<<std::endl;
                       //std::cout<<"PtDifference: "<<PtDifference.at(0).first<<std::endl;
                       difference = abs(jet.pt()-PFCand.pt());
                       PtDiffChHadr.push_back({difference,&PFCand});}

                     if(abs(PFCand.pdgId())== 130) {
                       //std::cout<<"genCand matching with NeutrHadr: "<<abs(genCand.pdgId())<<std::endl;
                       //std::cout<<"PtDifference: "<<PtDifference.at(0).first<<std::endl;
                       difference = abs(jet.pt()-PFCand.pt());
                       PtDiffNeutrHadr.push_back(std::make_pair(difference,&PFCand));}

                     if(abs(PFCand.pdgId())== 22) {
                       //std::cout<<"genCand matching with recoPhotons: "<<abs(genCand.pdgId())<<std::endl;
                       //std::cout<<"PtDifference: "<<PtDifference.at(0).first<<std::endl;
                       difference = abs(jet.pt()-PFCand.pt());
                       PtDiffPhotons.push_back(std::make_pair(difference,&PFCand));}
      }

                }
                  }

              }}}}   
    
        //Sort these PtDifference vectors by ascending order in PtDifference from the Jet
          std::sort(PtDiffChHadr.begin(),PtDiffChHadr.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
          //std::sort(PtDiffChHadr.begin(),PtDiffChHadr.end(),pairCompare);
          std::sort(PtDiffNeutrHadr.begin(),PtDiffNeutrHadr.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;}); 
          std::sort(PtDiffPhotons.begin(),PtDiffPhotons.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
          //sort the tracks

          std::sort(PtDiffPosTrk.begin(),PtDiffPosTrk.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;}); 
          std::sort(PtDiffNegTrk.begin(),PtDiffNegTrk.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});

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
          //Charged Hadron begin
          nChHadrj1 = PtDiffChHadr.size();
          for (int i=0;i<nChHadrj1;i++){
            const reco::PFCandidate &ChHadr = *(PtDiffChHadr.at(i).second);
            //std::cout<<"ChHadr: "<<ChHadr.charge();
            ChHadrTotalCharge= ChHadrTotalCharge + ChHadr.charge();
            if(ChHadr.bestTrack()!=nullptr){//method points to null if there is no track 
              const reco::Track &Trk = *(ChHadr.bestTrack());
              recoTrkj1++;
              //std::cout<<"TrkCharge: "<<Trk.charge()<<std::endl;
              TrkTotalCharge= TrkTotalCharge + Trk.charge();
            }}
          if(nChHadrj1>0){
            const reco::PFCandidate &ChHadr1 = *(PtDiffChHadr.at(0).second);
            if(ChHadr1.bestTrack()!=nullptr){//lets store some info about Trk1 associated with highestPt Hadron
              const reco::Track &Trk = *(ChHadr1.bestTrack());
              Trk1Pt = Trk.pt();
              Trk1_Charge = Trk.charge();
              if(Trk.charge()==1){
                Trk1PosCharge = Trk.charge();
                nTrk1Pos = nTrk1Pos+1;
              }
              else{
                Trk1NegCharge = Trk.charge();
                nTrk1Neg = nTrk1Neg+1;
              }
              Trk1PtFrac = (Trk1Pt/jet.pt());
              Trk1Eta = Trk.eta();
              Trk1Phi = Trk.phi();
              }
            ChHadr1_pt = ChHadr1.pt();
            ChHadr1_Charge = ChHadr1.charge();
            ChHadr1_PtDiff = PtDiffChHadr.at(0).first;
            ChHadr1_ptfrac = (ChHadr1_pt/jet.pt());
            ChHadr1_eta = ChHadr1.eta();
            ChHadr1_phi = ChHadr1.phi();
            ChHadr1_vx = ChHadr1.vx();
            ChHadr1_vy = ChHadr1.vy();
            ChHadr1_vz = ChHadr1.vz();
          }//Close the loop for first Charged Hadron 
          //Start to look at second charged Hadron in the event  
          if(PtDiffChHadr.size()>=2){
            const reco::PFCandidate &ChHadr2 = *(PtDiffChHadr.at(1).second);
            if(ChHadr2.bestTrack()!=nullptr){//lets store some info about Trk2 associated with second highestPt charged Hadron
              const reco::Track &Trk = *(ChHadr2.bestTrack());
              Trk2Pt = Trk.pt();
              Trk2_Charge = Trk.charge();
              if(Trk.charge()==1){
                Trk2PosCharge = Trk.charge(); 
                 nTrk2Pos = nTrk2Pos+1;
              }
              else{
                Trk2NegCharge = Trk.charge();
                nTrk2Neg = nTrk2Neg+1;
              }
              Trk2PtFrac = (Trk2Pt/jet.pt());
              Trk2Eta = Trk.eta();
              Trk2Phi = Trk.phi();
              }
            ChHadr2_pt = ChHadr2.pt();
            ChHadr2_Charge = ChHadr2.charge();
            ChHadr2_PtDiff = PtDiffChHadr.at(1).first;
            ChHadr2_ptfrac = (ChHadr2_pt/jet.pt());
            ChHadr2_eta = ChHadr2.eta();
            ChHadr2_phi = ChHadr2.phi();
            ChHadr2_vx = ChHadr2.vx();
            ChHadr2_vy = ChHadr2.vy();
            ChHadr2_vz = ChHadr2.vz();
            } 
            dRTrk12 = reco::deltaR(Trk1Eta,Trk1Phi,Trk2Eta,Trk2Phi);//Reco ChargeHadron analysis ends
            
            dRChHadr12 = reco::deltaR(ChHadr1_eta,ChHadr1_phi,ChHadr2_eta,ChHadr2_phi);//Reco ChargeHadron analysis ends
            ChHadr12_Ptfrac = (ChHadr1_pt+ChHadr2_pt)/(j1Pt);
            ChHadr12_Charge = ChHadr1_Charge + ChHadr2_Charge;

            Trk12_Ptfrac = (Trk1Pt+Trk2Pt)/(j1Pt);
            Trk12_Charge = Trk1_Charge + Trk2_Charge;
        
            //Reco Neutral Hadrons(130) inside the jet that match to genlevel hadrons coming from Zprime. Expect ~0
            nNeutrHadrj1 = PtDiffNeutrHadr.size();
            
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
              //const reco::Candidate &Pho1_Ancestor = *(Pho1.mother(0));
              //std::cout<<"Pho1_Ancestor: "<<abs(Pho1_Ancestor.pdgId())<<std::endl;
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
              //const reco::Candidate &Pho2_Ancestor = *(Pho2.mother(0));
              //std::cout<<"Pho2_Ancestor: "<<abs(Pho2_Ancestor.pdgId())<<std::endl;
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

          } //closing the leading jet
        }//closing the ak4jets loop
           
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
