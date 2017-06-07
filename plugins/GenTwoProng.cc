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
#include "Genutil.h"

// class declaration
class GenTwoProng : public edm::EDAnalyzer {
    public:
        explicit GenTwoProng(const edm::ParameterSet&);
        ~GenTwoProng();

       	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters);
        const reco::GenParticle* findThirdHadron(std::vector<const reco::GenParticle*> hadronCands, reco::CandidatePtrVector signalCands, std::vector<const reco::GenParticle*> daughters, int *recoTrack, const reco::GenParticle* *genTrack3, double *trackDR);
        int GetDecayMode(std::vector<const reco::GenParticle*>& daughters);
        reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
        void eraseHadronCands(std::vector<const reco::GenParticle*>& hadronCands, reco::CandidatePtrVector signalCands);
        bool isNeutrino(const reco::Candidate* daughter);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);           
           
    private:
		    virtual void beginJob() override;
		    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		    virtual void endJob() override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken;
        edm::EDGetTokenT<std::vector<reco::GenJet> > jetsToken;
        edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken;
        edm::EDGetTokenT<std::vector <pat::PackedGenParticle> >packedGenToken;

        TTree* GenTree;
        int     run_;
        int  event_;
        int     lumis_;
        int nvtx;
        
        double j1Pt;
        double j2Pt;
        double j1Eta;
        double j1Phi;
        float j1ConsEtaPhiSpread;
        
        //Charged Hadron Variables
        int nChHadrj1;
        double ChHadrTotalCharge; 

        double ChHadr1_pt;
        double ChHadr1_Charge;
        double ChHadr1_PtDiff;
        double ChHadr1_ptfrac;
        double ChHadr1_eta;
        double ChHadr1_phi;
        double ChHadr1_vx;//returns the position of the point of closest approach to the PV 
        double ChHadr1_vy;
        double ChHadr1_vz;
        double ChHadr1_vertexRef;//reference to the PV itself


        double ChHadr12_Ptfrac;
        int ChHadr12_Charge;
        double dRChHadr12; //deltaR between two highest PtTracks
        
        double ChHadr2_pt;
        double ChHadr2_Charge;
        double ChHadr2_PtDiff;
        double ChHadr2_ptfrac;
        double ChHadr2_eta;
        double ChHadr2_phi;
        double ChHadr2_vx;//returns the x-coordinate of vertex position
        double ChHadr2_vy;
        double ChHadr2_vz;

        //Neutral Hadron Variables
        int nNeutrHadrj1;
        double NeutrHadrTotalCharge;

        double NeutrHadr1_pt;
        double NeutrHadr1_Charge;
        double NeutrHadr1_PtDiff;
        double NeutrHadr1_ptfrac;
        double NeutrHadr1_eta;
        double NeutrHadr1_phi;
        double NeutrHadr1_vx;//returns the position of the point of closest approach to the PV 
        double NeutrHadr1_vy;
        double NeutrHadr1_vz;
        
        double NeutrHadr2_pt;
        double NeutrHadr2_Charge;
        double NeutrHadr2_PtDiff;
        double NeutrHadr2_ptfrac;
        double NeutrHadr2_eta;
        double NeutrHadr2_phi;
        double NeutrHadr2_vx;//returns the x-coordinate of vertex position
        double NeutrHadr2_vy;
        double NeutrHadr2_vz;

        double NeutrHadr12_Ptfrac; 
        double dRNH12; //deltaR between two highest PtTracks

        //Photons in packedGenParticles
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

GenTwoProng::GenTwoProng(const edm::ParameterSet& iConfig):
    vtxToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    jetsToken(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    prunedGenToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken(consumes<std::vector<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
    edm::Service<TFileService> fs;

    GenTree = fs->make<TTree>("GenTree", "GenTree");
    GenTree->Branch("nvtx",&nvtx,"nvtx/I");

    GenTree->Branch("j1Pt",&j1Pt,"j1Pt/D");
    GenTree->Branch("j2Pt",&j2Pt,"j2Pt/D");
    GenTree->Branch("j1Eta",&j1Eta,"j1Eta/D");
    GenTree->Branch("j1Phi",&j1Phi,"j1Phi/D");
    GenTree->Branch("j1ConsEtaPhiSpread",&j1ConsEtaPhiSpread,"j1ConsEtaPhiSpread/F");
    
    GenTree->Branch("nChHadrj1",&nChHadrj1,"nChHadrj1/I");
    GenTree->Branch("ChHadrTotalCharge",&ChHadrTotalCharge,"ChHadrTotalCharge/D");
    
    
    GenTree->Branch("ChHadr1_pt",&ChHadr1_pt,"ChHadr1_pt/D");
    GenTree->Branch("ChHadr1_Charge",&ChHadr1_Charge,"ChHadr1_Charge/D");
    GenTree->Branch("ChHadr1_PtDiff",&ChHadr1_PtDiff,"ChHadr1_PtDiff/D");
    GenTree->Branch("ChHadr1_ptfrac",&ChHadr1_ptfrac,"ChHadr1_ptfrac/D");
    GenTree->Branch("ChHadr1_eta",&ChHadr1_eta,"ChHadr1_eta/D");
    GenTree->Branch("ChHadr1_phi",&ChHadr1_phi,"ChHadr1_phi/D");
    GenTree->Branch("ChHadr1_vx",&ChHadr1_vx,"ChHadr1_vx/D");
    GenTree->Branch("ChHadr1_vy",&ChHadr1_vy,"ChHadr1_vy/D");
    GenTree->Branch("ChHadr1_vz",&ChHadr1_vz,"ChHadr1_vz/D");
 
    GenTree->Branch("ChHadr2_pt",&ChHadr2_pt,"ChHadr2_pt/D");
    GenTree->Branch("ChHadr2_Charge",&ChHadr2_Charge,"ChHadr2_Charge/D");
    GenTree->Branch("ChHadr2_PtDiff",&ChHadr2_PtDiff,"ChHadr2_PtDiff/D");
    GenTree->Branch("ChHadr2_ptfrac",&ChHadr2_ptfrac,"ChHadr2_ptfrac/D");
    GenTree->Branch("ChHadr2_eta",&ChHadr2_eta,"ChHadr2_eta/D");
    GenTree->Branch("ChHadr2_phi",&ChHadr2_phi,"ChHadr2_phi/D");
    GenTree->Branch("ChHadr2_vx",&ChHadr2_vx,"ChHadr2_vx/D");
    GenTree->Branch("ChHadr2_vy",&ChHadr2_vy,"ChHadr2_vy/D");
    GenTree->Branch("ChHadr2_vz",&ChHadr2_vz,"ChHadr2_vz/D");

    GenTree->Branch("ChHadr12_Ptfrac",&ChHadr12_Ptfrac,"ChHadr12_Ptfrac/D"); 
    GenTree->Branch("ChHadr12_Charge",&ChHadr12_Charge,"ChHadr12_Charge/I");
    GenTree->Branch("dRChHadr12",&dRChHadr12,"dRChHadr12/D");

    //Neutral Hadrons
    GenTree->Branch("nNeutrHadrj1",&nNeutrHadrj1,"nNeutrHadrj1/I");
    GenTree->Branch("NeutrHadrTotalCharge",&NeutrHadrTotalCharge,"NeutrHadrTotalCharge/I");
     
    GenTree->Branch("NeutrHadr1_pt",&NeutrHadr1_pt,"NeutrHadr1_pt/D"); 
    GenTree->Branch("NeutrHadr1_Charge",&NeutrHadr1_Charge,"NeutrHadr1_Charge/D");
    GenTree->Branch("NeutrHadr1_PtDiff",&NeutrHadr1_PtDiff,"NeutrHadr1_PtDiff/D");
    GenTree->Branch("NeutrHadr1_ptfrac",&NeutrHadr1_ptfrac,"NeutrHadr1_ptfrac/D");
    GenTree->Branch("NeutrHadr1_eta",&NeutrHadr1_eta,"NeutrHadr1_eta/D");
    GenTree->Branch("NeutrHadr1_phi",&NeutrHadr1_phi,"NeutrHadr1_phi/D");
    GenTree->Branch("NeutrHadr1_vx",&NeutrHadr1_vx,"NeutrHadr1_vx/D");
    GenTree->Branch("NeutrHadr1_vy",&NeutrHadr1_vy,"NeutrHadr1_vy/D");
    GenTree->Branch("NeutrHadr1_vz",&NeutrHadr1_vz,"NeutrHadr1_vz/D");
 
    GenTree->Branch("NeutrHadr2_pt",&NeutrHadr2_pt,"NeutrHadr2_pt/D");
    GenTree->Branch("NeutrHadr2_Charge",&NeutrHadr2_Charge,"NeutrHadr2_Charge/D");
    GenTree->Branch("NeutrHadr2_PtDiff",&NeutrHadr2_PtDiff,"NeutrHadr2_PtDiff/D");
    GenTree->Branch("NeutrHadr2_ptfrac",&NeutrHadr2_ptfrac,"NeutrHadr2_ptfrac/D");
    GenTree->Branch("NeutrHadr2_eta",&NeutrHadr2_eta,"NeutrHadr2_eta/D");
    GenTree->Branch("NeutrHadr2_phi",&NeutrHadr2_phi,"NeutrHadr2_phi/D");
    GenTree->Branch("NeutrHadr2_vx",&NeutrHadr2_vx,"NeutrHadr2_vx/D");
    GenTree->Branch("NeutrHadr2_vy",&NeutrHadr2_vy,"NeutrHadr2_vy/D");
    GenTree->Branch("NeutrHadr2_vz",&NeutrHadr2_vz,"NeutrHadr2_vz/D");

    GenTree->Branch("NeutrHadr12_Ptfrac",&NeutrHadr12_Ptfrac,"NeutrHadr12_Ptfrac/D");
    GenTree->Branch("dRNH12",&dRNH12,"dRNH12/D");

    //Photons from pi0

    GenTree->Branch("nPhoj1",&nPhoj1,"nPhoj1/I");
    GenTree->Branch("PhoTotalCharge",&PhoTotalCharge,"PhoTotalCharge/I");

    GenTree->Branch("Pho1_pt",&Pho1_pt,"Pho1_pt/D"); 
    GenTree->Branch("Pho1_Charge",&Pho1_Charge,"Pho1_Charge/D");
    GenTree->Branch("Pho1_PtDiff",&Pho1_PtDiff,"Pho1_PtDiff/D");
    GenTree->Branch("Pho1_ptfrac",&Pho1_ptfrac,"Pho1_ptfrac/D");
    GenTree->Branch("Pho1_eta",&Pho1_eta,"Pho1_eta/D");
    GenTree->Branch("Pho1_phi",&Pho1_phi,"Pho1_phi/D");
    GenTree->Branch("Pho1_vx",&Pho1_vx,"Pho1_vx/D");
    GenTree->Branch("Pho1_vy",&Pho1_vy,"Pho1_vy/D");
    GenTree->Branch("Pho1_vz",&Pho1_vz,"Pho1_vz/D");
 
    GenTree->Branch("Pho2_pt",&Pho2_pt,"Pho2_pt/D");
    GenTree->Branch("Pho2_Charge",&Pho2_Charge,"Pho2_Charge/D");
    GenTree->Branch("Pho2_PtDiff",&Pho2_PtDiff,"Pho2_PtDiff/D");
    GenTree->Branch("Pho2_ptfrac",&Pho2_ptfrac,"Pho2_ptfrac/D");
    GenTree->Branch("Pho2_eta",&Pho2_eta,"Pho2_eta/D");
    GenTree->Branch("Pho2_phi",&Pho2_phi,"Pho2_phi/D");
    GenTree->Branch("Pho2_vx",&Pho2_vx,"Pho2_vx/D");
    GenTree->Branch("Pho2_vy",&Pho2_vy,"Pho2_vy/D");
    GenTree->Branch("Pho2_vz",&Pho2_vz,"Pho2_vz/D");


    GenTree->Branch("Pho12_Ptfrac",&Pho12_Ptfrac,"Pho12_Ptfrac/D");
    GenTree->Branch("dRPho12",&dRPho12,"dRPho12/D");
}

GenTwoProng::~GenTwoProng()
{
}

void GenTwoProng::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken, vertices);
    nvtx=vertices->size();
    //const reco::Vertex &PV = vertices->front();
    
    edm::Handle<std::vector<reco::GenJet> > ak4jets;
    iEvent.getByToken(jetsToken, ak4jets);

    
    edm::Handle<std::vector<reco::GenParticle> > prunedGenParticles;
    iEvent.getByToken(prunedGenToken, prunedGenParticles);
    
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticles;
    iEvent.getByToken(packedGenToken, packedGenParticles);
    
    
    ChHadrTotalCharge = 0;
    int nZprime = 0;
    //NeutralHadrons
    NeutrHadrTotalCharge = 0;

    PhoTotalCharge = 0;
    
    run_    = iEvent.id().run();
    event_  = iEvent.id().event();
    lumis_  = iEvent.luminosityBlock();
    
    std::cout<<"NewEvent"<<std::endl;
    std::cout<<run_<<":"<<event_<<":"<<lumis_<<std::endl;

    std::vector<const pat::PackedGenParticle*> packedHadrons;
    for(std::vector<pat::PackedGenParticle>::const_iterator genParticle = packedGenParticles->begin(); genParticle != packedGenParticles->end(); genParticle++){
        //std::cout<<"PackedGenParticles: "<<abs(genParticle->pdgId())<<std::endl;
        if(TMath::Abs(genParticle->pdgId()) == 211) packedHadrons.push_back(&(*genParticle));
    }
    //std::cout<<"PackedGenHadrons: "<<packedHadrons.size()<<std::endl;
    //Make list of possible PF charged pions to be used in tau reconstruction
    std::vector<const reco::GenParticle*> hadronCandidates;
    for(std::vector<reco::GenParticle>::const_iterator candidate= prunedGenParticles->begin(); candidate != prunedGenParticles->end(); candidate++){
        //std::cout<<"PrunedGenParticles: "<<abs(candidate->pdgId())<<std::endl;
        if(TMath::Abs(candidate->pdgId()) == 211) hadronCandidates.push_back(&(*candidate));
    }

    //std::cout<<"PrunedGenParticles: "<<prunedGenParticles->size()<<std::endl;
    for(uint32_t i=0; i < ak4jets->size(); i++){
        const reco::GenJet &jet = (*ak4jets)[i];
        if (i==0){
        j1Pt = jet.pt();
        j1Eta = jet.eta();
        j1Phi = jet.phi();
        j1ConsEtaPhiSpread = jet.constituentEtaPhiSpread();
        double difference = 0;
        //These vectors will store the PtDifference between a packed PFCandidate and a nearbyJet
        std::vector<std::pair<double,const pat::PackedGenParticle*>> PtDiffChHadr;//pdgId = abs(211) 
        std::vector<std::pair<double,const pat::PackedGenParticle*>> PtDiffNeutrHadr;//pdgId = abs(130)
        std::vector<std::pair<double,const pat::PackedGenParticle*>> PtDiffPhotons;//pdgId = abs(22)
        
        //Let's only save particles in these categories which originate from Z'
        for(uint32_t j=0;j<prunedGenParticles->size();j++){
          if(abs((*prunedGenParticles)[j].pdgId())==600001){
            const reco::Candidate * Zprime = &(*prunedGenParticles)[j];
            nZprime++;
            //Loop over all packedGenParticles
            if(nZprime<2){
            for(uint32_t j = 0; j < packedGenParticles->size(); j++) {
              const pat::PackedGenParticle &genCand = (*packedGenParticles)[j]; 
              //Find ancestor of this genCand
              const reco::Candidate * genCand_Ancestor = &(*(genCand.mother(0)));
              if(genCand_Ancestor != nullptr && isAncestor(Zprime,genCand_Ancestor)){
                //std::cout<<"genCand_Ancestor is Zprime: "<<std::endl;  
                if(reco::deltaR(jet.eta(),jet.phi(),genCand.eta(),genCand.phi())< 0.2){
                    if(abs(genCand.pdgId())== 211) {
                      difference = abs(jet.pt()-genCand.pt());
                      PtDiffChHadr.push_back({difference,&genCand});}

                    if(abs(genCand.pdgId())== 130) {
                      difference = abs(jet.pt()-genCand.pt());
                      PtDiffNeutrHadr.push_back(std::make_pair(difference,&genCand));}

                    if(abs(genCand.pdgId())== 22) {
                      difference = abs(jet.pt()-genCand.pt());
                      PtDiffPhotons.push_back(std::make_pair(difference,&genCand));}
                }
      }
    }
            }//nZprime Condition is to avoid duplication
          }}
        //Sort these PtDifference vectors by ascending order in PtDifference from the Jet
        std::sort(PtDiffChHadr.begin(),PtDiffChHadr.end(),[](const auto& p1, const auto& p2){return p1.first<p2.first;});
        //std::sort(PtDiffChHadr.begin(),PtDiffChHadr.end(),pairCompare);
        std::sort(PtDiffNeutrHadr.begin(),PtDiffNeutrHadr.end(),pairCompare); //or you can use the specific pairCompare function. Sort has some issues with template that I cannot fix right now?
        std::sort(PtDiffPhotons.begin(),PtDiffPhotons.end(),pairCompare);//pairCompare is defined in util.h
          //Charged Hadron begin
          nChHadrj1 = PtDiffChHadr.size();
          for (int i=0;i<nChHadrj1;i++){
            const pat::PackedGenParticle &ChHadr = *(PtDiffChHadr.at(i).second);
            //std::cout<<"ChHadr: "<<ChHadr.charge();
            ChHadrTotalCharge= ChHadrTotalCharge + ChHadr.charge();
            //const reco::Candidate &Ancestor = *(ChHadr.mother(0));
            //std::cout<<"AncestorID: "<<abs(Ancestor.pdgId())<<std::endl;
            }
          //std::cout<<"TotalChHadrons in dR<0.1: "<<nChHadrj1<<std::endl;
          //std::cout<<"ChHadrTotalCharge: "<<ChHadrTotalCharge<<std::endl;
          if(nChHadrj1>0){
            const pat::PackedGenParticle &ChHadr1 = *(PtDiffChHadr.at(0).second);
            ChHadr1_pt = ChHadr1.pt(); //Does this Pt match with Trk1Pt at recoLevel
            ChHadr1_Charge = ChHadr1.charge();
            //std::cout<<"ChHadr1_Charge: "<<ChHadr1_Charge<<std::endl;
            ChHadr1_PtDiff = PtDiffChHadr.at(0).first;
            ChHadr1_ptfrac = (ChHadr1_pt/jet.pt());
            ChHadr1_eta = ChHadr1.eta();
            ChHadr1_phi = ChHadr1.phi();
            ChHadr1_vx = ChHadr1.vx();
            ChHadr1_vy = ChHadr1.vy();
            ChHadr1_vz = ChHadr1.vz();
          }//Close the loop for first Charged Hadron
            //std::cout<<"nZprime: "<<nZprime<<std::endl;
          //Start to look at second charged Hadron in the event  
          if(PtDiffChHadr.size()>=2){
            const pat::PackedGenParticle &ChHadr2 = *(PtDiffChHadr.at(1).second);
            ChHadr2_pt = ChHadr2.pt();
            ChHadr2_Charge = ChHadr2.charge();
            //const reco::Candidate &ChHadr2_Ancestor = *(ChHadr2.mother(2));
            //std::cout<<"ChHadr2_Ancestor: "<<abs(ChHadr2_Ancestor.pdgId())<<std::endl;
            ChHadr2_PtDiff = PtDiffChHadr.at(1).first;
            ChHadr2_ptfrac = (ChHadr2_pt/j1Pt);
            ChHadr2_eta = ChHadr2.eta();
            ChHadr2_phi = ChHadr2.phi();
            ChHadr2_vx = ChHadr2.vx();
            ChHadr2_vy = ChHadr2.vy();
            ChHadr2_vz = ChHadr2.vz();
            } 
            dRChHadr12 = reco::deltaR(ChHadr1_eta,ChHadr1_phi,ChHadr2_eta,ChHadr2_phi);//Reco ChargeHadron analysis ends
            ChHadr12_Ptfrac = (ChHadr1_pt+ChHadr2_pt)/(j1Pt);
            ChHadr12_Charge = ChHadr1_Charge + ChHadr2_Charge;
          //Reco Neutral Hadron Analysis begins
          nNeutrHadrj1 = PtDiffNeutrHadr.size();
          for (int i=0;i<nNeutrHadrj1;i++){
            const pat::PackedGenParticle &NeutrHadr = *(PtDiffNeutrHadr.at(i).second);
            NeutrHadrTotalCharge = NeutrHadrTotalCharge + NeutrHadr.charge();
          }
          if(nNeutrHadrj1>0){
            const pat::PackedGenParticle &NeutrHadr1 = *(PtDiffNeutrHadr.at(0).second);
            NeutrHadr1_pt = NeutrHadr1.pt();
            NeutrHadr1_Charge = NeutrHadr1.charge();
            NeutrHadr1_PtDiff = PtDiffNeutrHadr.at(0).first;
            NeutrHadr1_ptfrac = (NeutrHadr1_pt/jet.pt());
            NeutrHadr1_eta = NeutrHadr1.eta();
            NeutrHadr1_phi = NeutrHadr1.phi();
            NeutrHadr1_vx = NeutrHadr1.vx();
            NeutrHadr1_vy = NeutrHadr1.vy();
            NeutrHadr1_vz = NeutrHadr1.vz();
          }
            //std::cout<<"Where is this error?"<<std::endl;
            if(PtDiffNeutrHadr.size()>=2){
            const pat::PackedGenParticle &NeutrHadr2 = *(PtDiffNeutrHadr.at(1).second);
            NeutrHadr2_pt = NeutrHadr2.pt();
            NeutrHadr2_Charge = NeutrHadr2.charge();
            NeutrHadr2_PtDiff = PtDiffNeutrHadr.at(1).first;
            NeutrHadr2_ptfrac = (NeutrHadr2_pt/jet.pt());
            NeutrHadr2_eta = NeutrHadr2.eta();
            NeutrHadr2_phi = NeutrHadr2.phi();
            NeutrHadr2_vx = NeutrHadr2.vx();
            NeutrHadr2_vy = NeutrHadr2.vy();
            NeutrHadr2_vz = NeutrHadr2.vz();
            } 
            dRNH12 = reco::deltaR(NeutrHadr1_eta,NeutrHadr1_phi,NeutrHadr2_eta,NeutrHadr2_phi); //Reco Neutral Hadron Analysis ends 
            NeutrHadr12_Ptfrac = (NeutrHadr1_pt+NeutrHadr2_pt)/(j1Pt);

            //Add genlevel photon information, especially photons that come from pi0 and Z' in prunedCollection
            
          nPhoj1 = PtDiffPhotons.size();
          for (int i=0;i<nPhoj1;i++){
            const pat::PackedGenParticle &Pho = *(PtDiffPhotons.at(i).second);
            PhoTotalCharge = PhoTotalCharge + Pho.charge();
          }
          if(nPhoj1>0){
            const pat::PackedGenParticle &Pho1 = *(PtDiffPhotons.at(0).second);
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
            const pat::PackedGenParticle &Pho2 = *(PtDiffPhotons.at(1).second);
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
        if(i==1){
          j2Pt = jet.pt();
        }
        }//closing the ak4jets loop
           
            GenTree->Fill();
}
bool GenTwoProng::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
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
GenTwoProng::beginJob()
{
}


// ------------ method called once each job just after ending the event loop  ------------
	void 
GenTwoProng::endJob() 
{
}

//void GenTwoProng::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
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
//int GenTwoProng::GetDecayMode(std::vector<const reco::GenParticle*>& daughters){
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
//bool GenTwoProng::isNeutrino(const reco::Candidate* daughter)
//{
//  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
//}
////Gets visible 4-momentum of a particle from list of daughters
//reco::Candidate::LorentzVector GenTwoProng::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
//	reco::Candidate::LorentzVector p4_vis(0,0,0,0);	
// 	for(size_t i = 0; i < daughters.size(); ++i){
// 		if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
// 			p4_vis += daughters[i]->p4();
//		}
// 	}
//	return p4_vis;
//}
//
//void GenTwoProng::eraseHadronCands(std::vector<const reco::GenParticle*>& hadronCands, reco::CandidatePtrVector signalCands){
//    for(size_t i=0; i<signalCands.size(); ++i){
//        for(size_t j=0; j<hadronCands.size(); ++j){
//            if(signalCands[i]->pt()==hadronCands[j]->pt()) hadronCands.erase(hadronCands.begin()+j);
//        }
//    }
//}
////Finds closest charged pion not already used in 2-prong tau reconstruction
//const reco::GenParticle* GenTwoProng::findThirdHadron(std::vector<const reco::GenParticle*> hadronCands, reco::CandidatePtrVector signalCands, std::vector<const reco::GenParticle*> daughters,int *recoTrack, const reco::GenParticle* *genTrack3, double *trackDR){
//    const reco::GenParticle* thirdHadron = hadronCands[0];
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
GenTwoProng::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	  edm::ParameterSetDescription desc;
	  desc.setUnknown();
	  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenTwoProng);
