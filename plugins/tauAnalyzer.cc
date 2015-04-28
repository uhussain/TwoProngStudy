// -*- C++ -*-
//
// Package:    RecoTauTag/tauAnalyzer
// Class:      tauAnalyzer
// 
/**\class tauAnalyzer tauAnalyzer.cc RecoTauTag/tauAnalyzer/plugins/tauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Sat, 25 Apr 2015 15:51:58 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TTree.h"
#include "TauTrigMatch.h"

#include "FWCore/Framework/interface/EventSetup.h"
//
// class declaration
//

class tauAnalyzer : public edm::EDAnalyzer {
   public:
      explicit tauAnalyzer(const edm::ParameterSet&);
      ~tauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::InputTag tauSrc_;
      edm::InputTag discriminatorSrc_;
      TTree* tree;
      std::vector<Float_t>* pts_;
      std::vector<Int_t>* passDiscr_;
      double maxDR_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
tauAnalyzer::tauAnalyzer(const edm::ParameterSet& cfg)
{
   //now do what ever initialization is needed
  tauSrc_              = cfg.getParameter<edm::InputTag>("recoTau");
  discriminatorSrc_    = cfg.getParameter<edm::InputTag>("recoTauDiscriminator");

  edm::Service<TFileService> fs;
  //ntuple additions
  tree = fs->make<TTree>("Ntuple", "Ntuple");
  pts_ = new std::vector<Float_t>();
  passDiscr_ = new std::vector<Int_t>();
  tree->Branch("pt", "std::vector<float>", &pts_);
  tree->Branch("passDiscr", "std::vector<int>", &passDiscr_);
  maxDR_ = 0.5;
}


tauAnalyzer::~tauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete pts_;
   delete passDiscr_;
}


//
// member functions
//
// Get collection of generator particles with status 2

std::vector<const reco::GenParticle*> getGenParticleCollection(const edm::Event& evt) {
    std::vector<const reco::GenParticle*> output;
    edm::Handle< std::vector<reco::GenParticle> > handle;
    evt.getByLabel("genParticles", handle);
    // Loop over objects in current collection
    for (size_t j = 0; j < handle->size(); ++j) {
      const reco::GenParticle& object = handle->at(j);
      //if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
      if(object.status() == 2) output.push_back(&object);
    }
  return output;
}


// Get collection of pat::taus
//
std::vector<const pat::Tau*> getRecoCandCollections(const edm::Event& evt, const edm::InputTag& collection) {
    std::cout<<"Here 1.a"<<std::endl;
    std::vector<const pat::Tau*> output;
    std::cout<<"Here 1.b"<<std::endl;
    edm::Handle<std::vector<pat::Tau> > handle;
    evt.getByLabel(collection, handle);
    std::cout<<"Here 1.c"<<std::endl;
    // Loop over objects in current collection
    for (size_t j = 0; j < handle->size(); ++j) {
      const pat::Tau& object = handle->at(j);
      if(object.pt()>15. && fabs(object.eta())< 3.){
        output.push_back(&object);
      }
    }
  return output;
}
// Method to find the best match between tag tau and gen object. The best matched gen tau object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
const reco::GenParticle* findBestGenMatch(const pat::Tau* TagTauObj,
    std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
  const reco::GenParticle* output = NULL;
  double bestDeltaR = -1;
  for (size_t i = 0; i < GenPart.size(); ++i) {
    double deltaR = reco::deltaR(*TagTauObj, *GenPart[i]);
    if (deltaR < maxDR) {
      if (!output || deltaR < bestDeltaR) {
        output = GenPart[i];
        bestDeltaR = deltaR;
      }
    }
  }
  return output;
}

// ------------ method called for each event  ------------
void
tauAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //Handle<pat::TauCollection> tauObjects;
  //Handle<reco::PFTauCollection> tauObjects;
  //evt.getByLabel(tauSrc_, tauObjects);
  
  edm::Handle<reco::PFTauDiscriminator> discriminator;
  evt.getByLabel(discriminatorSrc_, discriminator);

/*  pts_->clear();
  passDiscr_->clear();
  for (unsigned int iTau = 0; iTau<taus->size() ; ++iTau){
        reco::PFTauRef tauCandidate(taus, iTau);
        pts_->push_back(tauCandidate->pt());
        // check if tau candidate has passed discriminator
        if( (*discriminator)[tauCandidate] > 0.5 ){
        // do something with your candidate
        passDiscr_->push_back(1);
        }
        else{
        passDiscr_->push_back(0);
        }
  }//end tau loop 
*/
  std::cout<<"Here 1"<<std::endl;
  std::vector<const pat::Tau*> tauObjects = getRecoCandCollections(evt, tauSrc_);
  
  std::cout<<"Here 2"<<std::endl;
  std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollection(evt);
  std::cout<<"Here 3"<<std::endl;
  for(unsigned int i = 0; i<tauObjects.size(); i++){
            const pat::Tau* TagTau = tauObjects[i];
            const reco::Candidate* bestGenMatch = findBestGenMatch(TagTau,GenObjects, maxDR_) ;
            std::cout<<"BestGenMatch Pt: "<<bestGenMatch->pt()<<std::endl;
   }

tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
tauAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
tauAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
tauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
tauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
tauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
tauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(tauAnalyzer);
