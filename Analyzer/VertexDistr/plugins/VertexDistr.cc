// -*- C++ -*-
//
// Package:    Analyzer/VertexDistr
// Class:      VertexDistr
//
/**\class VertexDistr VertexDistr.cc Analyzer/VertexDistr/plugins/VertexDistr.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nicola Minafra
//         Created:  Tue, 16 Oct 2018 13:25:00 GMT
//
//


// system include files
#include <memory>

#include "TProfile2D.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
// Particle flow collection
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class VertexDistr : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit VertexDistr(const edm::ParameterSet&);
      ~VertexDistr();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      edm::EDGetTokenT<std::vector< reco::PFCandidate > > pfToken_;

      TFileDirectory histo_dir_;

      TProfile2D* vxProfile_;
      TProfile2D* vyProfile_;
      TProfile2D* vzProfile_;
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
VertexDistr::VertexDistr(const edm::ParameterSet& iConfig) :
 pfToken_            ( consumes<std::vector< reco::PFCandidate> >     ( iConfig.getParameter<edm::InputTag>( "pfTag" ) ) )

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   histo_dir_ = fs->mkdir("histograms");

   vxProfile_ = histo_dir_.make<TProfile2D>("vxProfile", "vxProfile vs BX vs LS; BX; LS", 3600, 0, 3600, 1000, 0, 1000, -20, 20);
   vyProfile_ = histo_dir_.make<TProfile2D>("vyProfile", "vyProfile vs BX vs LS; BX; LS", 3600, 0, 3600, 1000, 0, 1000, -20, 20);
   vzProfile_ = histo_dir_.make<TProfile2D>("vzProfile", "vzProfile vs BX vs LS; BX; LS", 3600, 0, 3600, 1000, 0, 1000, -20, 20);

}


VertexDistr::~VertexDistr()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
VertexDistr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle <std::vector<reco::PFCandidate>> PFCandidates;
  iEvent.getByToken(pfToken_,PFCandidates);

  math::XYZTLorentzVector allCands(0.,0.,0.,0.);
  for( unsigned int i = 0; i < PFCandidates->size(); ++i ) {
    const reco::PFCandidate* pfAll = &((*PFCandidates)[i]);
    // std::cout<< pfAll->particleId() << std::endl;
    // std::cout<< pfAll->vz() << std::endl;
    vxProfile_->Fill(iEvent.bunchCrossing(), iEvent.luminosityBlock(), pfAll->vx());
    vyProfile_->Fill(iEvent.bunchCrossing(), iEvent.luminosityBlock(), pfAll->vy());
    vzProfile_->Fill(iEvent.bunchCrossing(), iEvent.luminosityBlock(), pfAll->vz());
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
VertexDistr::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
VertexDistr::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexDistr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexDistr);
