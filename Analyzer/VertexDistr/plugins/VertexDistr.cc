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
#include "TProfile.h"
#include "TGraph.h"
#include "TH2F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// LHCInfo
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"
#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
// Particle flow collection
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
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
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > barrelEcalHits_token_;
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > endcapEcalHits_token_;
      edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > ppsTimingHits_token_;

      std::vector<float>  ave;
      std::vector<float>  dif;
      std::vector<float>  vc1;
      std::vector<float>  vc2;
      std::vector<float>  aveRF;
      std::vector<float>  difRF;
      std::vector<float>  vc1RF;
      std::vector<float>  vc2RF;

      std::vector<unsigned int> train_zero;
      std::vector<unsigned int> train_notzero;
      std::vector<unsigned int> long_train_notzero;
      std::vector<unsigned int> subtrain_num;
      std::vector<unsigned int> train_num;

      unsigned int lastLumiSection;


      TFileDirectory histo_dir_;

      TProfile2D* vxProfile_;
      TProfile2D* vyProfile_;
      TProfile2D* vzProfile_;
      TProfile2D* timeProfile_;

      TProfile2D* vxVsTime_;
      TProfile2D* vyVsTime_;
      TProfile2D* vzVsTime_;
      unsigned int t0;


      TProfile* aveVsPFx_;
      TProfile* aveVsPFy_;
      TProfile* aveVsPFz_;
      TProfile* difVsPFx_;
      TProfile* difVsPFy_;
      TProfile* difVsPFz_;
      TProfile* vc1VsPFx_;
      TProfile* vc1VsPFy_;
      TProfile* vc1VsPFz_;
      TProfile* vc2VsPFx_;
      TProfile* vc2VsPFy_;
      TProfile* vc2VsPFz_;

      TProfile* aveRFVsPFx_;
      TProfile* aveRFVsPFy_;
      TProfile* aveRFVsPFz_;
      TProfile* difRFVsPFx_;
      TProfile* difRFVsPFy_;
      TProfile* difRFVsPFz_;
      TProfile* vc1RFVsPFx_;
      TProfile* vc1RFVsPFy_;
      TProfile* vc1RFVsPFz_;
      TProfile* vc2RFVsPFx_;
      TProfile* vc2RFVsPFy_;
      TProfile* vc2RFVsPFz_;

      TH2F* avgVsPFz_2D_;
      TH2F* difVsPFz_2D_;

      TProfile* ppsDiamond45VsBx_;
      TProfile* ppsDiamond56VsBx_;
      TProfile* ppsDiamond45VsBx_low_;
      TProfile* ppsDiamond56VsBx_low_;
      TProfile* ppsDiamond45VsBx_corrected_;
      TProfile* ppsDiamond56VsBx_corrected_;
      TProfile* ppsDiamond45VsBx_corrected_low_;
      TProfile* ppsDiamond56VsBx_corrected_low_;
      TProfile* ppsDiamond45VsBx_correctedRF_;
      TProfile* ppsDiamond56VsBx_correctedRF_;
      TProfile* ppsDiamond45VsBx_correctedRF_low_;
      TProfile* ppsDiamond56VsBx_correctedRF_low_;

      TProfile* aveVsBx_;
      TProfile* difVsBx_;
      TProfile* vc1VsBx_;
      TProfile* vc2VsBx_;

      TProfile* aveRFVsBx_;
      TProfile* difRFVsBx_;
      TProfile* vc1RFVsBx_;
      TProfile* vc2RFVsBx_;

      TGraph aveVsaveRF_;
      TGraph difVsdifRF_;
      TGraph vc1Vsvc1RF_;
      TGraph vc2Vsvc2RF_;
      unsigned int graph_ctr_1_;


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
 pfToken_            ( consumes<std::vector< reco::PFCandidate> >     ( iConfig.getParameter<edm::InputTag>( "pfTag" ) ) ),
 barrelEcalHits_token_ ( consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("barrelEcalHits")) ),
 endcapEcalHits_token_ ( consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("endcapEcalHits")) ),
 ppsTimingHits_token_  ( consumes< edm::DetSetVector<CTPPSDiamondRecHit> >    ( iConfig.getParameter<edm::InputTag>( "tagDiamondRecHits" ) ) ),
 t0(0), graph_ctr_1_(0)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   histo_dir_ = fs->mkdir("histograms");

   vxProfile_ = histo_dir_.make<TProfile2D>("vxProfile", "vxProfile vs BX vs LS; BX; LS", 3600, 0, 3600, 1000, 0, 1000, -20, 20);
   vyProfile_ = histo_dir_.make<TProfile2D>("vyProfile", "vyProfile vs BX vs LS; BX; LS", 3600, 0, 3600, 1000, 0, 1000, -20, 20);
   vzProfile_ = histo_dir_.make<TProfile2D>("vzProfile", "vzProfile vs BX vs LS; BX; LS", 3600, 0, 3600, 1000, 0, 1000, -20, 20);
   timeProfile_ = histo_dir_.make<TProfile2D>("timeProfile", "timeProfile vs BX vs LS; BX; LS", 3600, 0, 3600, 1000, 0, 1000, -20, 20);

   vxVsTime_ = histo_dir_.make<TProfile2D>("vxVsTime", "vxVsTime vs BX vs time; BX; s", 1000, 0, 1000, 1000, 0, 1000, -20, 20);
   vyVsTime_ = histo_dir_.make<TProfile2D>("vyVsTime", "vyVsTime vs BX vs time; BX; s", 1000, 0, 1000, 1000, 0, 1000, -20, 20);
   vzVsTime_ = histo_dir_.make<TProfile2D>("vzVsTime", "vzVsTime vs BX vs time; BX; s", 1000, 0, 1000, 1000, 0, 1000, -20, 20);

   aveVsPFx_ = histo_dir_.make<TProfile>("aveVsPFx", "aveVsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   aveVsPFy_ = histo_dir_.make<TProfile>("aveVsPFy", "aveVsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   aveVsPFz_ = histo_dir_.make<TProfile>("aveVsPFz", "aveVsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);
   difVsPFx_ = histo_dir_.make<TProfile>("difVsPFx", "difVsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   difVsPFy_ = histo_dir_.make<TProfile>("difVsPFy", "difVsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   difVsPFz_ = histo_dir_.make<TProfile>("difVsPFz", "difVsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc1VsPFx_ = histo_dir_.make<TProfile>("vc1VsPFx", "vc1VsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc1VsPFy_ = histo_dir_.make<TProfile>("vc1VsPFy", "vc1VsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc1VsPFz_ = histo_dir_.make<TProfile>("vc1VsPFz", "vc1VsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc2VsPFx_ = histo_dir_.make<TProfile>("vc2VsPFx", "vc2VsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc2VsPFy_ = histo_dir_.make<TProfile>("vc2VsPFy", "vc2VsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc2VsPFz_ = histo_dir_.make<TProfile>("vc2VsPFz", "vc2VsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);

   avgVsPFz_2D_ = histo_dir_.make<TH2F>("avgVsPFz_2D", "avgVsPFz_2D; ns; cm", 100, -0.05, 0.05, 30, -1.5, 1.5);
   difVsPFz_2D_ = histo_dir_.make<TH2F>("difVsPFz_2D", "difVsPFz_2D; ns; cm", 100, -0.05, 0.05, 30, -1.5, 1.5);

   aveRFVsPFx_ = histo_dir_.make<TProfile>("aveRFVsPFx", "aveRFVsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   aveRFVsPFy_ = histo_dir_.make<TProfile>("aveRFVsPFy", "aveRFVsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   aveRFVsPFz_ = histo_dir_.make<TProfile>("aveRFVsPFz", "aveRFVsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);
   difRFVsPFx_ = histo_dir_.make<TProfile>("difRFVsPFx", "difRFVsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   difRFVsPFy_ = histo_dir_.make<TProfile>("difRFVsPFy", "difRFVsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   difRFVsPFz_ = histo_dir_.make<TProfile>("difRFVsPFz", "difRFVsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc1RFVsPFx_ = histo_dir_.make<TProfile>("vc1RFVsPFx", "vc1RFVsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc1RFVsPFy_ = histo_dir_.make<TProfile>("vc1RFVsPFy", "vc1RFVsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc1RFVsPFz_ = histo_dir_.make<TProfile>("vc1RFVsPFz", "vc1RFVsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc2RFVsPFx_ = histo_dir_.make<TProfile>("vc2RFVsPFx", "vc2RFVsPFx; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc2RFVsPFy_ = histo_dir_.make<TProfile>("vc2RFVsPFy", "vc2RFVsPFy; ns; cm", 5000, -0.5, 0.5, -20, 20);
   vc2RFVsPFz_ = histo_dir_.make<TProfile>("vc2RFVsPFz", "vc2RFVsPFz; ns; cm", 5000, -0.5, 0.5, -20, 20);

   ppsDiamond45VsBx_ = histo_dir_.make<TProfile>("ppsDiamond45VsBx", "ppsDiamond45VsBx; BX; ns", 3600, 0, 3600, -20, 20);
   ppsDiamond56VsBx_ = histo_dir_.make<TProfile>("ppsDiamond56VsBx", "ppsDiamond56VsBx; BX; ns", 3600, 0, 3600, -20, 20);
   ppsDiamond45VsBx_corrected_ = histo_dir_.make<TProfile>("ppsDiamond45VsBx_corrected", "ppsDiamond45VsBx_corrected; BX; ns", 3600, 0, 3600, -20, 20);
   ppsDiamond56VsBx_corrected_ = histo_dir_.make<TProfile>("ppsDiamond56VsBx_corrected", "ppsDiamond56VsBx_corrected; BX; ns", 3600, 0, 3600, -20, 20);
   ppsDiamond45VsBx_correctedRF_ = histo_dir_.make<TProfile>("ppsDiamond45VsBx_correctedRF", "ppsDiamond45VsBx_correctedRF; BX; ns", 3600, 0, 3600, -20, 20);
   ppsDiamond56VsBx_correctedRF_ = histo_dir_.make<TProfile>("ppsDiamond56VsBx_correctedRF", "ppsDiamond56VsBx_correctedRF; BX; ns", 3600, 0, 3600, -20, 20);
   ppsDiamond45VsBx_low_ = histo_dir_.make<TProfile>("ppsDiamond45VsBx_low", "ppsDiamond45VsBx_low; BX; ns", 360, 0, 3600, -20, 20);
   ppsDiamond56VsBx_low_ = histo_dir_.make<TProfile>("ppsDiamond56VsBx_low", "ppsDiamond56VsBx_low; BX; ns", 360, 0, 3600, -20, 20);
   ppsDiamond45VsBx_corrected_low_ = histo_dir_.make<TProfile>("ppsDiamond45VsBx_corrected_low", "ppsDiamond45VsBx_corrected_low; BX; ns", 360, 0, 3600, -20, 20);
   ppsDiamond56VsBx_corrected_low_ = histo_dir_.make<TProfile>("ppsDiamond56VsBx_corrected_low", "ppsDiamond56VsBx_corrected_low; BX; ns", 360, 0, 3600, -20, 20);
   ppsDiamond45VsBx_correctedRF_low_ = histo_dir_.make<TProfile>("ppsDiamond45VsBx_correctedRF_low", "ppsDiamond45VsBx_correctedRF_low; BX; ns", 360, 0, 3600, -20, 20);
   ppsDiamond56VsBx_correctedRF_low_ = histo_dir_.make<TProfile>("ppsDiamond56VsBx_correctedRF_low", "ppsDiamond56VsBx_correctedRF_low; BX; ns", 360, 0, 3600, -20, 20);

   aveVsBx_ = histo_dir_.make<TProfile>("aveVsBx", "aveVsBx; BX; ns", 3600, 0, 3600, -2, 2);
   difVsBx_ = histo_dir_.make<TProfile>("difVsBx", "difVsBx; BX; ns", 3600, 0, 3600, -2, 2);
   vc1VsBx_ = histo_dir_.make<TProfile>("vc1VsBx", "vc1VsBx; BX; ns", 3600, 0, 3600, -2, 2);
   vc2VsBx_ = histo_dir_.make<TProfile>("vc2VsBx", "vc2VsBx; BX; ns", 3600, 0, 3600, -2, 2);

   aveRFVsBx_ = histo_dir_.make<TProfile>("aveRFVsBx", "aveRFVsBx; BX; ns", 3600, 0, 3600, -2, 2);
   difRFVsBx_ = histo_dir_.make<TProfile>("difRFVsBx", "difRFVsBx; BX; ns", 3600, 0, 3600, -2, 2);
   vc1RFVsBx_ = histo_dir_.make<TProfile>("vc1RFVsBx", "vc1RFVsBx; BX; ns", 3600, 0, 3600, -2, 2);
   vc2RFVsBx_ = histo_dir_.make<TProfile>("vc2RFVsBx", "vc2RFVsBx; BX; ns", 3600, 0, 3600, -2, 2);

   aveVsaveRF_.SetName("aveVsaveRF");
   aveVsaveRF_.SetTitle("aveVsaveRF");
   difVsdifRF_.SetName("difVsdifRF");
   difVsdifRF_.SetTitle("difVsdifRF");
   vc1Vsvc1RF_.SetName("vc1Vsvc1RF");
   vc1Vsvc1RF_.SetTitle("vc1Vsvc1RF");
   vc2Vsvc2RF_.SetName("vc2Vsvc2RF");
   vc2Vsvc2RF_.SetTitle("vc2Vsvc2RF");

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

  edm::Handle<EcalRecHitCollection> barrelHitHandle;
  EcalRecHitCollection barrelRecHits;
  iEvent.getByToken(barrelEcalHits_token_, barrelHitHandle);
  edm::Handle<EcalRecHitCollection> endcapHitHandle;
  iEvent.getByToken(endcapEcalHits_token_, endcapHitHandle);
  EcalRecHitCollection endcapRecHits;

  edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > diamondRecHits;
  iEvent.getByToken( ppsTimingHits_token_, diamondRecHits );


  if (iEvent.bunchCrossing()==0) std::cout<<"FOUND BX==0"<<std::endl;
  unsigned int bx = iEvent.bunchCrossing() - 1;
  unsigned int lumiBlock = iEvent.luminosityBlock();
  if (lumiBlock<250) return;

  edm::ESHandle<LHCInfo> lhcInfoHnd;
  iSetup.get<LHCInfoRcd>().get(lhcInfoHnd);

  if( lastLumiSection != lumiBlock ){
    lastLumiSection = lumiBlock;

    const LHCInfo* lhcInfo = lhcInfoHnd.product();

    bool first_zero( true );
    bool first_notzero( true );
    unsigned int zero(0);
    unsigned int notzero(0);
    unsigned int longnotzero(0);
    unsigned int count(0);
    unsigned int longcount(0);

    unsigned int BXInfoSize = lhcInfo->beam1VC().size();
    if (BXInfoSize != ave.size()){
      ave.resize(BXInfoSize);
      dif.resize(BXInfoSize);
      vc1.resize(BXInfoSize);
      vc2.resize(BXInfoSize);
      aveRF.resize(BXInfoSize);
      difRF.resize(BXInfoSize);
      vc1RF.resize(BXInfoSize);
      vc2RF.resize(BXInfoSize);

      train_zero.resize(BXInfoSize);
      train_notzero.resize(BXInfoSize);
      long_train_notzero.resize(BXInfoSize);
      subtrain_num.resize(BXInfoSize);
      train_num.resize(BXInfoSize);
    }

    for( unsigned int i  =  0; i < BXInfoSize; i++ ){

      ave.at(i) = ((lhcInfo->beam1VC()[i]+lhcInfo->beam2VC()[i]))*(2.5/360.0);
      dif.at(i) = (lhcInfo->beam2VC()[i]-lhcInfo->beam1VC()[i])*(2.5/360.0);
      vc1.at(i) = (lhcInfo->beam1VC()[i])*(2.5/360.0);
      vc2.at(i) = (lhcInfo->beam2VC()[i])*(2.5/360.0);


      aveRF.at(i) = ((lhcInfo->beam1RF()[i]+lhcInfo->beam2RF()[i]))*(2.5/360.0);
      difRF.at(i) = (lhcInfo->beam2RF()[i]-lhcInfo->beam1RF()[i])*(2.5/360.0);
      vc1RF.at(i) = (lhcInfo->beam1RF()[i])*(2.5/360.0);
      vc2RF.at(i) = (lhcInfo->beam2RF()[i])*(2.5/360.0);



      if( lhcInfo->beam1VC()[i] == 0.0 ) {
        if( first_zero == true ) {
                first_zero = false;
                first_notzero = true;
                zero = 0;
        }
        zero++;
      }
      else {
        if( first_notzero == true ) {
                first_notzero = false;
                first_zero = true;
                notzero = 0;
                if( zero > 10 ){
                        longnotzero = 0;
                        longcount++;
                }
                count++;
        }
        notzero++;
        longnotzero++;
      }
      train_zero.at(i) = zero;
      train_notzero.at(i) = notzero;
      long_train_notzero.at(i) = longnotzero;
      subtrain_num.at(i) = count;
      train_num.at(i) = longcount;
    }
  }

  //subtrain_position = train_notzero[bx];
  //train_position = long_train_notzero[bx];
  //subtrain_number = subtrain_num[bx];
  //train_number = train_num[bx];

  float ave_phase = ave[bx];
  float dif_phase = dif[bx];
  float vc1_phase = vc1[bx];
  float vc2_phase = vc2[bx];

  aveVsBx_->Fill(iEvent.bunchCrossing(), ave_phase);
  difVsBx_->Fill(iEvent.bunchCrossing(), dif_phase);
  vc1VsBx_->Fill(iEvent.bunchCrossing(), vc1_phase);
  vc2VsBx_->Fill(iEvent.bunchCrossing(), vc2_phase);


  float aveRF_phase = aveRF[bx];
  float difRF_phase = difRF[bx];
  float vc1RF_phase = vc1RF[bx];
  float vc2RF_phase = vc2RF[bx];

  aveRFVsBx_->Fill(iEvent.bunchCrossing(), aveRF_phase);
  difRFVsBx_->Fill(iEvent.bunchCrossing(), difRF_phase);
  vc1RFVsBx_->Fill(iEvent.bunchCrossing(), vc1RF_phase);
  vc2RFVsBx_->Fill(iEvent.bunchCrossing(), vc2RF_phase);

  aveVsaveRF_.SetPoint(graph_ctr_1_, ave_phase, aveRF_phase);
  difVsdifRF_.SetPoint(graph_ctr_1_, dif_phase, difRF_phase);
  vc1Vsvc1RF_.SetPoint(graph_ctr_1_, vc1_phase, vc1RF_phase);
  vc2Vsvc2RF_.SetPoint(graph_ctr_1_, vc2_phase, vc2RF_phase);
  ++graph_ctr_1_;



  if (t0 == 0) t0 = iEvent.time().unixTime();
  // std::cout<< "Starting for time: " << t0 << std::endl;

  math::XYZTLorentzVector allCands(0.,0.,0.,0.);
  for( unsigned int i = 0; i < PFCandidates->size(); ++i ) {
    const reco::PFCandidate* pfAll = &((*PFCandidates)[i]);
    if (pfAll->particleId() != 3) continue; // only muons


    // std::cout<< pfAll->particleId() << std::endl;
    // std::cout<< pfAll->vz() << std::endl;
    vxProfile_->Fill(iEvent.bunchCrossing(), iEvent.luminosityBlock(), pfAll->vx());
    vyProfile_->Fill(iEvent.bunchCrossing(), iEvent.luminosityBlock(), pfAll->vy());
    vzProfile_->Fill(iEvent.bunchCrossing(), iEvent.luminosityBlock(), pfAll->vz());

    vxVsTime_->Fill(iEvent.bunchCrossing(), iEvent.time().unixTime() - t0, pfAll->vx());
    vyVsTime_->Fill(iEvent.bunchCrossing(), iEvent.time().unixTime() - t0, pfAll->vy());
    vzVsTime_->Fill(iEvent.bunchCrossing(), iEvent.time().unixTime() - t0, pfAll->vz());

    aveVsPFx_->Fill(1*ave_phase, pfAll->vx());
    aveVsPFy_->Fill(1*ave_phase, pfAll->vy());
    aveVsPFz_->Fill(1*ave_phase, pfAll->vz());
    difVsPFx_->Fill(1*dif_phase, pfAll->vx());
    difVsPFy_->Fill(1*dif_phase, pfAll->vy());
    difVsPFz_->Fill(1*dif_phase, pfAll->vz());

    avgVsPFz_2D_->Fill(1*ave_phase, pfAll->vz());
    difVsPFz_2D_->Fill(1*dif_phase, pfAll->vz());

    vc1VsPFx_->Fill(1*vc1_phase, pfAll->vx());
    vc1VsPFy_->Fill(1*vc1_phase, pfAll->vy());
    vc1VsPFz_->Fill(1*vc1_phase, pfAll->vz());
    vc2VsPFx_->Fill(1*vc2_phase, pfAll->vx());
    vc2VsPFy_->Fill(1*vc2_phase, pfAll->vy());
    vc2VsPFz_->Fill(1*vc2_phase, pfAll->vz());

    aveRFVsPFx_->Fill(1*aveRF_phase, pfAll->vx());
    aveRFVsPFy_->Fill(1*aveRF_phase, pfAll->vy());
    aveRFVsPFz_->Fill(1*aveRF_phase, pfAll->vz());
    difRFVsPFx_->Fill(1*difRF_phase, pfAll->vx());
    difRFVsPFy_->Fill(1*difRF_phase, pfAll->vy());
    difRFVsPFz_->Fill(1*difRF_phase, pfAll->vz());
    vc1RFVsPFx_->Fill(1*vc1RF_phase, pfAll->vx());
    vc1RFVsPFy_->Fill(1*vc1RF_phase, pfAll->vy());
    vc1RFVsPFz_->Fill(1*vc1RF_phase, pfAll->vz());
    vc2RFVsPFx_->Fill(1*vc2RF_phase, pfAll->vx());
    vc2RFVsPFy_->Fill(1*vc2RF_phase, pfAll->vy());
    vc2RFVsPFz_->Fill(1*vc2RF_phase, pfAll->vz());
  }

  for ( const auto& rechits : *diamondRecHits ) {
    const CTPPSDiamondDetId detId( rechits.detId() );
    for ( const auto& rechit : rechits ) {
      if (rechit.getOOTIndex() == 0 ) {
        if ( detId.arm() == 0 ) {
          ppsDiamond45VsBx_->Fill(iEvent.bunchCrossing(), rechit.getT()-5.351);
          ppsDiamond45VsBx_low_->Fill(iEvent.bunchCrossing(), rechit.getT()-5.351);
          ppsDiamond45VsBx_corrected_->Fill(iEvent.bunchCrossing(), rechit.getT()-5.351 + ave_phase);
          ppsDiamond45VsBx_corrected_low_->Fill(iEvent.bunchCrossing(), rechit.getT()-5.351 + ave_phase);
          ppsDiamond45VsBx_correctedRF_->Fill(iEvent.bunchCrossing(), rechit.getT()-5.351 + aveRF_phase);
          ppsDiamond45VsBx_correctedRF_low_->Fill(iEvent.bunchCrossing(), rechit.getT()-5.351 + aveRF_phase);
        }
        else {
          ppsDiamond56VsBx_->Fill(iEvent.bunchCrossing(), rechit.getT()-4.789);
          ppsDiamond56VsBx_low_->Fill(iEvent.bunchCrossing(), rechit.getT()-4.789);
          ppsDiamond56VsBx_corrected_->Fill(iEvent.bunchCrossing(), rechit.getT()-4.789 + ave_phase);
          ppsDiamond56VsBx_corrected_low_->Fill(iEvent.bunchCrossing(), rechit.getT()-4.789 + ave_phase);
          ppsDiamond56VsBx_correctedRF_->Fill(iEvent.bunchCrossing(), rechit.getT()-4.789 + aveRF_phase);
          ppsDiamond56VsBx_correctedRF_low_->Fill(iEvent.bunchCrossing(), rechit.getT()-4.789 + aveRF_phase);
        }
      }
    }
  }

  // for(auto itb=barrelHitHandle->begin(); itb!=barrelHitHandle->end(); ++itb){
  //   // EBDetId id(itb->id());
  //   double time = itb->time();
  //
  // }
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
  aveVsaveRF_.Write();
  difVsdifRF_.Write();
  vc1Vsvc1RF_.Write();
  vc2Vsvc2RF_.Write();
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
