// -*- C++ -*-
//
// Package:    NickAnlzr/MyAnalyzer
// Class:      MyAnalyzer
// 
/**\class MyAnalyzer MyAnalyzer.cc NickAnlzr/MyAnalyzer/plugins/MyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nicola Minafra
//         Created:  Thu, 05 Oct 2017 10:14:10 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/ServiceRegistry/interface/Service.h"
 #include "CommonTools/UtilAlgos/interface/TFileService.h"
 
 #include "TH1.h"
 #include "TF1.h"
 #include "TH2.h"
 #include "TGraph.h"
//  #include<map.h>
 
#include "DataFormats/Provenance/interface/EventRange.h"
// #include "DataFormats/CTPPSDigi/interface/TotemVFATStatus.h"
#include "DataFormats/CTPPSDigi/interface/TotemFEDInfo.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"

#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"

#define COS_8_DEG 0.990268
#define SIN_8_DEG -0.139173
#define STRIP_DIAM_SHIFT_MM 4.9
#define PIXEL_DIAM_SHIFT_MM 43.8
 
class MyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MyAnalyzer(const edm::ParameterSet&);
      ~MyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
//       edm::EDGetTokenT< edm::DetSetVector<TotemVFATStatus> > tokenStatus_;
      edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > tokenLocalTrack_;
      edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondDigi> > tokenDigi_;
      edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenDiamondHit_;
      edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenDiamondTrack_;
//       edm::EDGetTokenT< std::vector<TotemFEDInfo> > tokenFEDInfo_;
      edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack> >  pixelLocalTrackToken_;
    
      
      TFileDirectory subdir_eff_, subdir_hits_, subdir_timing_, subdir_tot_, subdir_tot_single_, subdir_tracks_, subdir_strips_, subdir_pixel_, subdir_le_;
    
      unsigned int verbosity_;
      std::unordered_map<unsigned int, TH1D* > channelPlots_00_, channelPlots_01_, channelPlots_10_, channelPlots_11_, channelPlots_all_;
      
      TH2D *histo_45_, *histo_56_;
      std::vector< TH2D* > histo_tarm_vec_;
      TH2D *correlation_45_, *correlation_56_;
      std::vector< TH2D* > correlation_vec_;
      TH1D *histo_delta_;
      
      std::map<int, TGraph*> tot_graphsMap_, tot_graphsMap_single_;
      std::map<int, int> tot_counter_;
      
      std::map< int, TH1D*> histo_eff_map_;
      std::map<int, int> eff_triplecounting_map_;
      std::map<int, int> eff_doublecounting_map_;
      
      std::map< int, TH1D*> histo_hits_map_;
      
      TGraph *strip_diam_corr_45_, *strip_diam_corr_56_;
      std::vector< TGraph* > strip_diam_corr_vec_;
      std::vector< int > strip_diam_corr_counter_vec_;
      
      TH1D *strip_diam_corr_histo_45_, *strip_diam_corr_histo_56_;
      std::vector< TH1D* > strip_diam_corr_histo_vec_;
      
      TH2D *strip_45_, *strip_56_;
      std::vector< TH2D* > strip_vec_;
      
      TH2D *strip_all_45_, *strip_all_56_;
      std::vector< TH2D* > strip_all_vec_;
      
      TH2D *strip_anti_45_, *strip_anti_56_;
      std::vector< TH2D* > strip_anti_vec_;
      
      TGraph *pixel_diam_corr_45_, *pixel_diam_corr_56_;
      std::vector< TGraph* > pixel_diam_corr_vec_;
      std::vector< int > pixel_diam_corr_counter_vec_;
      
      TH1D *pixel_diam_corr_histo_45_, *pixel_diam_corr_histo_56_;
      std::vector< TH1D* > pixel_diam_corr_histo_vec_;
      
      TH2D *pixel_45_, *pixel_56_;
      std::vector< TH2D* > pixel_vec_;
      
      TH2D *pixel_all_45_, *pixel_all_56_;
      std::vector< TH2D* > pixel_all_vec_;
      
      TH2D *pixel_anti_45_, *pixel_anti_56_;
      std::vector< TH2D* > pixel_anti_vec_;
      
      int trackDisplayCounter_;
      std::string trackDisplayStr_;
      
      std::map< int, TH1D*> time_delay_map_;
      std::map< int, TH1D*> time_delay_cumulative_map_;
      int selectedOOTIndex_;
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
MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig) :
  tokenLocalTrack_  ( consumes< edm::DetSetVector<TotemRPLocalTrack> >     ( iConfig.getParameter<edm::InputTag>( "tagLocalTrack" ) ) ),
  tokenDigi_        ( consumes< edm::DetSetVector<CTPPSDiamondDigi> >      ( iConfig.getParameter<edm::InputTag>( "tagDigi" ) ) ),
  tokenDiamondHit_  ( consumes< edm::DetSetVector<CTPPSDiamondRecHit> >    ( iConfig.getParameter<edm::InputTag>( "tagDiamondRecHits" ) ) ),
  tokenDiamondTrack_( consumes< edm::DetSetVector<CTPPSDiamondLocalTrack> >( iConfig.getParameter<edm::InputTag>( "tagDiamondLocalTracks" ) ) ),
  pixelLocalTrackToken_ ( consumes<edm::DetSetVector<CTPPSPixelLocalTrack> >( iConfig.getParameter<edm::InputTag>("ctppsPixelLocalTracks" ) ) ),
  verbosity_                     ( iConfig.getUntrackedParameter<unsigned int>( "verbosity", 0 ) ),
  trackDisplayCounter_(0),
  trackDisplayStr_("track_"),
  selectedOOTIndex_( iConfig.getParameter<int>( "selectedOOTIndex" ) )
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   subdir_eff_ = fs->mkdir( "Efficiency" );
   subdir_hits_ = fs->mkdir( "HitDistributions" );
   subdir_strips_ = fs->mkdir( "Strips" );
   subdir_pixel_ = fs->mkdir( "Pixel" );
   subdir_timing_ = fs->mkdir( "Timing" );
   subdir_tracks_ = fs->mkdir( "Tracks" );
   subdir_tot_ = subdir_timing_.mkdir( "Tot" );
   subdir_tot_single_ = subdir_timing_.mkdir( "Tot_single" );
   subdir_le_ = subdir_timing_.mkdir( "Leading_Edges" );
   
   histo_45_ = subdir_timing_.make<TH2D>("tracks_45_t" , "tracks_45_t" , 10 , 0 , 10, 401, -2, 2 );
   histo_tarm_vec_.emplace_back( histo_45_ );
   histo_56_ = subdir_timing_.make<TH2D>("tracks_56_t" , "tracks_56_t" , 10 , 0 , 10, 401, -2, 2 );
   histo_tarm_vec_.emplace_back( histo_56_ );
   
   correlation_45_ = subdir_tracks_.make<TH2D>("correlation_45_t" , "tracks pl vs hits 45" , 10 , 0 , 10, 10, 0, 10 );
   correlation_45_->GetXaxis()->SetTitle( "number of planes per track" );
   correlation_45_->GetYaxis()->SetTitleOffset(1);
   correlation_45_->GetYaxis()->SetTitle( "number of hits per track" );
   correlation_vec_.emplace_back( correlation_45_ );
   correlation_56_ = subdir_tracks_.make<TH2D>("correlation_56_t" , "tracks pl vs hits 56" , 10 , 0 , 10, 10, 0, 10 );
   correlation_56_->GetXaxis()->SetTitle( "number of planes per track" );
   correlation_56_->GetYaxis()->SetTitleOffset(1);
   correlation_56_->GetYaxis()->SetTitle( "number of hits per track" );
   correlation_vec_.emplace_back( correlation_56_ );
   
   histo_delta_ = subdir_timing_.make<TH1D>("histo_delta_t" , "histo_delta_t" , 401, -2, 2 );
   histo_delta_->GetXaxis()->SetTitle( "t_56 - t_45 (ns)" );
   histo_delta_->GetYaxis()->SetTitleOffset(1);
   
   strip_diam_corr_45_ = subdir_strips_.make<TGraph>();
   strip_diam_corr_45_->SetName("strip_diam_corr_45");
   strip_diam_corr_45_->SetTitle("strip diam corr 45");
   strip_diam_corr_45_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_diam_corr_45_->GetYaxis()->SetTitleOffset(1);
   strip_diam_corr_45_->GetYaxis()->SetTitle( "diamond x (mm)" );
   strip_diam_corr_vec_.emplace_back( strip_diam_corr_45_ );
   strip_diam_corr_counter_vec_.emplace_back(0);
   
   strip_diam_corr_56_ = subdir_strips_.make<TGraph>();
   strip_diam_corr_56_->SetName("strip_diam_corr_56");
   strip_diam_corr_56_->SetTitle("strip diam corr 56");
   strip_diam_corr_56_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_diam_corr_56_->GetYaxis()->SetTitleOffset(1);
   strip_diam_corr_56_->GetYaxis()->SetTitle( "diamond x (mm)" );
   strip_diam_corr_vec_.emplace_back( strip_diam_corr_56_ );
   strip_diam_corr_counter_vec_.emplace_back(0);
   
   strip_diam_corr_histo_45_ = subdir_strips_.make<TH1D>("strip_diam_corr_histo_45", "strip_diam_corr_histo_45", 200,-10,10);
   strip_diam_corr_histo_45_->GetXaxis()->SetTitle( "strip - shift - diamond (mm)" );
   strip_diam_corr_histo_45_->GetYaxis()->SetTitleOffset(1);
   strip_diam_corr_histo_vec_.emplace_back( strip_diam_corr_histo_45_ );
   
   strip_diam_corr_histo_56_ = subdir_strips_.make<TH1D>("strip_diam_corr_histo_56", "strip_diam_corr_histo_56", 200,-10,10);
   strip_diam_corr_histo_56_->GetXaxis()->SetTitle( "strip - shift - diamond (mm)" );
   strip_diam_corr_histo_56_->GetYaxis()->SetTitleOffset(1);
   strip_diam_corr_histo_vec_.emplace_back( strip_diam_corr_histo_56_ );
   
   strip_all_45_ = subdir_strips_.make<TH2D>("strip_all_45_","Strips tr all distribution 45", 50,0,20,50,-10,10 );
   strip_all_45_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_all_45_->GetYaxis()->SetTitleOffset(1);
   strip_all_45_->GetYaxis()->SetTitle( "strip y (mm)" );
   strip_all_vec_.emplace_back( strip_all_45_ );
   
   strip_all_56_ = subdir_strips_.make<TH2D>("strip_all_56_","Strips tr all distribution 56", 50,0,20,50,-10,10 );
   strip_all_56_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_all_56_->GetYaxis()->SetTitleOffset(1);
   strip_all_56_->GetYaxis()->SetTitle( "strip y (mm)" );
   strip_all_vec_.emplace_back( strip_all_56_ );
   
   strip_45_ = subdir_strips_.make<TH2D>("strip_45","Strips tr distribution 45", 50,0,20,50,-10,10 );
   strip_45_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_45_->GetYaxis()->SetTitleOffset(1);
   strip_45_->GetYaxis()->SetTitle( "strip y (mm)" );
   strip_vec_.emplace_back( strip_45_ );
   
   strip_56_ = subdir_strips_.make<TH2D>("strip_56","Strips tr distribution 56", 50,0,20,50,-10,10 );
   strip_56_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_56_->GetYaxis()->SetTitleOffset(1);
   strip_56_->GetYaxis()->SetTitle( "strip y (mm)" );
   strip_vec_.emplace_back( strip_56_ );
   
   strip_anti_45_ = subdir_strips_.make<TH2D>("strip_anti_45","Strips tr anti distribution 45", 50,0,20,50,-10,10 );
   strip_anti_45_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_anti_45_->GetYaxis()->SetTitleOffset(1);
   strip_anti_45_->GetYaxis()->SetTitle( "strip y (mm)" );
   strip_anti_vec_.emplace_back( strip_anti_45_ );
   
   strip_anti_56_ = subdir_strips_.make<TH2D>("strip_anti_56","Strips tr anti distribution 56", 50,0,20,50,-10,10 );
   strip_anti_56_->GetXaxis()->SetTitle( "strip x (mm)" );
   strip_anti_56_->GetYaxis()->SetTitleOffset(1);
   strip_anti_56_->GetYaxis()->SetTitle( "strip y (mm)" );
   strip_anti_vec_.emplace_back( strip_anti_56_ );
   
   
   
   pixel_diam_corr_45_ = subdir_pixel_.make<TGraph>();
   pixel_diam_corr_45_->SetName("pixel_diam_corr_45");
   pixel_diam_corr_45_->SetTitle("pixel diam corr 45");
   pixel_diam_corr_45_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_diam_corr_45_->GetYaxis()->SetTitleOffset(1);
   pixel_diam_corr_45_->GetYaxis()->SetTitle( "diamond x (mm)" );
   pixel_diam_corr_vec_.emplace_back( pixel_diam_corr_45_ );
   pixel_diam_corr_counter_vec_.emplace_back(0);
   
   pixel_diam_corr_56_ = subdir_pixel_.make<TGraph>();
   pixel_diam_corr_56_->SetName("pixel_diam_corr_56");
   pixel_diam_corr_56_->SetTitle("pixel diam corr 56");
   pixel_diam_corr_56_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_diam_corr_56_->GetYaxis()->SetTitleOffset(1);
   pixel_diam_corr_56_->GetYaxis()->SetTitle( "diamond x (mm)" );
   pixel_diam_corr_vec_.emplace_back( pixel_diam_corr_56_ );
   pixel_diam_corr_counter_vec_.emplace_back(0);
   
   pixel_diam_corr_histo_45_ = subdir_pixel_.make<TH1D>("pixel_diam_corr_histo_45", "pixel_diam_corr_histo_45", 200,-10,10);
   pixel_diam_corr_histo_45_->GetXaxis()->SetTitle( "pixel - shift - diamond (mm)" );
   pixel_diam_corr_histo_45_->GetYaxis()->SetTitleOffset(1);
   pixel_diam_corr_histo_vec_.emplace_back( pixel_diam_corr_histo_45_ );
   
   pixel_diam_corr_histo_56_ = subdir_pixel_.make<TH1D>("pixel_diam_corr_histo_56", "pixel_diam_corr_histo_56", 200,-10,10);
   pixel_diam_corr_histo_56_->GetXaxis()->SetTitle( "pixel - shift - diamond (mm)" );
   pixel_diam_corr_histo_56_->GetYaxis()->SetTitleOffset(1);
   pixel_diam_corr_histo_vec_.emplace_back( pixel_diam_corr_histo_56_ );
   
   pixel_45_ = subdir_pixel_.make<TH2D>("pixel_45","Pixel tr distribution 45", 50,0,20,50,-10,10 );
   pixel_45_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_45_->GetYaxis()->SetTitleOffset(1);
   pixel_45_->GetYaxis()->SetTitle( "pixel y (mm)" );
   pixel_vec_.emplace_back( pixel_45_ );
   
   pixel_56_ = subdir_pixel_.make<TH2D>("pixel_56","Pixel tr distribution 56", 50,0,20,50,-10,10 );
   pixel_56_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_56_->GetYaxis()->SetTitleOffset(1);
   pixel_56_->GetYaxis()->SetTitle( "pixel y (mm)" );
   pixel_vec_.emplace_back( pixel_56_ );
   
   pixel_all_45_ = subdir_pixel_.make<TH2D>("pixel_all_45","Pixel tr distribution 45", 50,0,20,50,-10,10 );
   pixel_all_45_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_all_45_->GetYaxis()->SetTitleOffset(1);
   pixel_all_45_->GetYaxis()->SetTitle( "pixel y (mm)" );
   pixel_all_vec_.emplace_back( pixel_all_45_ );
   
   pixel_all_56_ = subdir_pixel_.make<TH2D>("pixel_all_56","Pixel tr distribution 56", 50,0,20,50,-10,10 );
   pixel_all_56_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_all_56_->GetYaxis()->SetTitleOffset(1);
   pixel_all_56_->GetYaxis()->SetTitle( "pixel y (mm)" );
   pixel_all_vec_.emplace_back( pixel_all_56_ );
   
   pixel_anti_45_ = subdir_pixel_.make<TH2D>("pixel_anti_45","Pixel tr anti distribution 45", 50,0,20,50,-10,10 );
   pixel_anti_45_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_anti_45_->GetYaxis()->SetTitleOffset(1);
   pixel_anti_45_->GetYaxis()->SetTitle( "pixel y (mm)" );
   pixel_anti_vec_.emplace_back( pixel_anti_45_ );
   
   pixel_anti_56_ = subdir_pixel_.make<TH2D>("pixel_anti_56","Pixel tr anti distribution 56", 50,0,20,50,-10,10 );
   pixel_anti_56_->GetXaxis()->SetTitle( "pixel x (mm)" );
   pixel_anti_56_->GetYaxis()->SetTitleOffset(1);
   pixel_anti_56_->GetYaxis()->SetTitle( "pixel y (mm)" );
   pixel_anti_vec_.emplace_back( pixel_anti_56_ );
   
     
   std::cout << "Counter \tNumOfHits \t\t x \t\tsigma" << std::endl;

}


MyAnalyzer::~MyAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{   
    edm::Service<TFileService> fs;
    
    // get event data
//     edm::Handle< edm::DetSetVector<TotemVFATStatus> > diamondVFATStatus;
//     event.getByToken( tokenStatus_, diamondVFATStatus );
// 
    edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > stripTracks;
    event.getByToken( tokenLocalTrack_, stripTracks );

    edm::Handle< edm::DetSetVector<CTPPSDiamondDigi> > diamondDigis;
    event.getByToken( tokenDigi_, diamondDigis );

//     edm::Handle< std::vector<TotemFEDInfo> > fedInfo;
//     event.getByToken( tokenFEDInfo_, fedInfo );
// 
    edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > diamondRecHits;
    event.getByToken( tokenDiamondHit_, diamondRecHits );

    edm::Handle< edm::DetSetVector<CTPPSDiamondLocalTrack> > diamondLocalTracks;
    event.getByToken( tokenDiamondTrack_, diamondLocalTracks );
    
    edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > pixelLocalTracks;
    event.getByToken(pixelLocalTrackToken_, pixelLocalTracks);

    // check validity
    bool valid = true;
//     valid &= diamondVFATStatus.isValid();
    valid &= diamondDigis.isValid();
//     valid &= fedInfo.isValid();

    if ( !valid ) {
      if ( verbosity_ ) {
        edm::LogProblem("CTPPSDiamondDQMSource")
          << "ERROR in TotemDQMModuleRP::analyze > some of the required inputs are not valid. Skipping this event.\n"
//           << "    diamondVFATStatus.isValid = " << diamondVFATStatus.isValid() << "\n"
          << "    diamondDigis.isValid = " << diamondDigis.isValid() << "\n";
//           << "    fedInfo.isValid = " << fedInfo.isValid();
      }

      return;
    }
    
    CTPPSDiamondLocalTrack last_track_45, last_track_56;
    int tracks_45_counter(0), tracks_56_counter(0);
      
    for (unsigned int arm=0; arm<=1; ++arm ) {
      
      TotemRPLocalTrack stripTrack;
      
      //Select only events with one track in the strips
      int strip_tracks_counter = 0;
      for ( const auto& ds1 : *stripTracks ) {
        for ( const auto& tr1 : ds1 ) {
          if ( ! tr1.isValid() )  continue;
          CTPPSDetId rpId1( ds1.detId() );
          if ( rpId1.arm() == arm ) {
            stripTrack = tr1;
            ++strip_tracks_counter;
          }
        }
      }
//       if ( strip_tracks_counter != 1 ) continue;     //Select only single tracks (stips) events. 
      
      
      CTPPSPixelLocalTrack pixelTrack;
      //Select only events with one track in the strips
      int pixel_tracks_counter = 0;
      for ( const auto& ds1 : *pixelLocalTracks ) {
        for ( const auto& tr1 : ds1 ) {
          if ( ! tr1.isValid() )  continue;
          CTPPSPixelDetId rpId1( ds1.detId() );
          if ( rpId1.arm() == arm ) {
            pixelTrack = tr1;
            ++pixel_tracks_counter;
            pixel_all_vec_[arm]->Fill( pixelTrack.getX0()-PIXEL_DIAM_SHIFT_MM, pixelTrack.getY0() );
          }
        }
      }
      if ( pixel_tracks_counter != 1 ) continue;     //Select only single tracks (pixels) events. 
      
        
      // Separate tracks based on num of hits
      int tracks_counter = 0;
      for ( const auto& tracks : *diamondLocalTracks ) {
        const CTPPSDiamondDetId detId( tracks.detId() );
        for ( const auto& track : tracks ) {
          if ( detId.arm() == arm ) {
            histo_tarm_vec_[arm]->Fill( track.getNumOfHits(), track.getT() );
            ++tracks_counter;
            
            // Strips
            double stripRotatedX = COS_8_DEG * stripTrack.getX0() - SIN_8_DEG * stripTrack.getY0();
            double stripRotatedY = SIN_8_DEG * stripTrack.getX0() + COS_8_DEG * stripTrack.getY0();
            
            if ( stripRotatedY > -20 && stripRotatedY < 20 ) {
              strip_diam_corr_vec_[arm]->SetPoint( (strip_diam_corr_counter_vec_[arm])++, stripRotatedX-STRIP_DIAM_SHIFT_MM, track.getX0() );
              strip_diam_corr_histo_vec_[arm]->Fill( stripRotatedX-STRIP_DIAM_SHIFT_MM - track.getX0() );
            }
            
            strip_all_vec_[arm]->Fill( stripRotatedX, stripRotatedY );
            
            if ( std::abs( stripRotatedX-STRIP_DIAM_SHIFT_MM - track.getX0() ) < 1 )
              strip_vec_[arm]->Fill( stripRotatedX, stripRotatedY );
            else
              strip_anti_vec_[arm]->Fill( stripRotatedX, stripRotatedY );
            
            
            // Pixels
            if ( pixelTrack.getY0() > -20 && pixelTrack.getY0() < 20 ) {
              pixel_diam_corr_vec_[arm]->SetPoint( (pixel_diam_corr_counter_vec_[arm])++, pixelTrack.getX0()-PIXEL_DIAM_SHIFT_MM, track.getX0() );
              pixel_diam_corr_histo_vec_[arm]->Fill( pixelTrack.getX0()-PIXEL_DIAM_SHIFT_MM - track.getX0() );
            }
            
            if ( std::abs( pixelTrack.getX0()-PIXEL_DIAM_SHIFT_MM - track.getX0() ) < 1 )
              pixel_vec_[arm]->Fill( pixelTrack.getX0()-PIXEL_DIAM_SHIFT_MM, pixelTrack.getY0() );
            else
              pixel_anti_vec_[arm]->Fill( pixelTrack.getX0()-PIXEL_DIAM_SHIFT_MM, pixelTrack.getY0() );
            
            
            
            if ( arm == 0 ) {
              last_track_45 = track;
              tracks_45_counter = tracks_counter;
            }
            else {
              last_track_56 = track;
              tracks_56_counter = tracks_counter;
            }
            correlation_vec_[arm]->Fill( track.getNumOfPlanes(), track.getNumOfHits() );
          }
        }
      }
      
      if ( tracks_counter > 1 ) continue;     //Select only single tracks (diamond) events. 
      
      for ( const auto& tracks : *diamondLocalTracks ) {
        const CTPPSDiamondDetId detId( tracks.detId() );
        if ( detId.arm() != arm ) continue;
        for ( const auto& track : tracks ) {
          
          // Efficiency
//           std::cout<<"Checking " << detId << std::endl;
          int width_eff=1;
          double tollerance=0.5;
          for (int x_eff=0; x_eff<15; x_eff+=width_eff) {
            if ( track.getMultipleHits()==0 && track.getNumOfHits() > 0 && track.getNumOfHits() <= 10 && track.getNumOfPlanes() > 1  && track.getX0() > x_eff && track.getX0() < x_eff+width_eff) {         
              // Check which planes participated
//               std::cout<<"Checking " << x_eff << std::endl;
              std::set<unsigned int> planes_set;
              for ( const auto& rechits : *diamondRecHits ) {
                CTPPSDiamondDetId detId_hit( rechits.detId() );
                for ( const auto& rechit : rechits ) {
                  if ( detId_hit.arm() == detId.arm() ) {
                    if ( CTPPSHitBelongsToTrack(track,rechit,tollerance) ) {
                      planes_set.insert(detId_hit.plane());
                    }
                  }
                }
              }
              if ( (int) planes_set.size() != track.getNumOfPlanes() ) std::cout << "ERRORR!! different numb of planes!!!!" << std::endl;
              
              // Allow reco with ufsd: if there are only two planes and one id ufsd
              if ( planes_set.find( 3 ) != planes_set.end() and planes_set.size() < 3 ) continue;
              
              for (unsigned int plane=0; plane<4; ++plane) {
                int map_index = detId.arm()*10 + plane + 1000*x_eff;
//                 std::cout<<"Looking " << map_index << std::endl;
                if ( eff_doublecounting_map_.find( map_index ) == eff_doublecounting_map_.end() ) {
                  eff_triplecounting_map_[map_index] = 0;
                  eff_doublecounting_map_[map_index] = 0;
                  std::cout<<"Created" << std::endl;
                }
                ++(eff_doublecounting_map_[map_index]);
              }
              
              for (auto& plane : planes_set ) {
                int map_index = detId.arm()*10 + plane + 1000*x_eff;
                ++(eff_triplecounting_map_[map_index]);
              }
              
              
/*              
              
              
              for (int plane_i=0; plane_i<4; ++plane_i) {
                unsigned int plane = plane_i;
                int map_index = detId.arm()*10 + plane + 1000*x_eff;
                

                ++(eff_doublecounting_map_[map_index]);

                bool plane_has_hit=false;
                for ( const auto& rechits : *diamondRecHits ) {
                  CTPPSDiamondDetId detId_hit( rechits.detId() );
                  for ( const auto& rechit : rechits ) {
                    if ( detId_hit.arm() == detId.arm() && detId_hit.plane() == plane &&  track.getNumOfPlanes() == 4 ) {
                      if ( CTPPSHitBelongsToTrack(track,rechit,1) ) {
                        plane_has_hit=true;
                      }
                    }
                  }
                }
                if ( plane_has_hit ) ++(eff_triplecounting_map_[map_index]);
              }*/
            }
          }
          
          //Time Over Threshold
          if ( track.getNumOfHits() >= 3 && track.getNumOfHits() <= 4 && track.getNumOfPlanes() > 2 ) {
            for ( const auto& rechits : *diamondRecHits ) {
              CTPPSDiamondDetId detId_hit( rechits.detId() );
              for ( const auto& rechit : rechits ) {
                if ( detId_hit.arm() == detId.arm() ) {
                  if ( CTPPSHitBelongsToTrack(track,rechit) &&  ( rechit.getToT() > 10 && rechit.getToT() < 20 ) ) {
                    int map_index = detId_hit.arm()*1000 + detId_hit.plane()*100 + detId_hit.channel();
                    if ( tot_graphsMap_.find( map_index ) == tot_graphsMap_.end() ) {
                      tot_graphsMap_[map_index] = subdir_tot_.make<TGraph>();
                      std::string title;
                      detId_hit.channelName( title, CTPPSDiamondDetId::nFull );
                      tot_graphsMap_[map_index]->SetName( title.c_str() );
                      tot_graphsMap_[map_index]->SetTitle( title.c_str() );
                      tot_graphsMap_[map_index]->GetXaxis()->SetTitle( "tot (ns)" );
                      tot_graphsMap_[map_index]->GetYaxis()->SetTitle( "t_hit - t_track (ns)" );
                      tot_graphsMap_[map_index]->GetYaxis()->SetTitleOffset(1);
                      tot_counter_[map_index] = 0;
                      
                      tot_graphsMap_single_[map_index] = subdir_tot_single_.make<TGraph>();
                      tot_graphsMap_single_[map_index]->SetName( title.c_str() );
                      tot_graphsMap_single_[map_index]->SetTitle( title.c_str() );
                      tot_graphsMap_single_[map_index]->GetXaxis()->SetTitle( "tot (ns)" );
                      tot_graphsMap_single_[map_index]->GetYaxis()->SetTitle( "t_hit - t_track (ns)" );
                      tot_graphsMap_single_[map_index]->GetYaxis()->SetTitleOffset(1);
                    }
                    else {
                      tot_graphsMap_[map_index]->SetPoint( tot_counter_[map_index], rechit.getToT(), rechit.getT() - track.getT() );
                      tot_graphsMap_single_[map_index]->SetPoint( (tot_counter_[map_index])++, rechit.getToT(), rechit.getT() );
                    }    
                  }
                }
              }
            }
          }
          
          
          // Track residuals
//           if ( track.getNumOfHits() >= 2) {  
//             std::string title(trackDisplayStr_);
//             if ( detId.arm() == 0 ) title+="45_";
//             else title+="56_";
//             title+=std::to_string(trackDisplayCounter_++);
//             TH2I* histo_tmp = subdir_tracks_.make<TH2I>( title.c_str(), title.c_str(), 4, -0.5, 3.5, 19.*10, -1, 18 );
//                         
//             for ( const auto& rechits : *diamondRecHits ) {
//               CTPPSDiamondDetId detId_hit( rechits.detId() );
//               if ( detId_hit.arm() == detId.arm() ) {
//                 for ( const auto& rechit : rechits ) {
//                   if ( CTPPSHitBelongsToTrack(track,rechit) ) {
//   //                   histo_tmp->Fill( detId_hit.plane(), detId_hit.channel() );
// //                     std::cout<<"\t\t\t\t\t"<<  detId_hit.plane() << "\t" << detId_hit.channel() << "\t" << rechit.getT() << std::endl;
//                     TAxis *hitHistoOOTTmpYAxis = histo_tmp->GetYaxis();
//                     int startBin = hitHistoOOTTmpYAxis->FindBin( rechit.getX() - 0.5*rechit.getXWidth() );
//                     int numOfBins = rechit.getXWidth()*10;
//                     for ( int i=0; i<numOfBins; ++i) {
//                       histo_tmp->Fill( detId_hit.plane(), hitHistoOOTTmpYAxis->GetBinCenter(startBin+i) );
//                     }
//                   }
//                 }
//               }
//             }
//           }
          
          
          // Track Display!
          if ( track.getNumOfPlanes() < 2 ) {
            std::cout << trackDisplayCounter_ << "\t\t" << track.getNumOfHits() << "\t\t" << track.getX0() << "\t\t" << track.getX0Sigma() << std::endl;
            
            std::string title(trackDisplayStr_);
            if ( detId.arm() == 0 ) title+="45_";
            else title+="56_";
            title+=std::to_string(trackDisplayCounter_++);
            TH2I* histo_tmp = subdir_tracks_.make<TH2I>( title.c_str(), title.c_str(), 4, -0.5, 3.5, 19.*10, -1, 18 );
            histo_tmp->GetXaxis()->SetTitle( "plane" );
            histo_tmp->GetYaxis()->SetTitle( "x (mm)" );
            histo_tmp->GetYaxis()->SetTitleOffset(1);
                        
            for ( const auto& rechits : *diamondRecHits ) {
              CTPPSDiamondDetId detId_hit( rechits.detId() );
              if ( detId_hit.arm() == detId.arm() ) {
                for ( const auto& rechit : rechits ) {
                  if ( CTPPSHitBelongsToTrack(track,rechit) ) {
  //                   histo_tmp->Fill( detId_hit.plane(), detId_hit.channel() );
//                     std::cout<<"\t\t\t\t\t"<<  detId_hit.plane() << "\t" << detId_hit.channel() << "\t" << rechit.getT() << std::endl;
                    TAxis *hitHistoOOTTmpYAxis = histo_tmp->GetYaxis();
                    int startBin = hitHistoOOTTmpYAxis->FindBin( rechit.getX() - 0.5*rechit.getXWidth() );
                    int numOfBins = rechit.getXWidth()*10;
                    for ( int i=0; i<numOfBins; ++i) {
                      histo_tmp->Fill( detId_hit.plane(), hitHistoOOTTmpYAxis->GetBinCenter(startBin+i) );
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      
    } // End of arm cycle
    
    // Leading edge plots for delay computation
    for ( const auto& rechits : *diamondRecHits ) {
      CTPPSDiamondDetId detId( rechits.detId() );
      for ( const auto& rechit : rechits ) {
        int tempId = 1000*detId.arm() + 100*detId.plane() + detId.channel();
        if ( time_delay_map_.find( tempId ) == time_delay_map_.end() ) {
          std::string title;
          detId.channelName( title, CTPPSDiamondDetId::nFull );
          time_delay_map_[ tempId ] = subdir_le_.make<TH1D>(title.c_str(), title.c_str(), 501,-25,25);
          time_delay_map_[ tempId ]->GetXaxis()->SetTitle( "t_hit (ns)" );
        }
        if ( time_delay_cumulative_map_.find( detId.arm() ) == time_delay_cumulative_map_.end() ) {
          std::string title("Cumulative_arm");
          title += std::to_string(detId.arm());
          time_delay_cumulative_map_[ detId.arm() ] = subdir_le_.make<TH1D>(title.c_str(), title.c_str(), 501,-25,25);
          time_delay_cumulative_map_[ detId.arm() ]->GetXaxis()->SetTitle( "t_hit (ns)" );
        }
        if ( rechit.getOOTIndex() == selectedOOTIndex_ ) {
          time_delay_map_[ tempId ]->Fill( rechit.getT()  );
          time_delay_cumulative_map_[ detId.arm() ]->Fill( rechit.getT()  );
        }
      }
    }
    
    //Hit maps per plane
    for ( const auto& rechits : *diamondRecHits ) {
      CTPPSDiamondDetId detId( rechits.detId() );
      int map_index = detId.arm()*10 + detId.plane();
      if (histo_hits_map_.find( map_index ) == histo_hits_map_.end() ) {
        std::string title("HitMap_arm");
        title += std::to_string( detId.arm() );
        title += "_plane";
        title += std::to_string( detId.plane() );
        histo_hits_map_[ map_index ] = subdir_hits_.make<TH1D>(title.c_str(), title.c_str(), 500, 0, 25);
        histo_hits_map_[ map_index ]->GetXaxis()->SetTitle( "x (mm)" );
        histo_hits_map_[ map_index ]->GetYaxis()->SetTitle( "# of hits" );
        histo_hits_map_[ map_index ]->GetYaxis()->SetTitleOffset(1);
      }
      
      for ( const auto& rechit : rechits ) {
        TAxis *hitHisto = histo_hits_map_[ map_index ]->GetXaxis();
        int startBin = hitHisto->FindBin( rechit.getX() - 0.5*rechit.getXWidth() );
        int stopBin = hitHisto->FindBin( rechit.getX() + 0.5*rechit.getXWidth() );
        for ( int i=0; i<stopBin-startBin; ++i) {
          histo_hits_map_[ map_index ]->Fill( hitHisto->GetBinCenter(startBin+i) );
        }
        
      }
    }
    
    
    
    
    

    if ( tracks_45_counter==1 && tracks_56_counter==1 && last_track_45.getNumOfHits() == 3 && last_track_56.getNumOfHits() == 3 )
      histo_delta_->Fill( last_track_56.getT() - last_track_45.getT() );
}


// ------------ method called once each job just before starting event loop  ------------
void 
MyAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyAnalyzer::endJob() 
{
  for ( auto& element : eff_triplecounting_map_ ) { 
    if ( eff_doublecounting_map_[ element.first ] > 0) {
      
//       int map_index = detId.arm()*10 + plane + 1000*x_eff;
      int arm = element.first % 100 / 10;
      int plane = element.first % 10;
      int x_eff = element.first / 1000;
      int map_index = arm*10 + plane;
      
      if (histo_eff_map_.find( map_index ) == histo_eff_map_.end() ) {
        std::string title("Efficiency_arm");
        title += std::to_string( arm );
        title += "_plane";
        title += std::to_string( plane );
        histo_eff_map_[ map_index ] = subdir_eff_.make<TH1D>(title.c_str(), title.c_str(), 25, 0, 25);
        histo_eff_map_[ map_index ]->GetXaxis()->SetTitle( "x (mm)" );
        histo_eff_map_[ map_index ]->GetYaxis()->SetTitle( "Efficiency %" );
        histo_eff_map_[ map_index ]->GetYaxis()->SetTitleOffset(1);
      }
      
      
      
      double counted = element.second;
      double total = eff_doublecounting_map_[ element.first ];
      double efficiency = counted / total;
      double error = std::sqrt( efficiency * ( 1 - efficiency ) / total );
      
      histo_eff_map_[ map_index ]->SetBinContent( x_eff + 1 , efficiency );
      histo_eff_map_[ map_index ]->SetBinError( x_eff + 1 , error );
      std::cout << "Arm: " << arm << "\tplane: " << plane << "\tx: " << x_eff << "\t\t" << counted << "\tout of " << total << "\t\t" << 100. * efficiency << "% +- " << 100. * error << std::endl;
    }
  }
  
  // Cosmetics for histo_hits
//   for ( auto& element : histo_hits_map_ ) { 
//     element.second->Scale( 1./element.second->GetMaximum() );
//   }
  
  //Print parameters for RecoCTPPS/TotemRPLocal/src/CTPPSDiamondTimingGetParameters.cc
//   TF1 gausForFit("gausForFit","gaus", -5, 20);
//   for ( auto& chan : time_delay_map_ ) {
//     chan.second->Fit("gausForFit","RQ");
//     std::cout << "t0_map_[ " << chan.first << " ] = " << gausForFit.GetParameter(1) << "; // Error " << gausForFit.GetParError(1) << std::endl;
//   }
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyAnalyzer);
