/****************************************************************************
 *
 * This is a part of CTPPS offline software.
 * Authors:
 *   Laurent Forthomme (laurent.forthomme@cern.ch)
 *   Nicola Minafra (nicola.minafra@cern.ch)
 *
 ****************************************************************************/

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "RecoCTPPS/TotemRPLocal/interface/CTPPSDiamondTrackRecognition.h"
#include "RecoCTPPS/TotemRPLocal/interface/CTPPSDiamondTimingCorrection.h"

class CTPPSDiamondLocalTrackFitter : public edm::stream::EDProducer<>
{
  public:
    explicit CTPPSDiamondLocalTrackFitter( const edm::ParameterSet& );
    ~CTPPSDiamondLocalTrackFitter();

    static void fillDescriptions( edm::ConfigurationDescriptions& );

  private:
    virtual void produce( edm::Event&, const edm::EventSetup& ) override;
    /// Check if one of the edges of the recHit is within the local track
    bool hitBelongsToTrack( const CTPPSDiamondLocalTrack& localTrack, const CTPPSDiamondRecHit& recHit ) const;

    edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > recHitsToken_;
    CTPPSDiamondTrackRecognition trk_algo_45_;
    CTPPSDiamondTrackRecognition trk_algo_56_;
};

CTPPSDiamondLocalTrackFitter::CTPPSDiamondLocalTrackFitter( const edm::ParameterSet& iConfig ) :
  recHitsToken_ ( consumes< edm::DetSetVector<CTPPSDiamondRecHit> >( iConfig.getParameter<edm::InputTag>( "recHitsTag" ) ) ),
  trk_algo_45_  ( iConfig.getParameter<edm::ParameterSet>( "trackingAlgorithmParams" ) ),
  trk_algo_56_  ( iConfig.getParameter<edm::ParameterSet>( "trackingAlgorithmParams" ) )
{
  produces< edm::DetSetVector<CTPPSDiamondLocalTrack> >();
}

CTPPSDiamondLocalTrackFitter::~CTPPSDiamondLocalTrackFitter()
{}

void
CTPPSDiamondLocalTrackFitter::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  std::unique_ptr< edm::DetSetVector<CTPPSDiamondLocalTrack> > pOut( new edm::DetSetVector<CTPPSDiamondLocalTrack> );

  edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > recHits;
  iEvent.getByToken( recHitsToken_, recHits );

  const CTPPSDiamondDetId id_45( 0, 1, 6, 0, 0 ), id_56( 1, 1, 6, 0, 0 );

  pOut->find_or_insert( id_45 ); // tracks in 4-5
  edm::DetSet<CTPPSDiamondLocalTrack>& tracks56 = pOut->find_or_insert( id_56 ); // tracks in 5-6

  // workaround to retrieve the detset for 4-5 without losing the reference
  edm::DetSet<CTPPSDiamondLocalTrack>& tracks45 = pOut->operator[]( id_45 );

  // feed hits to the track producers
  for ( const auto& vec : *recHits ) {
    const CTPPSDiamondDetId detid( vec.detId() );

    if ( detid.arm() == 0 ) {
      for ( const auto& hit : vec ) trk_algo_45_.addHit( hit );
    }
    else if ( detid.arm() == 1 ) {
      for ( const auto& hit : vec ) trk_algo_56_.addHit( hit );
    }
  }

  // retrieve the tracks for both arms
  trk_algo_45_.produceTracks( tracks45 );
  trk_algo_56_.produceTracks( tracks56 );
  
  // Timing corrections
  float weightedAvgDen, weightedAvgNum, weightTmp;
  int mhTmp, counterTmp;
  for ( auto& localtrack : tracks45 ) {
    weightedAvgNum = .0;
    weightedAvgDen = .0;
    mhTmp = 0;
    counterTmp = 0;
    
    for ( const auto& vec : *recHits ) {
      const CTPPSDiamondDetId detid( vec.detId() );
      if ( detid.arm() != 0 ) continue;

      for ( const auto& hit : vec ) {
        // first check if the hit contributes to the track
        if ( !hitBelongsToTrack( localtrack, hit ) ) continue;

        weightTmp = pow( 1./hit.getTPrecision(), 2 );
        weightedAvgNum += hit.getT() * weightTmp;
        weightedAvgDen += weightTmp;
        ++counterTmp;
        if ( hit.getMultipleHits() ) ++mhTmp;
      }
    }
    
    localtrack.setT( weightedAvgNum/weightedAvgDen );
    localtrack.setTSigma( std::sqrt( weightedAvgDen ) );
  }

  for ( auto& localtrack : tracks56 ) {
    weightedAvgNum = .0;
    weightedAvgDen = .0;
    mhTmp = 0;
    counterTmp = 0;
    
    for ( const auto& vec : *recHits ) {
      const CTPPSDiamondDetId detid( vec.detId() );
      if ( detid.arm() != 1 ) continue;

      for ( const auto& hit : vec ) {
        // first check if the hit contributes to the track
        if ( !hitBelongsToTrack( localtrack, hit ) ) continue;

        weightTmp = pow( 1./hit.getTPrecision(), 2 );
        weightedAvgNum += hit.getT() * weightTmp;
        weightedAvgDen += weightTmp;
        ++counterTmp;
        if ( hit.getMultipleHits() ) ++mhTmp;
      }
    }
    
    localtrack.setT( weightedAvgNum/weightedAvgDen );
    localtrack.setTSigma( std::sqrt( weightedAvgDen ) );
    localtrack.setNumOfHits( counterTmp );
    localtrack.setMultipleHits( mhTmp );
  }

  iEvent.put( std::move( pOut ) );

  // remove all hits from the track producers to prepare for the next event
  trk_algo_45_.clear();
  trk_algo_56_.clear();
}

bool
CTPPSDiamondLocalTrackFitter::hitBelongsToTrack( const CTPPSDiamondLocalTrack& localTrack, const CTPPSDiamondRecHit& recHit ) const
{
  return ( recHit.getX() + recHit.getXWidth() > localTrack.getX0() - localTrack.getX0Sigma()
        && recHit.getX() + recHit.getXWidth() < localTrack.getX0() + localTrack.getX0Sigma() )
      || ( recHit.getX() - recHit.getXWidth() > localTrack.getX0() + localTrack.getX0Sigma()
        && recHit.getX() - recHit.getXWidth() < localTrack.getX0() - localTrack.getX0Sigma() );
};

void
CTPPSDiamondLocalTrackFitter::fillDescriptions( edm::ConfigurationDescriptions& descr )
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>( "recHitsTag", edm::InputTag( "ctppsDiamondRecHits" ) )
    ->setComment( "input rechits collection to retrieve" );
  desc.add<int>( "verbosity", 0 )
    ->setComment( "general verbosity of this module" );

  edm::ParameterSetDescription trackingAlgoParams;
  trackingAlgoParams.add<double>( "threshold", 1.5 )
    ->setComment( "minimal number of rechits to be observed before launching the track recognition algorithm" );
  trackingAlgoParams.add<double>( "thresholdFromMaximum", 0.5 );
  trackingAlgoParams.add<double>( "resolution", 0.01 /* mm */ )
    ->setComment( "spatial resolution on the horizontal coordinate (in mm)" );
  trackingAlgoParams.add<double>( "sigma", 0.1 );
  trackingAlgoParams.add<double>( "startFromX", -0.5 /* mm */ )
    ->setComment( "starting horizontal coordinate of rechits for the track recognition" );
  trackingAlgoParams.add<double>( "stopAtX", 19.5 /* mm */ )
    ->setComment( "ending horizontal coordinate of rechits for the track recognition" );

  trackingAlgoParams.add<std::string>( "pixelEfficiencyFunction", "(TMath::Erf((x-[0]+0.5*[1])/([2]/4)+2)+1)*TMath::Erfc((x-[0]-0.5*[1])/([2]/4)-2)/4" )
    ->setComment( "efficiency function for single pixel\n"
                  "can be defined as:\n"
                  " * Precise: (TMath::Erf((x-[0]+0.5*[1])/([2]/4)+2)+1)*TMath::Erfc((x-[0]-0.5*[1])/([2]/4)-2)/4\n"
                  " * Fast: (x>[0]-0.5*[1])*(x<[0]+0.5*[1])+((x-[0]+0.5*[1]+[2])/[2])*(x>[0]-0.5*[1]-[2])*(x<[0]-0.5*[1])+(2-(x-[0]-0.5*[1]+[2])/[2])*(x>[0]+0.5*[1])*(x<[0]+0.5*[1]+[2])\n"
                  " * Legacy: (1/(1+exp(-(x-[0]+0.5*[1])/[2])))*(1/(1+exp((x-[0]-0.5*[1])/[2])))\n"
                  "with:\n"
                  "  [0]: centre of pad\n"
                  "  [1]: width of pad\n"
                  "  [2]: sigma: distance between efficiency ~100 -> 0 outside width" );

  trackingAlgoParams.add<double>( "yPosition", 0.0 )
    ->setComment( "vertical offset of the outcoming track centre" );
  trackingAlgoParams.add<double>( "yWidth", 0.0 )
    ->setComment( "vertical track width" );

  desc.add<edm::ParameterSetDescription>( "trackingAlgorithmParams", trackingAlgoParams )
    ->setComment( "list of parameters associated to the track recognition algorithm" );

  descr.add( "ctppsDiamondLocalTracks", desc );
}

DEFINE_FWK_MODULE( CTPPSDiamondLocalTrackFitter );
