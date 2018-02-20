import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(0) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root',',' with the source file you want to use
    fileNames = cms.untracked.vstring(
    *(
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/1E2FB550-5190-E711-8010-02163E01A731.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/ECE07B13-5C90-E711-AA7C-02163E011DDC.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/AABBE824-5990-E711-BE37-02163E0144E2.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/562C1C92-5A90-E711-904D-02163E01A2C7.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/D23D5509-6390-E711-B82E-02163E014785.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/88377025-5990-E711-8E26-02163E01A4E6.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/FE7699DC-5990-E711-8DF5-02163E011D9E.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/600BB493-5A90-E711-84EE-02163E01A1EB.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/AA6E5A99-6090-E711-B9AA-02163E012B3A.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/F45AD155-5B90-E711-8ACE-02163E019C49.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/909B3BB3-5C90-E711-8649-02163E01A1EB.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/ACC8407D-5890-E711-A88A-02163E0134FB.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/A4DF7114-5E90-E711-9BD7-02163E011D9E.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/9EF2EBDC-5990-E711-BC87-02163E0138EC.root',
'/store/data/Run2017D/ZeroBias/AOD/PromptReco-v1/000/302/159/00000/DE7748D6-5C90-E711-A008-02163E011B3C.root',
  
     )
    ),
    #inputCommands = cms.untracked.vstring(
       #'drop CTPPSPixelCluseredmDetSetVector_ctppsPixelClusters__RECO'
    #)
)
    
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_relval', '')

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

#process.ctppsDiamondRecHits.planeInversion = cms.int32(1)               # Put >0 for Runs < 300670     #FIXME
#process.ctppsDiamondRecHits.coarseCorrection = cms.int32(1)               # Put >0 for Runs < ??       #FIXME

# rechits production
process.load('Geometry.VeryForwardGeometry.geometryRP_cfi')
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondRecHits_cfi')

# local tracks fitter
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondLocalTracks_cfi')

process.ctppsDiamondLocalTracks.trackingAlgorithmParams.threshold = cms.double(1.5)
process.ctppsDiamondLocalTracks.trackingAlgorithmParams.sigma = cms.double(0.1)
process.ctppsDiamondLocalTracks.trackingAlgorithmParams.resolution = cms.double(0.01) # in mm
#process.ctppsDiamondLocalTracks.trackingAlgorithmParams.pixelEfficiencyFunction = cms.string("(TMath::Erf((x-[0]+0.5*[1])/([2]/4)+2)+1)*TMath::Erfc((x-[0]-0.5*[1])/([2]/4)-2)/4")

# pixel
process.load('RecoCTPPS.PixelLocal.ctppsPixelLocalTracks_cfi')

#process.load("MyAnalyzer.CfiFile_cfi.py")
process.demo = cms.EDAnalyzer('MyAnalyzer',
 tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
 tagDiamondRecHits = cms.InputTag("ctppsDiamondRecHits"),
 tagDiamondLocalTracks = cms.InputTag("ctppsDiamondLocalTracks"),
 tagLocalTrack = cms.InputTag("totemRPLocalTrackFitter"),
 ctppsPixelLocalTracks = cms.InputTag("ctppsPixelLocalTracks"),
 selectedOOTIndex = cms.int32(4)                                                                        #FIXME
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('302159_15_01_001.root')                                                   #FIXME
)

process.p = cms.Path(
    #process.ctppsRawToDigi *
    process.recoCTPPS *
    #process.ctppsDiamondRawToDigi *
    process.ctppsDiamondRecHits *
    process.ctppsDiamondLocalTracks *
    process.ctppsPixelLocalTracks *
    process.demo
    )

