import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(0) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root',',' with the source file you want to use
    fileNames = cms.untracked.vstring(
    *(
      '/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/02448DBC-0A77-E711-8240-02163E01A2ED.root',
#'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/106FE4C4-0D77-E711-8785-02163E01373C.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/469D8C89-1477-E711-A6A4-02163E01190C.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/4C729FDF-0C77-E711-A078-02163E014548.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/664E5FE8-0C77-E711-917D-02163E01469B.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/6A5C61DE-0C77-E711-A08C-02163E013720.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/6A76D2EC-0E77-E711-8903-02163E01A7A5.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/82937E4A-0C77-E711-8725-02163E01421E.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/8A8595B1-1177-E711-B1E8-02163E014761.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/A2BDBF69-1677-E711-83E2-02163E01A4ED.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/A4F9404C-0B77-E711-ADB0-02163E0134D7.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/A615970E-1177-E711-9AA8-02163E0139D9.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/ACBCD4BB-0F77-E711-9AB2-02163E014342.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/C2D92F59-0C77-E711-B3AD-02163E01A5B4.root',
'/store/data/Run2017C/ZeroBias/AOD/PromptReco-v2/000/300/088/00000/E041BD4B-0B77-E711-A90A-02163E012140.root',
  
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

process.ctppsDiamondRecHits.planeInversion = cms.int32(1)               # Put >0 for Runs < 300670     #FIXME
process.ctppsDiamondRecHits.coarseCorrection = cms.int32(1)               # Put >0 for Runs < ??       #FIXME

# rechits production
process.load('Geometry.VeryForwardGeometry.geometryRP_cfi')
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondRecHits_cfi')

# local tracks fitter
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondLocalTracks_cfi')

process.ctppsDiamondLocalTracks.trackingAlgorithmParams.threshold = cms.double(1.9)
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
     fileName = cms.string('300088_allBunches_updated3.root')                                                   #FIXME
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

