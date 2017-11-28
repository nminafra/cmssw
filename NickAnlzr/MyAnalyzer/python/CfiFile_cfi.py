import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer("MyAnalyzer",
    #tagStatus = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
    tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
    #tagFEDInfo = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
    #tagDiamondRecHits = cms.InputTag("ctppsDiamondRecHits"),
    #tagDiamondLocalTracks = cms.InputTag("ctppsDiamondLocalTracks"),
    #tagLocalTrack = cms.InputTag("totemRPLocalTrackFitter"),
    
     
    verbosity = cms.untracked.uint32(0),
)