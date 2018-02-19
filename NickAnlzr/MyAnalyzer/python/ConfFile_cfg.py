import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(0) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      *(
'/store/data/Run2017E/ZeroBias/RAW/v1/000/304/739/00000/A4A8636E-9AAD-E711-BE7D-02163E013645.root',
'/store/data/Run2017E/ZeroBias/RAW/v1/000/304/739/00000/C0CE5580-9AAD-E711-8DC8-02163E01A554.root',
    )
      )
)
    
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_relval', '')

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

# rechits production
process.load('Geometry.VeryForwardGeometry.geometryRP_cfi')
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondRecHits_cfi')

# local tracks fitter
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondLocalTracks_cfi')

#process.load("MyAnalyzer.CfiFile_cfi.py")
process.demo = cms.EDAnalyzer('MyAnalyzer',
 tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
 tagDiamondRecHits = cms.InputTag("ctppsDiamondRecHits"),
 tagDiamondLocalTracks = cms.InputTag("ctppsDiamondLocalTracks"),
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('304739.root')
)

process.p = cms.Path(
    process.ctppsRawToDigi *
    process.recoCTPPS *
    process.ctppsDiamondRawToDigi *
    process.ctppsDiamondRecHits *
    process.ctppsDiamondLocalTracks *
    process.demo
    )

