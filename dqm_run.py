import FWCore.ParameterSet.Config as cms
import string

process = cms.Process('RECODQM')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"

process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
      FILELIST
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_relval', '')

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

# load local geometry to avoid GT
process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2018_cfi')
process.load('RecoCTPPS.TotemRPLocal.totemTimingLocalReconstruction_cff')

# CTPPS DQM modules
process.load("DQM.CTPPS.ctppsDQM_cff")

process.path = cms.Path(
    process.ctppsRawToDigi *
    process.recoCTPPS *
    process.ctppsDQM
)

process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver
)

# output configuration
# process.output = cms.OutputModule("PoolOutputModule",
#   fileName = cms.untracked.string("file:./AOD.root"),
#   outputCommands = cms.untracked.vstring(
#     'drop *',
#     'keep *_*_TotemTiming_*',
#   )
# )

# process.outpath = cms.EndPath(process.output)

process.schedule = cms.Schedule(
process.path,
process.end_path,
# process.outpath
)
