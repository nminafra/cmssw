import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('VertexDistr',
 pfTag = cms.InputTag('particleFlow'),
)
