import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:ALCARECOTkAlCosmics0T.root')
)

process.MessageLogger = cms.Service("MessageLogger")

process.demo = cms.EDFilter("Weber",
#                            src = cms.InputTag("generalTracks")
#                            src = cms.InputTag("ALCARECOTkAlCosmicsCTF0T")
#                            src = cms.InputTag("ALCARECOTkAlCosmicsCosmicTF0T")
                            src = cms.InputTag("ALCARECOTkAlCosmicsRS0T")
#                            src = cms.InputTag("generalTracks")
		
                            )

process.p = cms.Path(process.demo)


