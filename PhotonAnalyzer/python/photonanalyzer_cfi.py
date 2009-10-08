import FWCore.ParameterSet.Config as cms

photonAnalyzer = cms.EDAnalyzer('PhotonAnalyzer',
                                photonCollection = cms.untracked.InputTag("ThesePhotonsDoNotExist"),
                                ptMin = cms.untracked.double(10)
)
