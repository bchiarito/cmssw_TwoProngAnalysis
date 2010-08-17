import FWCore.ParameterSet.Config as cms

diphotonAnalyzer = cms.EDAnalyzer('ExoDiPhotonAnalyzer',
                                  photonCollection = cms.untracked.InputTag("photons"),
                                  ptMin = cms.untracked.double(10),
                                  hltResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                  L1Results = cms.untracked.InputTag("gtDigis"),
                                  removeSpikes = cms.untracked.bool(False),
                                  requireTightPhotons = cms.untracked.bool(True)
)
