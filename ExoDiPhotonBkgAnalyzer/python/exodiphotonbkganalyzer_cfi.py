import FWCore.ParameterSet.Config as cms

diphotonBkgAnalyzer = cms.EDAnalyzer('ExoDiPhotonBkgAnalyzer',
                                      photonCollection = cms.untracked.InputTag("photons"),
                                      ptMin = cms.untracked.double(10)
)
