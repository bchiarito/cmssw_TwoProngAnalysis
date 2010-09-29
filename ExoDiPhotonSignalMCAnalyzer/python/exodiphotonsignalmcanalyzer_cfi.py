import FWCore.ParameterSet.Config as cms

diphotonSignalMCAnalyzer = cms.EDAnalyzer('ExoDiPhotonSignalMCAnalyzer',
                                          photonCollection = cms.untracked.InputTag("photons"),
                                          ptMin = cms.untracked.double(10)
)
