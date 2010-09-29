import FWCore.ParameterSet.Config as cms

diphotonSignalMCAnalyzer = cms.EDAnalyzer('ExoDiPhotonSignalMCAnalyzer',
                                          photonCollection = cms.untracked.InputTag("photons"),
                                          ptMin = cms.untracked.double(10),
            # careful with the HLT process name for MC samples!
            # it changes every time there is a re-reco done on the same RAW files!
            # eg for Spring10 35X samples, I believe the name is "REDIGI"
                                     hltResults = cms.untracked.InputTag("TriggerResults","","HLT"),
)
