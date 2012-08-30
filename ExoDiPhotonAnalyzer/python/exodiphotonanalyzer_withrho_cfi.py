import FWCore.ParameterSet.Config as cms

diphotonAnalyzer = cms.EDAnalyzer('ExoDiPhotonAnalyzer',
                                  photonCollection = cms.untracked.InputTag("photons"),
                                  ptMin = cms.untracked.double(10),
                                  hltResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                  L1Results = cms.untracked.InputTag("gtDigis"),
                                  rho25Correction = cms.InputTag("kt6PFJets25","rho"),
                                  pileupCorrection = cms.untracked.InputTag("addPileupInfo"),
                                  removeSpikes = cms.untracked.bool(False),
                                  requireTightPhotons = cms.untracked.bool(True),
                                  isMC = cms.untracked.bool(True),
                                  PUMCFileName = cms.untracked.string("PileUpMC.root"),
                                  PUDataFileName = cms.untracked.string("PileupDataAug10thHistogram.root"),
                                  PUMCHistName = cms.untracked.string("MCPileUpHistogram"),
                                  PUDataHistName = cms.untracked.string("pileup")
                                  )
