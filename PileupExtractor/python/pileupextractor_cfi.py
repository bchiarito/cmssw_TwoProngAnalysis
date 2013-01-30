import FWCore.ParameterSet.Config as cms

diphotonPileUpExtractor = cms.EDAnalyzer('PileupExtractor',
                                  pileupCollection = cms.untracked.InputTag("addPileupInfo"),
)
