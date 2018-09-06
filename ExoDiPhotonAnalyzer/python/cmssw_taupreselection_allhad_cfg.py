import FWCore.ParameterSet.Config as cms

Taufilter = cms.EDFilter('ZtoTauHadTruthSelector',
  filterByTruthDecayType = cms.untracked.vdouble(5.3,5.4,5.2,5.1),
  )

out = cms.OutputModule("PoolOutputModule", 
  fileName = cms.untracked.string("MINIAOD.root"),
  SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring("p") )
  )
