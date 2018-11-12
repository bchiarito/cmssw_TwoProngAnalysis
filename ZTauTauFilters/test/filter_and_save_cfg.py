import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
  "file:file:/cms/chiarito/samples/miniaod/dysig/rootfile_DYJetsToLL_M-50_80k.root"
  ))

process.filt = cms.EDFilter('ZtoTauHadRecoSelector',
  dumpCutflow = cms.untracked.bool(True)
  )
process.p = cms.Path(process.filt)

process.out = cms.OutputModule("PoolOutputModule", 
  fileName = cms.untracked.string("filtered_events.root"),
  SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring("p") )
  )

process.e = cms.EndPath(process.out)
