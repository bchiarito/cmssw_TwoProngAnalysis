import FWCore.ParameterSet.Config as cms

process = cms.Process("TAUFILT")

process.source = cms.Source ("PoolSource",
  fileNames = cms.untracked.vstring (
  "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00099D43-77ED-E611-8889-5065F381E1A1.root"
  ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (10000))
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.filt = cms.EDFilter('ZtoTauHad',
  tree = cms.untracked.bool(False),
  filterByTruthDecayType = cms.untracked.vdouble(5.4,5.3),
  ptMin = cms.untracked.double(40.0),
  absEtaMax = cms.untracked.double(2.5),
  )

process.out = cms.OutputModule("PoolOutputModule", 
  fileName = cms.untracked.string("MINIAOD.root"),
  SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring("p") )
  )

process.p = cms.Path(process.filt)

process.end = cms.EndPath(process.out)
