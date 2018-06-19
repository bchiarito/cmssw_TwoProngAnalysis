import FWCore.ParameterSet.Config as cms

process = cms.Process("GENXSANA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000
process.MessageLogger.cerr.FwkReport.limit = 10

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/00BD0FD2-FCB2-E611-BA46-0025905A609A.root')
)

process.genxs = cms.EDAnalyzer("GenXSecAnalyzer")

process.path = cms.Path(process.genxs)
