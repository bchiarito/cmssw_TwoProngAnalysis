import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(1);

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.source = cms.Source("PoolSource",
fileNames=cms.untracked.vstring(
'dcap://pnfs/cms/WAX/resilient/duong/Output_313/Sample/Pat/QcdPt30_EXODipho_7Tev_1.root',
) ,
skipEvents = cms.untracked.uint32(0) 
)

process.demo = cms.EDAnalyzer('ADDrootTreeMaker',
		outFileName = cms.string("file:ADDrootTreeMaker.root")
)

process.p = cms.Path(process.demo)
