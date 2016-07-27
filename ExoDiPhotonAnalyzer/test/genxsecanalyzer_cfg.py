import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
                    input = cms.untracked.int32(10000)
)
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#'file:/uscms_data/d3/skaplan/JetFlavorStudies/CMSSW_5_3_13/src/Analyzers/JetMatchingAnalyzer/chadronghost_lighthadrondr_events.root')
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall15DR76/GGJets_M-500To1000_Pt-50_13TeV-sherpa/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/46954485-C3D7-E511-99E1-B083FECFF297.root')
)
process.printEventNumber = cms.OutputModule("AsciiOutputModule")

# process.printList = cms.EDAnalyzer("ParticleListDrawer",
#                     src = cms.InputTag("genParticles"),
#                     maxEventsToPrint  = cms.untracked.int32(-1)
# )
process.g = cms.EDAnalyzer("GenXSecAnalyzer")

process.p = cms.Path(process.g)

# process.outpath = cms.EndPath(process.printEventNumber)
process.MessageLogger.destinations = cms.untracked.vstring('cout','cerr')
