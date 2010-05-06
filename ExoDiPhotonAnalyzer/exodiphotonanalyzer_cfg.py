import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root',
    'file:/tmp/chenders/02071949-FA40-DF11-9990-001A64789DEC.root',
    )
)

# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('diphoton_tree.root')
)

#load diphoton analyzer
process.load("DiPhotonAnalysis.ExoDiPhotonAnalyzer.exodiphotonanalyzer_cfi")
#process.photonAnalyzer.photonCollection = "photons"
process.diphotonAnalyzer.ptMin = 0.0 # jsut to get some entries

process.path  = cms.Path(process.diphotonAnalyzer)

