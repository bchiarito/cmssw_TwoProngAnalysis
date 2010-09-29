import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonSignalMCAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
    'file:RSGravToGG_kMpl001_M-1000_7TeV-pythia6_Spring10-START3X_V26-v1_GEN-SIM-RECO_50events.root'
    )
)

# global tag for MC because now we need geometry
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START3X_V27::All'

# geometry for ecal 
process.load("Configuration.StandardSequences.Geometry_cff")



# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('diphoton_tree.root')
)


#process.diphotonAnalysis = cms.EDAnalyzer('ExoDiPhotonSignalMCAnalyzer'
#)
process.load("DiPhotonAnalysis/ExoDiPhotonSignalMCAnalyzer/exodiphotonsignalmcanalyzer_cfi");
process.diphotonSignalMCAnalyzer.ptMin = 0.0 


process.p = cms.Path(process.diphotonSignalMCAnalyzer)
