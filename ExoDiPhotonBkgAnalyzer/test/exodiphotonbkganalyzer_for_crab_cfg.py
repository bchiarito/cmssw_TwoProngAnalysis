import FWCore.ParameterSet.Config as cms

process = cms.Process("DiPhotonBkgAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    'file:/tmp/chenders/DiPhotonBorn_Pt25to250_Spring10_GEN-SIM-RECO_200events.root'
#    'file:/tmp/chenders/PhotonJet_Pt50to80_Spring10_GEN-SIM_RECO_200events.root'
#    'file:PhotonJet_Pt50to80_Spring10_GEN-SIM_RECO_200events.root'
#    'file:DiPhotonBorn_Pt25to250_Spring10_GEN-SIM-RECO_200events.root'
    )
)

process.load("Configuration.StandardSequences.Geometry_cff")

# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('diphoton_tree.root')
)


process.load("DiPhotonAnalysis/ExoDiPhotonBkgAnalyzer/exodiphotonbkganalyzer_cfi")
process.diphotonBkgAnalyzer.ptMin = 0.0
# important to use correct HLT process name for MC samples which are re-reco's of
# existing raw files (the HLT process cant be re-run on the same file...)
# an example Spring10 35X sample contains the following
# edmDumpEventContent /tmp/chenders/DiPhotonBorn_Pt25to250_Spring10_GEN-SIM-RECO_200events.root | grep TriggerResults
#edm::TriggerResults               "TriggerResults"        ""            "HLT8E29."
#edm::TriggerResults               "TriggerResults"        ""            "HLT." 
#edm::TriggerResults               "TriggerResults"        ""            "REDIGI."
#edm::TriggerResults               "TriggerResults"        ""            "RECO."

# so, it's REDIGI in this case
process.diphotonBkgAnalyzer.hltResults = cms.untracked.InputTag("TriggerResults","","REDIGI")



# no longer want to use diphoton filter for bkg,
# since the analyser needs to run on every event for proper normalisation purposes


process.p = cms.Path(process.diphotonBkgAnalyzer)
