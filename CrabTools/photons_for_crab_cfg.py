#config file for diphoton physics plot for Exotica October Exercise
# this version configured for crab use
# Conor Henderson, 5 October 2009

#Import configuration python objects
import FWCore.ParameterSet.Config as cms

#there must be exactly one 'cms.Process', called 'process'
process = cms.Process("ConorsAnalysis")

#input file to process - for crab, this should be empty
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#     'file:/afs/cern.ch/user/c/chenders/scratch0/CMSSW_3_1_1/src/relval_gammajets_80_120.root'
    )
)

# set the number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photon_hists.root')
)

#load photon analyzer
process.load("DiPhotonAnalysis.PhotonAnalyzer.photonanalyzer_cfi")
process.photonAnalyzer.photonCollection = "photons"
process.photonAnalyzer.ptMin = 10.0;

process.path  = cms.Path(process.photonAnalyzer)

# dont need output module unless I want output events


