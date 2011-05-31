import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# just for testing!
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/405/362A2A54-D64F-E011-9AC5-003048F11CF0.root',
    'file:/mnt/tier2/store/mc/Spring11/WH_ZH_TTH_HToGG_M-120_7TeV-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0013/CE7189BA-C353-E011-A55A-E41F13181B60.root'
#'file:/mnt/tier2/store/mc/Spring11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0001/44C6CDAE-574E-E011-AA68-00237DA15F56.root'
        )
)

# need to introduce the global tag now
# because the L1GtUtils method needs to fetch records...
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#use the right global tag!
# global tag for prompt reco with 3_11X (as of Mar25th, 2011)
process.GlobalTag.globaltag = 'GR_P_V14::All'
# global tag for prompt reco with 38X (as of Sept30)
#process.GlobalTag.globaltag = 'GR_R_38X_V14::All'
# this is global tag for PromptReco with 36X
#process.GlobalTag.globaltag = 'GR10_P_V7::All'
# this is the tag claimed for data reprocessing with 35X
#process.GlobalTag.globaltag = 'GR_R_35X_V8::All'
# and this is the tag for prompt reco with 35X
#process.GlobalTag.globaltag = 'GR10_P_V5::All'

# geometry for ecal 
process.load("Configuration.StandardSequences.Geometry_cff")

# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('diphoton_tree3.root')
)

# filter on good vertex
# based on example in CMSSW/HeavyFlavorAnalysis/Onia2MuMu/test/onia2MuMuPATData_cfg.py
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
)
#process.primaryVertexPath = cms.Path(process.primaryVertexFilter)


# filter out scraping
# based on Onia example, and CMSSW/DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py for the GOODCOLL skim defn
# this requires that if there is >10 tracks,
# then at least 0.25 fraction of them must be 'high purity'

process.noScraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
)


process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)


#load diphoton analyzer
process.load("DiPhotonAnalysis.ExoDiPhotonAnalyzer.exodiphotonanalyzer_withrho_cfi")
#process.photonAnalyzer.photonCollection = "photons"

process.diphotonAnalyzer.rhoCorrection = cms.InputTag("kt6PFJets","rho")
process.diphotonAnalyzer.ptMin = 30.0 # pt cut on all photons
process.diphotonAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonAnalyzer.requireTightPhotons = False # ie only tight photons will be written 



# include all the filters as well as the analyzer
process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets+process.diphotonAnalyzer)

