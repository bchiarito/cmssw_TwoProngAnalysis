import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# just for testing!
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:HighMassSkims_April27_Mgreater500.root',
#    'file:/tmp/chenders/02071949-FA40-DF11-9990-001A64789DEC.root',
#'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/torimoto/0065C919-F53B-DF11-8BF5-001D09F29146.root'
        #direct from castor
#        '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/601/02071949-FA40-DF11-9990-001A64789DEC.root'
#    '/store/data/Run2010A/EG/RECO/v4/000/142/040/FE902386-B79C-DF11-88BF-00304879FC6C.root'
    # a 383 prompt reco file on castor
#    '/store/data//Run2010B//Photon/RECO/PromptReco-v2/000/146/944/00056B4C-7CCC-DF11-BAF2-00304879FA4C.root'
#    '/store/data//Run2010B/Photon/RECO/PromptReco-v2/000/148/819/02B8C9B5-45E0-DF11-8310-001D09F295A1.root',
#    '/store/data//Run2010B/Photon/RECO/PromptReco-v2/000/148/819/0636566B-59E0-DF11-A0C2-0019B9F70607.root'
#    'rfio:/castor/cern.ch/user/c/chenders/DiPhotonBorn_Pt25to250_Spring10_GEN-SIM-RECO_200events.root',
#    'rfio:/castor/cern.ch/user/y/yma/RSGravitons/Nov2010/RShighest.root',
#    '/store/data/Run2010B/Photon/RECO/PromptReco-v2/000/148/953/642845D7-C0E1-DF11-B968-0030487A18F2.root'
#  '/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/684/F652691B-0383-E111-9821-001D09F241F0.root'
#    'rfio:/castor/cern.ch/cms/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/738/E4AD6153-0684-E111-AAE3-001D09F24DA8.root'
     'rfio:/castor/cern.ch/user/j/jcarson/DiPhotonTrees/Data/HighMassSkims_April27_Mgreater500.root'      
#      'HighMassSkims_April27_Mgreater500.root'
     )
)

# need to introduce the global tag now
# because the L1GtUtils method needs to fetch records...
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#use the right global tag!
#global tag fro >382 data (Nov22)
process.GlobalTag.globaltag = 'GR_R_52_V7::All'
# global tag for prompt reco with 38X (as of Sept30)
#process.GlobalTag.globaltag = 'GR10_P_V9::All'
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
    fileName = cms.string('April27HighMassSkimsRhoStudy.root')
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


process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)


#load diphoton analyzer
process.load("DiPhotonAnalysis.ExoDiPhotonAnalyzer.exodiphotonanalyzer_withrho_cfi")
#process.photonAnalyzer.photonCollection = "photons"
process.diphotonAnalyzer.rho25Correction = cms.InputTag("kt6PFJets25","rho") 
process.diphotonAnalyzer.ptMin = 70 # pt cut on all photons
process.diphotonAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonAnalyzer.requireTightPhotons = False # ie only tight photons will be written 


# precede with a diphoton filter, to speed things up
#process.load("DiPhotonAnalysis.DiPhotonFilter.diphotonfilter_cfi")

#process.diphotonFilter.ptMin_photon1 = 20.0 
#process.diphotonFilter.ptMin_photon2 = 20.0 

# include all the filters as well as the analyzer
#process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.diphotonFilter+process.diphotonAnalyzer)
#process.path =cms.Path(process.diphotonAnalyzer)
process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets25+process.diphotonAnalyzer)
#process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets25+process.diphotonAnalyzer)
