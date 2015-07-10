import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'root://eoscms//eos/cms/store/data/Run2012C/DoublePhotonHighPt/AOD/PromptReco-v2/000/201/602/CAAA8533-D1EF-E111-87A1-0025901D5C80.root'
     )
)

# need to introduce the global tag now
# because the L1GtUtils method needs to fetch records...
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#use the right global tag!
## #global tag for 2012 A+B rereco within 532patch4
## process.GlobalTag.globaltag = 'FT_53_V6_AN3::All'
## #global tag for 2012 A+B ecal recover within 532patch4
## process.GlobalTag.globaltag = 'FT_53_V6C_AN3::All'
## #global tag for 2012C v1 within 532patch4
## process.GlobalTag.globaltag = 'FT_53_V10_AN3::All'
## #global tag for 2012C v2 within 532patch4
## process.GlobalTag.globaltag = 'GR_P_V41_AN3::All'
## #global tag for 2012C ecal recover within 532patch4
## process.GlobalTag.globaltag = 'GR_P_V42C::All'
#global tag for 2012D within 532patch4
process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'

# geometry for ecal 
#When in 5_3_X Need to use diff GeometryDB
##process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")

# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ExoDiPhotonAnalyzer.root')
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
process.diphotonAnalyzer.isMC = False #MC = True or  Data = False
process.diphotonAnalyzer.IDMethod = cms.untracked.string("ParticleFlow")
process.diphotonAnalyzer.PFIDCategory = cms.untracked.string("Loose")

# If running on data the following four entries should not be changed. They are loaded into the analyzer as strings but in the case isMC = False then all the both old_pu_n and pu_n will both be filled with -9999.99

process.diphotonAnalyzer.PUDataFileName = 'PileupDataAug10thHistogram.root' #DataPileUp
process.diphotonAnalyzer.PUMCFileName = 'PileUpMC.root'  #"MC PileUP"
process.diphotonAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonAnalyzer.PUMCHistName = "pu_n_BeforeCuts" #Name of histogram in PUMCFileName  Need to be binned to 80
#precede with a diphoton filter, to speed things up
#process.load("DiPhotonAnalysis.DiPhotonFilter.diphotonfilter_cfi")

#process.diphotonFilter.ptMin_photon1 = 20.0 
#process.diphotonFilter.ptMin_photon2 = 20.0

# include all the filters as well as the analyzer
#process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.diphotonFilter+process.diphotonAnalyzer)
#process.path =cms.Path(process.diphotonAnalyzer)
process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets25+process.diphotonAnalyzer)
#process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets25+process.diphotonAnalyzer)
