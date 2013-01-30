import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonSignalMCAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet( threshold = cms.untracked.string('ERROR') ),
    destinations = cms.untracked.vstring('cout')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'root://eoscms//eos/cms/store/user/charaf/RSGravToGG_kMpl-001_M-1000_TuneZ2star_8TeV-pythia6/EXOMCRECO_Summer12_DR53X_PU_S10_START53_V7A-v0/06e6bebb6525c1c46ccfcc56d82513c0/RSGravToGG_kMpl-001_M-1000_TuneZ2star_8TeV-pythia6_Summer12_DR53X_PU_S10_START53_V7A-v0_10_1_BrT.root'
    )
)

# global tag for MC because now we need geometry
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#use the right global tag!
#global tag for 2012 A+B rereco within 532patch4
process.GlobalTag.globaltag = 'FT_53_V6_AN1::All'
## #global tag for 2012C v1+v2 within 532patch4
## process.GlobalTag.globaltag = 'GR_P_V40_AN1::All'
##process.GlobalTag.globaltag = 'START53_V11::All
##process.GlobalTag.globaltag = 'START53_V7A::All'

# geometry for ecal 
#outdated as of 5_3_X
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")

# load jets.
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)


# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ExoDiPhotonAnalyzerNoSel.root')
)

# filter on good vertex
# based on example in CMSSW/HeavyFlavorAnalysis/Onia2MuMu/test/onia2MuMuPATData\
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )
#process.primaryVertexPath = cms.Path(process.primaryVertexFilter)

process.noScraping = cms.EDFilter("FilterOutScraping",
applyfilter = cms.untracked.bool(True),
debugOn = cms.untracked.bool(False),
numtrack = cms.untracked.uint32(10),
thresh = cms.untracked.double(0.25)
)


#process.diphotonAnalysis = cms.EDAnalyzer('ExoDiPhotonSignalMCAnalyzer'
#)
process.load("DiPhotonAnalysis/ExoDiPhotonSignalMCAnalyzer/exodiphotonsignalmcanalyzer_cfi");
process.diphotonSignalMCAnalyzer.ptMin = 0.0 
process.diphotonSignalMCAnalyzer.rho25Correction = cms.InputTag("kt6PFJets25","rho")
process.diphotonSignalMCAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonSignalMCAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonSignalMCAnalyzer.requireTightPhotons = False # ie only tight photons will be written
process.diphotonSignalMCAnalyzer.PUDataFileName = 'PileupDataAug10thHistogram.root'
process.diphotonSignalMCAnalyzer.PUMCFileName = 'PileUpMC.root'
process.diphotonSignalMCAnalyzer.isMC = False # MC = True or  Data = False
process.diphotonSignalMCAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonSignalMCAnalyzer.PUMCHistName = "MCPileUpHisto"  #Name of histogram in PUMCFileName  Need to be binned to 80


process.p = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets25+process.diphotonSignalMCAnalyzer)
