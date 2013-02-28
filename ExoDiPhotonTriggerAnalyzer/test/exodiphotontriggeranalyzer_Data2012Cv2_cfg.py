import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

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
#global tag for 2012C v2 within 532patch4
process.GlobalTag.globaltag = 'GR_P_V41_AN3::All'
## #global tag for 2012C ecal recover within 532patch4
## process.GlobalTag.globaltag = 'GR_P_V42C::All'
## #global tag for 2012D within 532patch4
## process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)

process.demo = cms.EDAnalyzer('ExoDiPhotonTriggerAnalyzer',
                                  photonCollection = cms.untracked.InputTag("photons"),
                                  ptMin = cms.untracked.double(80),
                                  hltResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                  L1Results = cms.untracked.InputTag("gtDigis"),
                                  rho25Correction = cms.InputTag("kt6PFJets25","rho"),
                                  pileupCorrection = cms.untracked.InputTag("addPileupInfo"),
                                  removeSpikes = cms.untracked.bool(False),
                                  requireTightPhotons = cms.untracked.bool(True),
                                  requireGenEventInfo = cms.untracked.bool(False),
                                  isMC = cms.untracked.bool(False),
                                  #If running on Data do not change. Only loaded as strings into the Analyzer and will have no effect. If on MC then change to the specific files you need lodaed for the MC   
                                  PUMCFileName = cms.untracked.string("PileUpMC.root"),
                                  PUDataFileName = cms.untracked.string("PileupDataAug10thHistogram.root"),
                                  PUMCHistName = cms.untracked.string("MCPileUpHistogram"),
                                  PUDataHistName = cms.untracked.string("pileup"),
                                  PFIDCategory = cms.untracked.string("Tight"),
                                  IDMethod = cms.untracked.string("ParticleFlow")
)

# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ExoDiPhotonTriggerAnalyzer.root')
)

process.p = cms.Path(process.kt6PFJets25+process.demo)
