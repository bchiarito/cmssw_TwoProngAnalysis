import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonSignalMCAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet( threshold = cms.untracked.string('WARNING') ),
    destinations = cms.untracked.vstring('cout')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
      # test official sample --> /store/mc/Summer12_DR53X/RSGravToGG_kMpl-001_M-1750_TuneZ2star_8TeV-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/
     # 'file:/data2/scooper/DiPhotons/RSGravToGG_kMpl-001_M-1750_Summer12_AODSIM/3AE503BF-E102-E211-ABF8-78E7D1E49636.root'
#        '
#    )
#)

inputFilesAOD = cms.untracked.vstring(
    # AOD test files
    #'root://cmsxrootd.fnal.gov///store/mc/RunIISpring15DR74/RSGravToGG_kMpl-01_M-1250_TuneCUEP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/18D3BCE6-6105-E511-8177-02163E010D77.root'
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files
    # none
    )

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = True

if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "photon_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "photon_ntuple_mini.root"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles ) 


# global tag for MC because now we need geometry
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All' #50 ns 
#'MCRUN2_74_V9::All' 25 ns

# geometry for ecal 
process.load("Configuration.Geometry.GeometryIdeal_cff")

# load jets.
#process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
#process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)


# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ExoDiPhotonSignalAnalyzer.root')
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



##-----------------taken from Ilya-----------------
#
# Several photon variables can not be found inside of a photon object
# and it is easiest to compute them upstream with a dedicated producer,
# such as this standard producer used for photon ID.
#    The producer computes full5x5 cluster shapes and PF isolation variables.
#
# Do not forget to add this producer to the path below!
#
##process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
##-----------------taken from Ilya-----------------
##process.photonIDValueMapProducer.src = cms.InputTag('gedPhotons')
##process.photonIDValueMapProducer.particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons")

#
# Set up photon ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
##-----------------taken from Ilya-----------------



#process.diphotonAnalysis = cms.EDAnalyzer('ExoDiPhotonSignalMCAnalyzer')
process.load("DiPhotonAnalysis/ExoDiPhotonSignalMCAnalyzer/exodiphotonsignalmcanalyzer_cfi")
process.diphotonSignalMCAnalyzer.ptMin = 0.0
#process.diphotonSignalMCAnalyzer.rho25Correction = cms.InputTag("kt6PFJets25","rho")
process.diphotonSignalMCAnalyzer.rho25Correction = cms.InputTag("fixedGridRhoFastjetAll","rho") 
process.diphotonSignalMCAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonSignalMCAnalyzer.requireTightPhotons = False # ie only tight photons will be written
process.diphotonSignalMCAnalyzer.requireGenEventInfo = True #write MC info when running on MC
process.diphotonSignalMCAnalyzer.PUDataFileName = 'PileupDataDec14thHistogram.root'
process.diphotonSignalMCAnalyzer.PUMCFileName = 'PileUpMC.root'
process.diphotonSignalMCAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonSignalMCAnalyzer.PUMCHistName = "MCPileUpHisto"  #Name of histogram in PUMCFileName  Need to be binned to 80

#process.diphotonSignalMCAnalyzer.isMC = False #MC = True or  Data = False
process.diphotonSignalMCAnalyzer.IDMethod = cms.untracked.string("ParticleFlow")
process.diphotonSignalMCAnalyzer.PFIDCategory = cms.untracked.string("Loose")
process.diphotonSignalMCAnalyzer.photonCollection = cms.untracked.InputTag("gedPhotons")


#process.p = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets25+process.diphotonSignalMCAnalyzer)

process.p = cms.Path(process.primaryVertexFilter+process.noScraping+process.egmPhotonIDSequence+process.diphotonSignalMCAnalyzer)
