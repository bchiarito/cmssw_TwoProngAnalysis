import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('pumcfilename','PileUpDiPhotonJetsBoxSherpaM500_1000_2015AOD_Asympt25ns_Pileup.root',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"MC Pileup file")
options.parseArguments()

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from a GJet PT40 dataset
    #'root://eoscms//eos/cms/store/mc/RunIISpring15DR74/GGJets_M-500To1000_Pt-50_13TeV-sherpa/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/20000/0A63F950-BC22-E511-9A17-02163E00B8D3.root'
    'root://cmsxrootd.fnal.gov///store/mc/RunIISpring15DR74/GGJets_M-1000To2000_Pt-50_13TeV-sherpa/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/40000/0416B71D-1D3A-E511-B439-C4346BB2E5F8.root'
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from a GJet PT40 dataset
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


# need to introduce the global tag now
# because the L1GtUtils method needs to fetch records...
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#use the right global tag!
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

# geometry for ecal 
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


##process.load('RecoJets.Configuration.RecoPFJets_cff')
##process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
##process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)



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

process.load("DiPhotonAnalysis.ExoDiPhotonAnalyzer.ana")
process.xsecanal = cms.EDAnalyzer("GenXSecAnalyzer")

#load diphoton analyzer
process.load("DiPhotonAnalysis.ExoDiPhotonAnalyzer.exodiphotonanalyzer_withrho_cfi")
#process.photonAnalyzer.photonCollection = "photons"
##process.diphotonAnalyzer.rho25Correction = cms.InputTag("kt6PFJets25","rho") 
process.diphotonAnalyzer.rho25Correction = cms.InputTag("fixedGridRhoFastjetAll","rho") 
process.diphotonAnalyzer.ptMin = 20 # pt cut on all photons
process.diphotonAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonAnalyzer.requireTightPhotons = False # ie only tight photons will be written 
process.diphotonAnalyzer.requireGenEventInfo = True #write MC info when running on MC

process.diphotonAnalyzer.isMC = True #MC = True or  Data = False
##process.diphotonAnalyzer.IDMethod = cms.untracked.string("highpt")
process.diphotonAnalyzer.IDMethod = cms.untracked.string("highpt")
process.diphotonAnalyzer.PFIDCategory = cms.untracked.string("Loose")
process.diphotonAnalyzer.photonCollection = cms.untracked.InputTag("gedPhotons")

# If running on data the following four entries should not be changed. They are loaded into the analyzer as strings but in the case isMC = False then all the both old_pu_n and pu_n will both be filled with -9999.99

process.diphotonAnalyzer.PUDataFileName = 'PileUpData2015DHighPtIDv2.root' #DataPileUp
##process.diphotonAnalyzer.PUMCFileName = 'PileUpDiPhotonJetsBoxSherpaM500_1000_2015AOD_Asympt25ns_Pileup.root'  #"MC PileUP"
process.diphotonAnalyzer.PUMCFileName = options.pumcfilename
process.diphotonAnalyzer.PUDataHistName = "pileupData" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonAnalyzer.PUMCHistName = "pileupMC" #Name of histogram in PUMCFileName  Need to be binned to 80

##process.diphotonAnalyzer.processBinID = options.binID

process.path  = cms.Path(process.xsecanal*process.primaryVertexFilter*process.noScraping*process.egmPhotonIDSequence*process.diphotonAnalyzer)
