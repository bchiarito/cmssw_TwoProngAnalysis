from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('globalTag',
                '76X_dataRun2_16Dec2015_v0',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "global tag to use when running")
options.register('useAOD',
                True,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "whether or not to use AOD")
options.register('isMC',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "whether to run over data or MC")
options.register('pumcfilename',
                'PileUpMC_DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V7C-v1_rebinned.root',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "MC Pileup Filename")
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from a GJet PT40 dataset
##'root://eoscms//eos/cms/store/data/Commissioning2015/MinimumBias/AOD/PromptReco-v1/000/232/881/00000/04078E81-49AB-E411-A8E9-02163E012295.root'
##'root://xrootd.unl.edu//store/mc/Phys14DR/RSGravToGG_kMpl-001_M-1500_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/10000/1A742A2B-5167-E411-A342-002590A831AA.root'
#'root://eoscms//eos/cms/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/246/865/00000/0EA17D6D-B609-E511-9404-02163E014682.root'
##'root://eoscms//eos/cms/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/246/908/00000/448D1972-E909-E511-A7E0-02163E0125CE.root'
#'/tmp/charaf/448D1972-E909-E511-A7E0-02163E0125CE.root'
# 'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/GGJets_M-200To500_Pt-50_13TeV-sherpa/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/20000/109053D5-A022-E511-AD36-001E67A42BA2.root'
'root://cmsxrootd.fnal.gov//store/data/Run2015D/DoubleEG/AOD/PromptReco-v3/000/256/630/00000/2AE2E490-235F-E511-B4F3-02163E01383F.root'
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from a GJet PT40 dataset
'root://eoscms.cern.ch//store/mc/RunIIFall15MiniAODv2/GGJets_M-1000To2000_Pt-50_13TeV-sherpa/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/04B53B17-24D9-E511-B1ED-00259075D72E.root'
# "root://cmsxrootd.fnal.gov//store/data/Run2015D/DoubleEG/MINIAOD/PromptReco-v3/000/257/969/00000/F24329DE-706A-E511-998A-02163E012B1A.root"
# "file:pickevents.root"
# 'root://cmsxrootd.fnal.gov//store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/000FAE50-82A6-E511-BC87-00261894397F.root'
    )

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed

if options.useAOD :
    inputFiles = inputFilesAOD
    outputFile = "photon_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "photon_ntuple_mini.root"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", 
                fileNames      = inputFiles,
                # debugVerbosity = cms.untracked.uint32(200),
                # debugFlag      = cms.untracked.bool(True)

 ) 


# need to introduce the global tag now
# because the L1GtUtils method needs to fetch records...
# if options.globalTag == 'MCRUN2_74_V9::All':
#   process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
# else:
#   process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')


#use the right global tag!
##process.GlobalTag.globaltag = 'GR_P_V54::All'
##process.GlobalTag.globaltag = 'GR_E_V48::All'
##process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
# process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'
process.GlobalTag.globaltag = options.globalTag

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
vtxCollName = 'offlinePrimaryVertices'
if not options.useAOD:
  vtxCollName = 'offlineSlimmedPrimaryVertices'
print "Using %s vtx collection"%vtxCollName
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag(vtxCollName),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
)
#process.primaryVertexPath = cms.Path(process.primaryVertexFilter)


# filter out scraping
# based on Onia example, and CMSSW/DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py for the GOODCOLL skim defn
# this requires that if there is >10 tracks,
# then at least 0.25 fraction of them must be 'high purity'


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
if options.useAOD == True :
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

## update AK4PFchs jet collection in MiniAOD JECs

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Do not forget 'L2L3Residual' on data!
)

#load diphoton analyzer
process.load("DiPhotonAnalysis.ExoDiPhotonAnalyzer.exodiphotonanalyzer_withrho_cfi")
#process.photonAnalyzer.photonCollection = "photons"
##process.diphotonAnalyzer.rho25Correction = cms.InputTag("kt6PFJets25","rho") 
process.diphotonAnalyzer.rho25Correction = cms.InputTag("fixedGridRhoFastjetAll","rho") 
process.diphotonAnalyzer.ptMin = 50 # pt cut on all photons
process.diphotonAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonAnalyzer.requireTightPhotons = False # ie only tight photons will be written 
process.diphotonAnalyzer.requireGenEventInfo = False #write MC info when running on MC

process.diphotonAnalyzer.isAOD = cms.bool(options.useAOD) # True=AOD, False=MiniAOD
process.diphotonAnalyzer.isMC = cms.untracked.bool(options.isMC) # False by default, run with isMC=True for MC
process.diphotonAnalyzer.IDMethod = cms.untracked.string("highpt")
process.diphotonAnalyzer.PFIDCategory = cms.untracked.string("Loose")
process.diphotonAnalyzer.photonCollection = cms.untracked.InputTag("gedPhotons")

# If running on data the following four entries should not be changed. They are loaded into the analyzer as strings but in the case isMC = False then all the both old_pu_n and pu_n will both be filled with -9999.99

process.diphotonAnalyzer.PUDataFileName = 'PileupDataAug10thHistogram.root' #DataPileUp
process.diphotonAnalyzer.PUMCFileName = cms.untracked.string(options.pumcfilename)  #"MC PileUP"
process.diphotonAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonAnalyzer.PUMCHistName = "MCPileUpHisto" #Name of histogram in PUMCFileName  Need to be binned to 80

##process.path  = cms.Path(process.primaryVertexFilter+process.noScraping+process.kt6PFJets25+process.diphotonAnalyzer)
##process.path  = cms.Path(process.kt6PFJets25+process.photonIDValueMapProducer*process.diphotonAnalyzer)
process.path  = cms.Path(process.primaryVertexFilter*process.egmPhotonIDSequence*process.diphotonAnalyzer)
##process.path  = cms.Path(process.egmPhotonIDSequence*process.diphotonAnalyzer)
##process.path  = cms.Path(process.egmPhotonIDSequence)
##process.path  = cms.Path(process.diphotonAnalyzer)
##process.path  = cms.Path(process.photonIDValueMapProducer)
