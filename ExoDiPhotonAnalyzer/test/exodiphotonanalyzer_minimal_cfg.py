# Command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('globalTag',
                '76X_dataRun2_16Dec2015_v0',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "global tag to use when running")
options.register('doLumis',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "let config file do lumi selection instead of CRAB - must be FALSE if using CRAB!")
options.register('debug',
                True,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "True includes all output, False removes most of the per event output")
options.register('fakeRateOnly',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Only fill the fake rate histos for Brandon's analysis, no tree")
options.register('pumcfilename',
                'PileUpMC_DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V7C-v1_rebinned.root',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "MC Pileup Filename")
options.register("sample",
                "",
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "which sample we want to run over")
options.setDefault('maxEvents', 1)
options.parseArguments()

# Begin configuration
import FWCore.ParameterSet.Config as cms
process = cms.Process("ExoDiPhotonAnalysis")

# Log messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# The source
sample = options.sample
readFiles = []
doLumis = options.doLumis
if sample == "local":
    isMC = True
    isSignal = True
    readFiles.extend( [
        'file:./MiniAODv2_Eta_generic.root' ] )
if sample == "jet":
    isMC = False
    isSignal = False
    readFiles.extend( [
       '/store/data/Run2016G/JetHT/MINIAOD/03Feb2017-v1/100000/006E7AF2-AEEC-E611-A88D-7845C4FC3B00.root' ] )
if sample == "photon":
    isMC = False
    isSignal = False
    readFiles.extend( [
       '/store/data/Run2016G/SinglePhoton/MINIAOD/03Feb2017-v1/110000/00F4619E-9BEB-E611-8D7A-002590494BE2.root' ] )
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring( readFiles ))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

# use JSON file to choose what lumis to process
if isMC==False and doLumis==True: 
    import FWCore.PythonUtilities.LumiList as LumiList
    goodlumis = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    process.source.lumisToProcess = LumiList.LumiList(filename = goodlumis).getVLuminosityBlockRange()

# Allow unscheduled
process.options.allowUnscheduled = cms.untracked.bool(True)

# Global Tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = options.globalTag

# Geometry for ecal 
process.load("Configuration.StandardSequences.GeometryDB_cff")

# Output file service
process.TFileService = cms.Service( "TFileService", fileName = cms.string('ExoDiPhotonAnalyzer.root') )

# filter on good vertex
# based on example in CMSSW/HeavyFlavorAnalysis/Onia2MuMu/test/onia2MuMuPATData_cfg.py
vtxCollName = 'offlineSlimmedPrimaryVertices'
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag(vtxCollName),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

# the ntuplizer
process.diphotonAnalyzer = cms.EDAnalyzer('ExoDiPhotonAnalyzer',
                                  photonCollection = cms.untracked.InputTag("photons"),
                                  ptMin = cms.untracked.double(10),
                                  hltResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                  L1Results = cms.untracked.InputTag("gtDigis"),
                                  rho25Correction = cms.InputTag("kt6PFJets25","rho"),
                                  pileupCorrection = cms.untracked.InputTag("addPileupInfo"),
                                  removeSpikes = cms.untracked.bool(False),
                                  requireTightPhotons = cms.untracked.bool(True),
                                  requireGenEventInfo = cms.untracked.bool(False),
                                  isMC = cms.untracked.bool(True),
                                  # jet cuts
                                  jetPtCut = cms.untracked.double(30.),
                                  jetEtaCut = cms.untracked.double(2.5),
                                  # charged decay configuration
                                  chargedHadronPairMinDeltaR = cms.untracked.double(0.05),
                                  chargedHadronMinPt = cms.untracked.double(10.0),
                                  isolationConeR = cms.untracked.double(0.3),
                                  photonPhiBoxSize = cms.untracked.double(0.8),
                                  photonEtaBoxSize = cms.untracked.double(0.087),
                                  photonPtCut = cms.untracked.double(10.0),
                                  chargedIsoCut = cms.untracked.double(0.1),
                                  chargedIsoFakeMax = cms.untracked.double(0.3),
                                  neutralIsoCut = cms.untracked.double(0.1),
                                  neutralIsoFakeMax = cms.untracked.double(0.3),
                                  egammaIsoCut = cms.untracked.double(0.1),
                                  egammaIsoFakeMax = cms.untracked.double(0.3),
                                  generatorEtaMatchDR = cms.untracked.double(0.1),
                                  chargedDecayCutflow = cms.untracked.bool(True),
                                  noTreeOnlyFakeRateHistos = cms.untracked.bool(False),
                                  # charged decay configuration end
                                  #If running on Data do not change. Only loaded as strings into the Analyzer and will have no effect. If on MC then change
                                  PUMCFileName = cms.untracked.string("PileUpMC.root"),
                                  PUDataFileName = cms.untracked.string("PileupDataAug10thHistogram.root"),
                                  PUMCHistName = cms.untracked.string("MCPileUpHistogram"),
                                  PUDataHistName = cms.untracked.string("pileup"),
                                  PFIDCategory = cms.untracked.string("Medium"),
                                  IDMethod = cms.untracked.string("Detector"),
                                  #Input taken from Ilya's code
                                  rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                  # Objects specific to AOD format
                                  photons = cms.InputTag("gedPhotons"),
                                  genParticles = cms.InputTag("genParticles"),
                                  # Objects specific to MiniAOD format
                                  photonsMiniAOD = cms.InputTag("slimmedPhotons"),
                                  genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                  # ValueMap names from the producer upstream
                                  full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                  phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                  phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                  phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                                  # Locations of files with the effective area constants.
                                  # The constants in these files below are derived for PHYS14 MC.
                                  effAreaChHadFile = cms.FileInPath
                                  ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
                                  effAreaNeuHadFile= cms.FileInPath
                                  ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
                                  effAreaPhoFile   = cms.FileInPath
                                  ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt"),
                                  # ID decisions (common to all formats)
                                  phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                                  phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
                                  phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
                                  # recHitsEB = cms.InputTag("reducedEcalRecHitsEB"),
                                  # recHitsEE = cms.InputTag("reducedEcalRecHitsEE"),
                                  # recHitsEBMiniAOD = cms.InputTag("reducedEgamma:reducedEBRecHits"),
                                  # recHitsEEMiniAOD = cms.InputTag("reducedEgamma:reducedEERecHits")
                                  )
process.diphotonAnalyzer.rho25Correction = cms.InputTag("fixedGridRhoFastjetAll","rho") 
process.diphotonAnalyzer.ptMin = 50 # pt cut on all photons
process.diphotonAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonAnalyzer.requireTightPhotons = False # ie only tight photons will be written 
process.diphotonAnalyzer.requireGenEventInfo = False #write MC info when running on MC
process.diphotonAnalyzer.isAOD = cms.bool(False) # True=AOD, False=MiniAOD
process.diphotonAnalyzer.isMC = cms.untracked.bool(isMC)
process.diphotonAnalyzer.isSignal = cms.untracked.bool(isSignal)
process.diphotonAnalyzer.debug = cms.untracked.bool(options.debug)
process.diphotonAnalyzer.noTreeOnlyFakeRateHistos = cms.untracked.bool(options.fakeRateOnly)
process.diphotonAnalyzer.IDMethod = cms.untracked.string("highpt")
process.diphotonAnalyzer.PFIDCategory = cms.untracked.string("Loose")
process.diphotonAnalyzer.photonCollection = cms.untracked.InputTag("gedPhotons")
process.diphotonAnalyzer.jetCollection    = cms.InputTag("slimmedJets")
process.diphotonAnalyzer.omitChargedDecayCode = cms.bool(False)
process.diphotonAnalyzer.skipPhotonMCCode = cms.bool(True)
# If running on data the following four entries should not be changed. They are loaded into the analyzer as strings but in the case isMC = False then all the both old_pu_n and pu_n will both be filled with -9999.99
process.diphotonAnalyzer.PUDataFileName = 'PileupDataAug10thHistogram.root' #DataPileUp
process.diphotonAnalyzer.PUMCFileName = cms.untracked.string(options.pumcfilename)  #"MC PileUP"
process.diphotonAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonAnalyzer.PUMCHistName = "MCPileUpHisto" #Name of histogram in PUMCFileName  Need to be binned to 80

# The full cmssw configuration path
process.path  = cms.Path(process.primaryVertexFilter*process.egmPhotonIDSequence*process.diphotonAnalyzer)
