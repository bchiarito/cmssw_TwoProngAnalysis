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
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "True includes all output, False removes most of the per event output")
options.register("sample",
                "local",
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "which sample we want to run over")
options.register("out",
                'Trial',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "output file name")
options.setDefault('maxEvents', 10)
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
elif sample == "jet":
    isMC = False
    isSignal = False
    readFiles.extend( [
       '/store/data/Run2016G/JetHT/MINIAOD/03Feb2017-v1/100000/006E7AF2-AEEC-E611-A88D-7845C4FC3B00.root' ] )
elif sample == "photon":
    isMC = False
    isSignal = False
    readFiles.extend( [
       '/store/data/Run2016G/SinglePhoton/MINIAOD/03Feb2017-v1/110000/00F4619E-9BEB-E611-8D7A-002590494BE2.root' ] )
else:
    print "Not a valid sample name!!"
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring( readFiles ))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

# JSON file to choose what lumis to process
if doLumis==True: 
    import FWCore.PythonUtilities.LumiList as LumiList
    goodlumis = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    process.source.lumisToProcess = LumiList.LumiList(filename = goodlumis).getVLuminosityBlockRange()

# Allow unscheduled
process.options.allowUnscheduled = cms.untracked.bool(True)

# Global Tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = options.globalTag

# Geometry for photon saturation 
process.load("Configuration.StandardSequences.GeometryDB_cff")

# Output file service
process.TFileService = cms.Service( "TFileService", fileName = cms.string( options.out + '_ExoDiPhotonAnalyzer.root') )

# filter on vertices
vtxCollName = 'offlineSlimmedPrimaryVertices'
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag(vtxCollName),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
)

### what is this for?
# Setup VID for EGM ID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

# the ntuplizer
process.diphotonAnalyzer = cms.EDAnalyzer('ExoDiPhotonAnalyzer',
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
                                  rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                  gedphotonsMiniAOD = cms.InputTag("slimmedPhotons"), # new block of code
                                  bits = cms.InputTag("TriggerResults","","HLT"),
                                  prescales = cms.InputTag("patTrigger"),
                                  objects = cms.InputTag("selectedPatTrigger"),
                                  )
process.diphotonAnalyzer.isMC = cms.untracked.bool(isMC)
process.diphotonAnalyzer.isSignal = cms.untracked.bool(isSignal)
process.diphotonAnalyzer.debug = cms.untracked.bool(options.debug)
process.diphotonAnalyzer.noTreeOnlyFakeRateHistos = cms.untracked.bool(False)
process.diphotonAnalyzer.omitChargedDecayCode = cms.bool(False)

# The full cmssw configuration path
process.path  = cms.Path(process.primaryVertexFilter * process.egmPhotonIDSequence * process.diphotonAnalyzer)
