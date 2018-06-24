# Command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register("sample", "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "which sample we want to run over")
options.register('data2015', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "running on 2015 data needs different geometry include")
options.register("out", '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "output file name")
options.register('debug', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "True includes all output, False removes most of the per event output")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "global tag to use when running")
options.register('isSignal', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Specify singal MC for looking for Phi and omega gen particles")
options.register('isTauTau', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Specify Z->ll MC")
options.register('mcInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "include mc weight in Ttree")
options.register("mcXS", 1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "mc cross section, if desired to be filled in trees")
options.register("mcN", 1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "mc number generated, if desired to be filled in trees")
# two-prong object definition
options.register('optionalExtraTrack', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register("trackDR", 0.05, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("minPt", 20.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("maxEta", 2.5, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("constituentMinPt", 3.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("trackMinPt", 3.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("photonMinPt", 3.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("trackAsym", 0.2, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("photonAsym", 0.15, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("photonBoxPhi", 0.2, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("photonBoxEta", 0.05, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("chargedIsoCut", 0.1, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("neutralIsoCut", 0.1, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
options.register("egammaIsoCut", 0.1, VarParsing.multiplicity.singleton, VarParsing.varType.float, "")
# photon definition
options.register('addConeHE', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Add cut to high-pt-photon-id: Cone based HE < 0.05")
# output specification
options.register('ntuples', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Add ntuples (Ttrees) to output")
options.register('includeCands', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Include all cand twoprongs in ntuple")
options.register('includeLoose', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Include Loose twoprongs in ntuple")
options.register('fakeRateHistos', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register('triggerEffHistos', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register('twoprongYieldHistos', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register('stackedDalitzHistos', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.setDefault('maxEvents', 10)
options.parseArguments()

# Begin configuration
import FWCore.ParameterSet.Config as cms
process = cms.Process("TwoProngNtuplizer")

# Log messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.debug) )

# Source
readFiles = []
readFiles.extend( [ options.sample ] )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring( readFiles ))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

# Output files
if options.out == "":
  pre = 'TwoProngNtuplizer'
else:
  pre = 'TwoProngNtuplizer_'
post = '.root'
process.TFileService = cms.Service( "TFileService", fileName = cms.string( pre + options.out + post ) )

# Global Tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = options.globalTag

# Allow unscheduled
process.options.allowUnscheduled = cms.untracked.bool(True)

# Geometry for photon saturation 
if not options.data2015:
  process.load("Configuration.StandardSequences.GeometryDB_cff")
else:
  process.load("Configuration.Geometry.GeometryECALHCAL_cff")

# filter on vertices
vtxCollName = 'offlineSlimmedPrimaryVertices'
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag(vtxCollName),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
)

# adds computation of more Photon ID decisions, this block comes from high-pt-id code, but these ids are not currently being use in high-pt-id
# included for reference and for agreement with high-pt-id framework
# Setup VID for EGM ID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

# the ntuplizer
process.diphotonAnalyzer = cms.EDAnalyzer('ExoDiPhotonAnalyzer',
                                  # two-prong object
                                  candidateMinPt = cms.untracked.double(options.minPt),
                                  candidateAbsMaxEta = cms.untracked.double(options.maxEta),
                                  candidateTrackAsymmetryCut = cms.untracked.double(options.trackAsym),
                                  candidatePhotonAsymmetryCut = cms.untracked.double(options.photonAsym),
                                  candidateOptionalExtraTrack = cms.untracked.bool(options.optionalExtraTrack),
                                  chargedHadronPairMinDR = cms.untracked.double(options.trackDR),
                                  chargedHadronMinPt = cms.untracked.double(options.constituentMinPt),
                                  photonPtCut = cms.untracked.double(options.constituentMinPt),
                                  photonPhiBoxSize = cms.untracked.double(options.photonBoxPhi),
                                  photonEtaBoxSize = cms.untracked.double(options.photonBoxEta),
                                  isolationConeR = cms.untracked.double(0.3),
                                  chargedIsoCut = cms.untracked.double(options.chargedIsoCut),
                                  chargedIsoLooseMax = cms.untracked.double(0.3),
                                  neutralIsoCut = cms.untracked.double(options.neutralIsoCut),
                                  neutralIsoLooseMax = cms.untracked.double(0.3),
                                  egammaIsoCut = cms.untracked.double(options.egammaIsoCut),
                                  egammaIsoLooseMax = cms.untracked.double(0.3),
                                  generatorMatchDR = cms.untracked.double(0.1),
                                  # high-pt-photon-id options
                                  rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                  # HLT paths
                                  bits = cms.InputTag("TriggerResults","","HLT"),
                                  prescales = cms.InputTag("patTrigger"),
                                  objects = cms.InputTag("selectedPatTrigger"),
                                  )
# Ntuplizer Options
process.diphotonAnalyzer.addPhotonCutDrConeHE = cms.untracked.bool(options.addConeHE)
process.diphotonAnalyzer.includeAllCandObjects = cms.untracked.bool(options.includeCands)
process.diphotonAnalyzer.includeOldPhotons = cms.untracked.bool(False)
process.diphotonAnalyzer.debug = cms.untracked.bool(options.debug)
process.diphotonAnalyzer.includeAllLooseObjects = cms.untracked.bool(options.includeLoose)
process.diphotonAnalyzer.includeSignalGenParticles = cms.untracked.bool(options.isSignal)
process.diphotonAnalyzer.makeTrees = cms.untracked.bool(options.ntuples)
process.diphotonAnalyzer.fakeRateHistos = cms.untracked.bool(options.fakeRateHistos)
process.diphotonAnalyzer.triggerEffHistos = cms.untracked.bool(options.triggerEffHistos)
process.diphotonAnalyzer.twoprongYieldHistos = cms.untracked.bool(options.twoprongYieldHistos)
process.diphotonAnalyzer.stackedDalitzHistos = cms.untracked.bool(options.stackedDalitzHistos)
process.diphotonAnalyzer.includeMCInfo = cms.untracked.bool(options.mcInfo)
process.diphotonAnalyzer.mcXS = cms.untracked.double(options.mcXS)
process.diphotonAnalyzer.mcN = cms.untracked.double(options.mcN)
process.diphotonAnalyzer.runningOnTauTauMC = cms.untracked.bool(options.isTauTau)

# The full cmssw configuration path
process.path  = cms.Path(process.egmPhotonIDSequence * process.diphotonAnalyzer)
