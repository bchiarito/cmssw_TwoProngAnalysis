# Command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ("python")
# required
options.register("sample", "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "which sample we want to run over")
options.register("globalTag", "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "global tag to use when running")
# usually required
options.register("isSignal", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Specify singal MC for looking for Phi and omega gen particles")
options.register("isTauTau", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Specify Z->ll MC")
options.register("mcInfo", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "include mc weight in Ttree")
# other
options.register("mcXS", 1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "mc cross section, if desired to be filled in trees")
options.register("mcN", 1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "mc number generated, if desired to be filled in trees")
options.register("out", "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "output file name")
options.register("debug", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "True includes all output, False removes most of the per event output")
options.register("doLumis", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "use a JSON file to specify lumis")
options.register("originalGeometry", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "use the original loads for geometry, taken from original diphoton ntuplizer")
# two-prong object definition
options.register("standardTwoProng", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register("tauModifiedTwoProng", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register("commandLineTwoProng", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
# two-prong object definition detailed
options.register("optionalExtraTrack", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
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
options.register("addConeHE", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Add cut to high-pt-photon-id: Cone based HE < 0.05")
# output specification
options.register("ntuples", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Add ntuples (Ttrees) to output")
options.register("includeCands", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Include all cand twoprongs in ntuple")
options.register("includeLoose", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Include Loose twoprongs in ntuple")
options.register("fakeRateHistos", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register("triggerEffHistos", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register("twoprongYieldHistos", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.register("stackedDalitzHistos", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "")
options.setDefault("maxEvents", 10)
options.parseArguments()

# shortcut settings
if options.sample == 'eta125':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi125_Eta/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185319/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'eta125'
if options.sample == 'eta300':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi300_Eta/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185421/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'eta300'
if options.sample == 'eta500':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi500_Eta/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185444/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'eta500'
if options.sample == 'eta750':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi750_Eta/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185603/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'eta750'
if options.sample == 'eta1000':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi1000_Eta/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185641/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'eta1000'
if options.sample == 'etaprime125':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi125_Etaprime/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_190054/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'etaprime125'
if options.sample == 'etaprime300':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi300_Etaprime/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185936/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'etaprime300'
if options.sample == 'etaprime500':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi500_Etaprime/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185904/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'etaprime500'
if options.sample == 'etaprime750':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi750_Etaprime/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185807/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'etaprime750'
if options.sample == 'etaprime1000':
  options.sample = 'file:/cms/jrj90/eos/twoprong_generation/etaetaprime_run/july132018/Phi1000_Etaprime/76X_mcRun2_asymptotic_v12_Run2_25ns_MINIAOD/180713_185720/0000/MINIAOD_1.root'
  options.globalTag = 'mc2016'
  options.isSignal = True
  options.out = 'etaprime1000'

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
if options.globalTag == 'mc2016':
  process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
elif options.globalTag == 'mc2015':
  process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif options.globalTag == 'data2015':
  process.GlobalTag.globaltag = '74X_dataRun2_HLT_v3'
elif options.globalTag == 'data2016':
  process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
else:
  process.GlobalTag.globaltag = options.globalTag
print "Using GlobalTag: ", process.GlobalTag.globaltag

# Allow unscheduled
process.options.allowUnscheduled = cms.untracked.bool(True)

# Geometry for photon saturation 
if not options.originalGeometry:
  process.load("Configuration.Geometry.GeometryECALHCAL_cff") # default
else:
  process.load("Configuration.StandardSequences.GeometryDB_cff")

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

# JSON file to choose what lumis to process
if not options.doLumis=="":
    import FWCore.PythonUtilities.LumiList as LumiList
    goodlumis = options.doLumis
    process.source.lumisToProcess = LumiList.LumiList(filename = goodlumis).getVLuminosityBlockRange()

# the ntuplizer
if options.commandLineTwoProng:
  process.twoprongNtuplizer = cms.EDAnalyzer('ExoDiPhotonAnalyzer',
                                  # two-prong object options
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
                                  neutralIsoCut = cms.untracked.double(options.neutralIsoCut),
                                  egammaIsoCut = cms.untracked.double(options.egammaIsoCut),
                                  chargedIsoLooseMax = cms.untracked.double(0.3),
                                  neutralIsoLooseMax = cms.untracked.double(0.3),
                                  egammaIsoLooseMax = cms.untracked.double(0.3),
                                  generatorMatchDR = cms.untracked.double(0.1),
                                  # high-pt-photon-id options
                                  rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                  addPhotonCutDrConeHE = cms.untracked.bool(options.addConeHE),
                                  includeOldPhotons = cms.untracked.bool(True),
                                  # HLT paths
                                  bits = cms.InputTag("TriggerResults","","HLT"),
                                  prescales = cms.InputTag("patTrigger"),
                                  objects = cms.InputTag("selectedPatTrigger"),
                                  # MC options
                                  includeSignalGenParticles = cms.untracked.bool(options.isSignal),
                                  includeMCInfo = cms.untracked.bool(options.mcInfo),
                                  mcXS = cms.untracked.double(options.mcXS),
                                  mcN = cms.untracked.double(options.mcN),
                                  runningOnTauTauMC = cms.untracked.bool(options.isTauTau),
                                  # control output
                                  makeTrees = cms.untracked.bool(options.ntuples),
                                  includeAllCandObjects = cms.untracked.bool(options.includeCands),
                                  includeAllLooseObjects = cms.untracked.bool(options.includeLoose),
                                  debug = cms.untracked.bool(options.debug),
                                  # histos
                                  fakeRateHistos = cms.untracked.bool(options.fakeRateHistos),
                                  triggerEffHistos = cms.untracked.bool(options.triggerEffHistos),
                                  twoprongYieldHistos = cms.untracked.bool(options.twoprongYieldHistos),
                                  stackedDalitzHistos = cms.untracked.bool(options.stackedDalitzHistos),
                                  )
elif options.tauModifiedTwoProng:
  process.load('DiPhotonAnalysis.ExoDiPhotonAnalyzer.cmssw_twoprongntuplizer_taumodified_cfi')
elif options.standardTwoProng:
  process.load('DiPhotonAnalysis.ExoDiPhotonAnalyzer.cmssw_twoprongntuplizer_standard_cfi')
else:
  print "must select one twoprong version"
# settings always overwritten by command line
process.twoprongNtuplizer.includeMCInfo = options.mcInfo
process.twoprongNtuplizer.includeSignalGenParticles = options.isSignal
process.twoprongNtuplizer.runningOnTauTauMC = options.isTauTau
process.twoprongNtuplizer.mcXS = options.mcXS
process.twoprongNtuplizer.mcN = options.mcN
process.twoprongNtuplizer.debug = options.debug

# The path
process.path  = cms.Path(process.egmPhotonIDSequence * process.twoprongNtuplizer)
