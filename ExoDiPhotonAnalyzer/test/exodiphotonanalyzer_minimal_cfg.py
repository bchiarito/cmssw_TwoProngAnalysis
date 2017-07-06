# Command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('globalTag',
                '',
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
                '',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "output file name")
options.register('isSignal',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Specify singal MC for looking for Phi and omega gen particles")
options.register('local',
                True,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Specify running on the command line, as opposed to crab/conder")
options.register('addConeHE',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Add cut to high-pt-photon-id: Cone based HE < 0.05")
options.register('includeLoose',
                True,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Include Loose twoprongs in ntuple")
options.register('ntuples',
                True,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Add ntuples (Ttrees) to output")
options.register('fakeRateHistos',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Add ntuples (Ttrees) to output")
options.register('triggerEffHistos',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Add ntuples (Ttrees) to output")
options.register('twoprongYieldHistos',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "Add ntuples (Ttrees) to output")
options.register("minPt",
                0,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.float,
                "")
options.register("maxEta",
                2.5,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.float,
                "")
options.register("constituentMinPt",
                10,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.float,
                "")
options.setDefault('maxEvents', 100)
options.parseArguments()

# set variables
isSignal = options.isSignal
doLumis = options.doLumis
sample = options.sample
globalTag = options.globalTag
outname = options.out

# Begin configuration
import FWCore.ParameterSet.Config as cms
process = cms.Process("ExoDiPhotonAnalysis")

# Log messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.debug) )

# Source
readFiles = []
if sample == "signal":
    readFiles.extend( [
        'file:/cms/chiarito/samples/signal/MiniAODv2_Eta_generic.root' ] )
    if options.local:
      isSignal = True
      doLumis = False
      globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
      if outname == "":
        outname = "signal"

elif sample == "signal125":
    readFiles.extend( [
        'file:/cms/chiarito/samples/signal/MiniAODv2_Eta_125_pipipi0.root' ] )
    if options.local:
      isSignal = True
      doLumis = False
      globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
      if outname == "":
        outname = "signal125"
elif sample == "signal300":
    readFiles.extend( [
        'file:/cms/chiarito/samples/signal/MiniAODv2_Eta_300_pipipi0.root' ] )
    if options.local:
      isSignal = True
      doLumis = False
      globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
      if outname == "":
        outname = "signal300"
elif sample == "signal500":
    readFiles.extend( [
        'file:/cms/chiarito/samples/signal/MiniAODv2_Eta_500_pipipi0.root' ] )
    if options.local:
      isSignal = True
      doLumis = False
      globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
      if outname == "":
        outname = "signal500"
elif sample == "signal750":
    readFiles.extend( [
        'file:/cms/chiarito/samples/signal/MiniAODv2_Eta_750_pipipi0.root' ] )
    if options.local:
      isSignal = True
      doLumis = False
      globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
      if outname == "":
        outname = "signal750"
elif sample == "signal1000":
    readFiles.extend( [
        'file:/cms/chiarito/samples/signal/MiniAODv2_Eta_1000_pipipi0.root' ] )
    if options.local:
      isSignal = True
      doLumis = False
      globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
      if outname == "":
        outname = "signal1000"

elif sample == "jet":
    # 100k events
    readFiles.extend( [
       '/store/data/Run2016G/JetHT/MINIAOD/03Feb2017-v1/100000/006E7AF2-AEEC-E611-A88D-7845C4FC3B00.root' ] )
    if options.local:
      isSignal = False
      doLumis = True
      globalTag = "80X_dataRun2_2016SeptRepro_v7"
      if outname == "":
        outname = "jet"
elif sample == "photon":
    # 100k events
    readFiles.extend( [
       '/store/data/Run2016G/SinglePhoton/MINIAOD/03Feb2017-v1/110000/00F4619E-9BEB-E611-8D7A-002590494BE2.root' ] )
    if options.local:
      isSignal = False
      doLumis = True
      globalTag = "80X_dataRun2_2016SeptRepro_v7"
      if outname == "":
        outname = "photon"
elif options.local:
    print "Not a valid sample name!!"
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring( readFiles ))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

# Output files
if outname == "":
  pre = 'TwoProngNtuplizer'
else:
  pre = 'TwoProngNtuplizer_'
post = '.root'
process.TFileService = cms.Service( "TFileService", fileName = cms.string( pre + outname + post ) )

# Global Tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = globalTag

# JSON file to choose what lumis to process
if options.doLumis==True: 
    import FWCore.PythonUtilities.LumiList as LumiList
    goodlumis = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    process.source.lumisToProcess = LumiList.LumiList(filename = goodlumis).getVLuminosityBlockRange()

# Allow unscheduled
process.options.allowUnscheduled = cms.untracked.bool(True)

# Geometry for photon saturation 
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

# the ntuplizer
process.diphotonAnalyzer = cms.EDAnalyzer('ExoDiPhotonAnalyzer',
                                  # two-prong object
                                  candidateMinPt = cms.untracked.double(options.minPt),
                                  candidateAbsMaxEta = cms.untracked.double(options.maxEta),
                                  chargedHadronPairMinDR = cms.untracked.double(0.05),
                                  chargedHadronMinPt = cms.untracked.double(options.constituentMinPt),
                                  photonPtCut = cms.untracked.double(options.constituentMinPt),
                                  photonPhiBoxSize = cms.untracked.double(0.8),
                                  photonEtaBoxSize = cms.untracked.double(0.087),
                                  isolationConeR = cms.untracked.double(0.3),
                                  chargedIsoCut = cms.untracked.double(0.1),
                                  chargedIsoLooseMax = cms.untracked.double(0.3),
                                  neutralIsoCut = cms.untracked.double(0.1),
                                  neutralIsoLooseMax = cms.untracked.double(0.3),
                                  egammaIsoCut = cms.untracked.double(0.1),
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
process.diphotonAnalyzer.includeAllCandObjects = cms.untracked.bool(False)
process.diphotonAnalyzer.includeOldPhotons = cms.untracked.bool(False)
process.diphotonAnalyzer.debug = cms.untracked.bool(options.debug)
process.diphotonAnalyzer.includeAllLooseObjects = cms.untracked.bool(options.includeLoose)
process.diphotonAnalyzer.includeSignalGenParticles = cms.untracked.bool(isSignal)
process.diphotonAnalyzer.makeTrees = cms.untracked.bool(options.ntuples)
process.diphotonAnalyzer.fakeRateHistos = cms.untracked.bool(False)
process.diphotonAnalyzer.triggerEffHistos = cms.untracked.bool(False)
process.diphotonAnalyzer.twoprongYieldHistos = cms.untracked.bool(True)

# The full cmssw configuration path
process.path  = cms.Path(process.primaryVertexFilter * process.egmPhotonIDSequence * process.diphotonAnalyzer)
