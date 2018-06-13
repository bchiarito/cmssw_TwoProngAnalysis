from FWCore.ParameterSet.VarParsing import VarParsing
# Command line options
options = VarParsing ('python')
options.setDefault('maxEvents', 100)
options.register("sample", "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "which sample we want to run over")
options.parseArguments()

# Begin configuration
import FWCore.ParameterSet.Config as cms
process = cms.Process("ExoDiPhotonAnalysis")

# Log messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# Source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring( options.sample ))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

# Output files
outname = "TwoProngNtuplizer.root"
process.TFileService = cms.Service( "TFileService", fileName = cms.string( outname ) )

# Global Tag
globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = globalTag

# Allow unscheduled
process.options.allowUnscheduled = cms.untracked.bool(True)

# Geometry for photon saturation 
process.load("Configuration.StandardSequences.GeometryDB_cff")

# Filter on vertices
#vtxCollName = 'offlineSlimmedPrimaryVertices'
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                           vertexCollection = cms.InputTag(vtxCollName),
#                                           minimumNDOF = cms.uint32(4),
#                                           maxAbsZ = cms.double(24),	
#                                           maxd0 = cms.double(2)	
#)

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
                                  # two-prong object options
                                  candidateMinPt = cms.untracked.double(20.0),
                                  candidateAbsMaxEta = cms.untracked.double(2.5),
                                  candidateTrackAsymmetryCut = cms.untracked.double(0.2),
                                  candidatePhotonAsymmetryCut = cms.untracked.double(0.15),
                                  candidateOptionalExtraTrack = cms.untracked.bool(False),
                                  chargedHadronPairMinDR = cms.untracked.double(0.05),
                                  chargedHadronMinPt = cms.untracked.double(3.0),
                                  photonPtCut = cms.untracked.double(3.0),
                                  photonPhiBoxSize = cms.untracked.double(0.2),
                                  photonEtaBoxSize = cms.untracked.double(0.05),
                                  isolationConeR = cms.untracked.double(0.3),
                                  chargedIsoCut = cms.untracked.double(0.1),
                                  neutralIsoCut = cms.untracked.double(0.1),
                                  egammaIsoCut = cms.untracked.double(0.1),
                                  chargedIsoLooseMax = cms.untracked.double(0.3),
                                  neutralIsoLooseMax = cms.untracked.double(0.3),
                                  egammaIsoLooseMax = cms.untracked.double(0.3),
                                  generatorMatchDR = cms.untracked.double(0.1),
                                  # high-pt-photon-id options
                                  rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                  addPhotonCutDrConeHE = cms.untracked.bool(False),
                                  includeOldPhotons = cms.untracked.bool(False),
                                  # HLT paths
                                  bits = cms.InputTag("TriggerResults","","HLT"),
                                  prescales = cms.InputTag("patTrigger"),
                                  objects = cms.InputTag("selectedPatTrigger"),
                                  # MC options
                                  includeSignalGenParticles = cms.untracked.bool(False),
                                  includeMCInfo = cms.untracked.bool(False),
                                  mcXS = cms.untracked.double(1.0),
                                  mcN = cms.untracked.double(1.0),
                                  runningOnTauTauMC = cms.untracked.bool(False),
                                  # control output
                                  makeTrees = cms.untracked.bool(True),
                                  includeAllCandObjects = cms.untracked.bool(True),
                                  includeAllLooseObjects = cms.untracked.bool(False),
                                  debug = cms.untracked.bool(False),
                                  # histos
                                  fakeRateHistos = cms.untracked.bool(False),
                                  triggerEffHistos = cms.untracked.bool(False),
                                  twoprongYieldHistos = cms.untracked.bool(False),
                                  stackedDalitzHistos = cms.untracked.bool(False),
                                  )

process.path  = cms.Path(process.primaryVertexFilter * process.egmPhotonIDSequence * process.diphotonAnalyzer)
