from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('globalTag',
                '76X_dataRun2_16Dec2015_v0',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "global tag to use when running")
options.register('useAOD',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "whether or not to use AOD")
options.register('isMC',
                True,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "whether to run over data or MC")
options.register('isSignal',
                False,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.bool,
                "whether MC is a signal or background (used for storing gen informaton)")
options.register('pumcfilename',
                'PileUpMC_DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V7C-v1_rebinned.root',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "MC Pileup Filename")
options.register("sample",
                "eta",
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "which sample we want to run over"
)
options.setDefault('maxEvents', 10)
options.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

import sys
etaCMSFileList = cms.untracked.vstring(
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_1.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_3.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_4.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_7.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_8.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_9.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_11.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_12.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_14.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_15.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_16.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_17.root',
    'root://cmseos.fnal.gov//store/user/bchiari1/noreplica/diphotonProject/lhe/step3_output_19.root',
)
etaLocalFileList = cms.untracked.vstring(
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_1.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_11.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_12.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_14.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_15.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_16.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_17.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_19.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_3.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_4.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_7.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_8.root',
    'file:/afs/cern.ch/user/b/bchiari1/samples/lhe/diphoton750gev/eta/step3_output_9.root',
)
jetHTFileList = cms.untracked.vstring(
    '/store/data/Run2015C_25ns/JetHT/MINIAOD/16Dec2015-v1/20000/0A98D31C-49B5-E511-A886-0CC47A4C8EEA.root',
    '/store/data/Run2015C_25ns/JetHT/MINIAOD/16Dec2015-v1/20000/1079AE90-45B5-E511-9827-0002C94CDAE2.root',
)
singlePhotonFileList = cms.untracked.vstring(
    '/store/data/Run2015C_25ns/SinglePhoton/MINIAOD/16Dec2015-v1/50000/02E07C2B-90B2-E511-96BB-0CC47A78A3F4.root',
    '/store/data/Run2015C_25ns/SinglePhoton/MINIAOD/16Dec2015-v1/50000/12FEA727-90B2-E511-8AC8-0CC47A4D75F4.root',
)
qcdFileList = cms.untracked.vstring(
    '/store/mc/RunIIFall15MiniAODv2/QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/00AB0982-0AB9-E511-9CAD-441EA158FECC.root'
)

sample = options.sample
fList = None
if sample == "eta":
  fList = etaCMSFileList
elif sample == "etaLocal":
  fList = etaLocalFileList
elif sample == "jetHT":
  fList = jetHTFileList
elif sample == "singlePhoton":
  fList = singlePhotonFileList
elif sample == "qcd":
  fList = qcdFileList
process.source = cms.Source ("PoolSource", 
                fileNames      = fList,
                # debugVerbosity = cms.untracked.uint32(200),
                # debugFlag      = cms.untracked.bool(True)
 ) 

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
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

#load diphoton analyzer
process.load("DiPhotonAnalysis.ExoDiPhotonAnalyzer.exodiphotonanalyzer_withrho_cfi")
process.diphotonAnalyzer.rho25Correction = cms.InputTag("fixedGridRhoFastjetAll","rho") 
process.diphotonAnalyzer.ptMin = 50 # pt cut on all photons
process.diphotonAnalyzer.removeSpikes = False # ie spikes will be exlcuded from tree
process.diphotonAnalyzer.requireTightPhotons = False # ie only tight photons will be written 
process.diphotonAnalyzer.requireGenEventInfo = False #write MC info when running on MC
process.diphotonAnalyzer.isAOD = cms.bool(options.useAOD) # True=AOD, False=MiniAOD
process.diphotonAnalyzer.isMC = cms.untracked.bool(options.isMC)
process.diphotonAnalyzer.isSignal = cms.untracked.bool(options.isSignal)
process.diphotonAnalyzer.IDMethod = cms.untracked.string("highpt")
process.diphotonAnalyzer.PFIDCategory = cms.untracked.string("Loose")
process.diphotonAnalyzer.photonCollection = cms.untracked.InputTag("gedPhotons")
process.diphotonAnalyzer.jetCollection    = cms.InputTag("slimmedJets")
# If running on data the following four entries should not be changed. They are loaded into the analyzer as strings but in the case isMC = False then all the both old_pu_n and pu_n will both be filled with -9999.99
process.diphotonAnalyzer.PUDataFileName = 'PileupDataAug10thHistogram.root' #DataPileUp
process.diphotonAnalyzer.PUMCFileName = cms.untracked.string(options.pumcfilename)  #"MC PileUP"
process.diphotonAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonAnalyzer.PUMCHistName = "MCPileUpHisto" #Name of histogram in PUMCFileName  Need to be binned to 80

process.path  = cms.Path(process.primaryVertexFilter*process.egmPhotonIDSequence*process.diphotonAnalyzer)
