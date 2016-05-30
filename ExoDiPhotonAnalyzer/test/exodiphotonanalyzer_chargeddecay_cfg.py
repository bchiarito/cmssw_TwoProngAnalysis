from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('globalTag',
                '74X_mcRun2_asymptotic_realisticBS_v1',
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
options.register('pumcfilename',
                'PileUpMC_DiPhotonJetsBox_M60_8TeV-sherpa_Summer12_DR53X-PU_S10_START53_V7C-v1_rebinned.root',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "MC Pileup Filename")
options.register("decayType",
                "eta",
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "which sample we want to run over"
)
options.setDefault('maxEvents', -1)
options.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process("ExoDiPhotonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

import sys
etaFileList = cms.untracked.vstring(
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_3.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_1.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_7.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_4.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_8.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_9.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_12.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_11.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_16.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_17.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_14.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_15.root',
    'file:/cms/data26/feigelis/SeesawProject/PostHadronizer/output/Scott/eta/step3_output_19.root',
)
sample = options.decayType
fList = None
if sample == "eta":
  fList = etaFileList
print "#########################################"
print "#"
print "#  Running over the %s signal sample"%sample
print "#"
print "#########################################"
process.source = cms.Source ("PoolSource", 
                fileNames      = fList,
                # debugVerbosity = cms.untracked.uint32(200),
                # debugFlag      = cms.untracked.bool(True)
 ) 

#use the right global tag!
##process.GlobalTag.globaltag = 'GR_P_V54::All'
##process.GlobalTag.globaltag = 'GR_E_V48::All'
##process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
# process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'
# process.GlobalTag.globaltag = options.globalTag

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
process.diphotonAnalyzer.isMC = cms.untracked.bool(options.isMC) # False by default, run with isMC=True for MC
process.diphotonAnalyzer.IDMethod = cms.untracked.string("highpt")
process.diphotonAnalyzer.PFIDCategory = cms.untracked.string("Loose")
process.diphotonAnalyzer.photonCollection = cms.untracked.InputTag("gedPhotons")
process.diphotonAnalyzer.jetCollection    = cms.string("slimmedJets")
# If running on data the following four entries should not be changed. They are loaded into the analyzer as strings but in the case isMC = False then all the both old_pu_n and pu_n will both be filled with -9999.99
process.diphotonAnalyzer.PUDataFileName = 'PileupDataAug10thHistogram.root' #DataPileUp
process.diphotonAnalyzer.PUMCFileName = cms.untracked.string(options.pumcfilename)  #"MC PileUP"
process.diphotonAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonAnalyzer.PUMCHistName = "MCPileUpHisto" #Name of histogram in PUMCFileName  Need to be binned to 80

process.path  = cms.Path(process.primaryVertexFilter*process.egmPhotonIDSequence*process.diphotonAnalyzer)
