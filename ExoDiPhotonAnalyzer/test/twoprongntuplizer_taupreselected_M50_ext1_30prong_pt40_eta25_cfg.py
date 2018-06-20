# Command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
# required
options.register("sample", "signal125", VarParsing.multiplicity.singleton, VarParsing.varType.string, "which sample we want to run over")
options.register("out", '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "output file name")
# other
options.register('debug', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "True includes all output, False removes most of the per event output")
options.register('taupreselection', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "add Z->tau tau filter to path")
options.register('doLumis', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "let config file do lumi selection instead of CRAB - must be FALSE if using CRAB!")
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
# for crab
options.register('local', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Specify running on the command line, as opposed to crab/conder")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "global tag to use when running")
options.register('isSignal', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Specify singal MC for looking for Phi and omega gen particles")
options.register('isTauTau', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Specify Z->ll MC")
options.register('mcInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "include mc weight in Ttree")
options.register("mcXS", 1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "mc cross section, if desired to be filled in trees")
options.register("mcN", 1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "mc number generated, if desired to be filled in trees")

options.setDefault('maxEvents', 10)
options.parseArguments()

# set variables
isSignal = options.isSignal
isTauTau = options.isTauTau
doLumis = options.doLumis
sample = options.sample
globalTag = options.globalTag
outname = options.out
mcInfo = options.mcInfo

# Begin configuration
import FWCore.ParameterSet.Config as cms
process = cms.Process("ExoDiPhotonAnalysis")

# Log messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.debug) )

# Source
readFiles = []
readFiles.extend( [
'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_1.root'
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_10.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_11.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_12.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_13.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_14.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_15.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_16.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_17.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_18.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_19.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_2.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_20.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_21.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_22.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_23.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_24.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_25.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_26.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_27.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_28.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_29.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_3.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_30.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_31.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_32.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_33.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_34.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_35.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_36.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_37.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_38.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_39.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_4.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_40.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_41.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_42.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_43.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_44.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_45.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_46.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_47.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_48.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_49.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_5.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_50.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_51.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_52.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_53.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_54.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_55.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_56.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_57.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_58.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_59.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_6.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_60.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_61.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_62.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_63.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_64.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_65.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_66.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_67.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_68.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_69.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_7.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_70.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_71.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_72.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_73.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_74.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_75.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_76.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_77.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_78.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_79.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_8.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_80.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_81.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_82.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_83.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_84.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_85.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_86.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_87.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_88.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_89.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_9.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_90.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_91.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_92.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_93.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_94.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_95.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_96.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_97.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_98.root',
#'file:/cms/chiarito/eos/twoprongstudies/miniaod/ztoll/taufilt_30prong_pt40_eta25/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_ext1/180315_190729/0000/MINIAOD_99.root'
] )
isSignal = False
doLumis = False
mcInfo = True
globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
if outname == "":
  outname = "taupreselected"

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
process.diphotonAnalyzer.includeSignalGenParticles = cms.untracked.bool(isSignal)
process.diphotonAnalyzer.makeTrees = cms.untracked.bool(options.ntuples)
process.diphotonAnalyzer.fakeRateHistos = cms.untracked.bool(options.fakeRateHistos)
process.diphotonAnalyzer.triggerEffHistos = cms.untracked.bool(options.triggerEffHistos)
process.diphotonAnalyzer.twoprongYieldHistos = cms.untracked.bool(options.twoprongYieldHistos)
process.diphotonAnalyzer.stackedDalitzHistos = cms.untracked.bool(options.stackedDalitzHistos)
process.diphotonAnalyzer.includeMCInfo = cms.untracked.bool(mcInfo)
process.diphotonAnalyzer.mcXS = cms.untracked.double(options.mcXS)
process.diphotonAnalyzer.mcN = cms.untracked.double(options.mcN)
process.diphotonAnalyzer.runningOnTauTauMC = cms.untracked.bool(isTauTau)

# the tau filter
if options.taupreselection:
  process.tauFilter = cms.EDFilter('ZtoTauHadTruthSelector',
        filterByTruthDecayType = cms.untracked.vdouble(5.4,5.3),
        ptMin = cms.untracked.double(40.0),
        absEtaMax = cms.untracked.double(2.5),
        )

# The full cmssw configuration path
if not options.taupreselection:
  process.path  = cms.Path(process.primaryVertexFilter * process.egmPhotonIDSequence * process.diphotonAnalyzer)
else:
  process.path  = cms.Path(process.primaryVertexFilter * process.tauFilter * process.egmPhotonIDSequence * process.diphotonAnalyzer)
