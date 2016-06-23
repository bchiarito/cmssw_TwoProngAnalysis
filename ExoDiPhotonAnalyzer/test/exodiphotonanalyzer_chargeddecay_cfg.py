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
genTestpythia6List = cms.untracked.vstring(
'/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_10.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_101.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_102.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_103.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_106.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_107.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_108.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_11.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_114.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_117.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_118.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_119.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_12.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_120.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_121.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_122.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_124.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_125.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_126.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_127.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_128.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_129.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_13.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_130.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_131.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_132.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_133.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_135.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_136.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_140.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_141.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_143.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_144.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_146.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_147.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_148.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_149.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_15.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_150.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_16.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_17.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_18.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_19.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_20.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_22.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_23.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_24.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_25.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_26.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_27.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_28.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_30.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_31.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_32.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_33.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_34.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_35.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_36.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_37.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_38.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_39.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_40.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_42.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_43.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_44.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_45.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_46.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_47.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_48.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_50.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_52.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_53.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_55.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_56.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_57.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_58.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_59.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_6.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_60.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_61.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_63.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_64.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_65.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_66.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_67.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_68.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_69.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_7.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_70.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_71.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_72.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_74.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_75.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_77.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_78.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_8.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_80.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_83.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_84.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_87.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia8_MiniAODv2/160621_132817/0000/MiniAODv2_9.root'
)

genTestpythia8List = cms.untracked.vstring(
    '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_1.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_10.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_11.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_12.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_13.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_14.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_15.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_16.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_17.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_18.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_19.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_2.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_20.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_21.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_22.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_23.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_24.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_25.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_26.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_27.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_28.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_29.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_3.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_30.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_31.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_32.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_33.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_34.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_35.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_36.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_37.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_38.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_39.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_4.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_40.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_41.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_42.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_43.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_44.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_45.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_46.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_47.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_48.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_49.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_5.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_50.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_51.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_52.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_53.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_54.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_55.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_56.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_57.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_58.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_59.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_6.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_60.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_61.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_62.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_63.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_64.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_65.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_66.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_67.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_68.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_69.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_7.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_70.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_71.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_72.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_73.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_74.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_75.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_76.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_77.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_78.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_79.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_8.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_80.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_81.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_82.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_83.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_84.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_85.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_86.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_87.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_88.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_89.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_9.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_90.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_91.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_92.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_93.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_94.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_95.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_96.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_97.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_98.root',
       '/store/user/bchiari1/noreplica/diphotonProject/miniaod/DiPhoton750GeVSignalEta/pythia6_Tauola_13TeV_forcedDecay_3pi0_MINIAOD/160620_201752/0000/MiniAODv2_99.root'
)

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
elif sample == "gen8":
  fList = genTestpythia8List
elif sample == "gen6":
  fList = genTestpythia6List
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

process.diphotonAnalyzer.debug = cms.untracked.bool(False)
# If running on data the following four entries should not be changed. They are loaded into the analyzer as strings but in the case isMC = False then all the both old_pu_n and pu_n will both be filled with -9999.99
process.diphotonAnalyzer.PUDataFileName = 'PileupDataAug10thHistogram.root' #DataPileUp
process.diphotonAnalyzer.PUMCFileName = cms.untracked.string(options.pumcfilename)  #"MC PileUP"
process.diphotonAnalyzer.PUDataHistName = "pileup" #Name of histogram in PUDataFileName Need to be binned to 80
process.diphotonAnalyzer.PUMCHistName = "MCPileUpHisto" #Name of histogram in PUMCFileName  Need to be binned to 80

process.path  = cms.Path(process.primaryVertexFilter*process.egmPhotonIDSequence*process.diphotonAnalyzer)
