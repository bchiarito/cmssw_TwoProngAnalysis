from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'crab3test'
config.section_('JobType')
config.JobType.psetName = 'exodiphotonanalyzer_Data2015Commissioning_cfg.py'
config.JobType.pyCfgParams = ['globalTag=76X_mcRun2_asymptotic_RunIIFall15DR76_v1','useAOD=False',"isMC=True"]
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['PileupDataAug10thHistogram.root', 'PileUpMC.root', 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt']
# config.JobType.outputFiles = ['ExoDiPhotonAnalyzer.root']
config.section_('Data')
config.Data.inputDataset = '/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.outputDatasetTag = ''
config.Data.outLFNDirBase = '/store/user/skaplan/noreplica/750GeVResonanceNtuples/BrandonSteveMerged/'
config.Data.publication = False
config.Data.unitsPerJob = 100000
config.Data.totalUnits = -1
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.lumiMask = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab3mcjobs'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = "GGJets_M-60To200_Pt-50"
    config.Data.inputDataset = "/GGJets_M-60To200_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "GGJets_M-200To500_Pt-50"
    config.Data.inputDataset = "/GGJets_M-200To500_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "GGJets_M-500To1000_Pt-50"
    config.Data.inputDataset = "/GGJets_M-500To1000_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "GGJets_M-1000To2000_Pt-50"
    config.Data.inputDataset = "/GGJets_M-1000To2000_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "GGJets_M-2000To4000_Pt-50"
    config.Data.inputDataset = "/GGJets_M-2000To4000_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "GGJets_M-4000To6000_Pt-50"
    config.Data.inputDataset = "/GGJets_M-4000To6000_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "GGJets_M-6000To8000_Pt-50"
    config.Data.inputDataset = "/GGJets_M-6000To8000_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "GGJets_M-8000To13000_Pt-50"
    config.Data.inputDataset = "/GGJets_M-8000To13000_Pt-50_13TeV-sherpa/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_5to10"
    config.Data.inputDataset = "/QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_10to15"
    config.Data.inputDataset = "/QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_15to30"
    config.Data.inputDataset = "/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_30to50"
    config.Data.inputDataset = "/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_50to80"
    config.Data.inputDataset = "/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_80to120"
    config.Data.inputDataset = "/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_120to170"
    config.Data.inputDataset = "/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"

    config.General.requestName = "QCD_Pt_170to300"
    config.Data.inputDataset = "/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_300to470"
    config.Data.inputDataset = "/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_470to600"
    config.Data.inputDataset = "/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_600to800"
    config.Data.inputDataset = "/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_800to1000"
    config.Data.inputDataset = "/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_1000to1400"
    config.Data.inputDataset = "/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_1400to1800"
    config.Data.inputDataset = "/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_1800to2400"
    config.Data.inputDataset = "/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_2400to3200"
    config.Data.inputDataset = "/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "QCD_Pt_3200toInf"
    config.Data.inputDataset = "/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
    submit(config)
