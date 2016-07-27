from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'crab3test'
config.section_('JobType')
config.JobType.psetName = 'exodiphotonanalyzer_Data2015Commissioning_cfg.py'
config.JobType.pyCfgParams = ['globalTag=76X_dataRun2_16Dec2015_v0','useAOD=False','isMC=False']
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['PileupDataAug10thHistogram.root', 'PileUpMC.root', 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt']
config.JobType.outputFiles = ['ExoDiPhotonAnalyzer.root']
config.section_('Data')
config.Data.inputDataset = '/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.outputDatasetTag = ''
config.Data.outLFNDirBase = '/store/user/skaplan/noreplica/ExoDiPhotonNtuples/Data2015/'
config.Data.publication = False
config.Data.unitsPerJob = 100
config.Data.totalUnits = -1
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab3data2015jobs'

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

    config.General.requestName = 'DoubleEG_RunC'
    config.Data.inputDataset = '/DoubleEG/Run2015C_25ns-16Dec2015-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleEG_RunD'
    config.Data.inputDataset = '/DoubleEG/Run2015D-16Dec2015-v2/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_RunC'
    config.Data.inputDataset = '/JetHT/Run2015C_25ns-16Dec2015-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_RunD'
    config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_RunC'
    config.Data.inputDataset = '/SinglePhoton/Run2015C_25ns-16Dec2015-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_RunD'
    config.Data.inputDataset = '/SinglePhoton/Run2015D-16Dec2015-v1/MINIAOD'
    submit(config)
