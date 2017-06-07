from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'trigeffcalc_jetdata2016_wConeHE_crab3jobs'
config.section_('JobType')
config.JobType.psetName = 'exodiphotonanalyzer_minimal_cfg.py'
config.JobType.pyCfgParams = ['local=False','globalTag=80X_dataRun2_2016SeptRepro_v7','TrigEffOnly=True','addConeHE=True']
config.JobType.pluginName = 'Analysis'
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/bchiari1/noreplica/diphotonProject/scratch/photon_eff_calc'
config.Data.publication = False
config.Data.unitsPerJob = 200
config.Data.totalUnits = -1
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).

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

    config.General.requestName = 'JetHT_Run2016H2'
    config.Data.inputDataset = '/JetHT/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_Run2016H3'
    config.Data.inputDataset = '/JetHT/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_Run2016G'
    config.Data.inputDataset = '/JetHT/Run2016G-03Feb2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_Run2016F'
    config.Data.inputDataset = '/JetHT/Run2016F-18Apr2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_Run2016E'
    config.Data.inputDataset = '/JetHT/Run2016E-18Apr2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_Run2016D'
    config.Data.inputDataset = '/JetHT/Run2016D-18Apr2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_Run2016C'
    config.Data.inputDataset = '/JetHT/Run2016C-18Apr2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'JetHT_Run2016B'
    config.Data.inputDataset = '/JetHT/Run2016B-18Apr2017_ver2-v1/MINIAOD'
    submit(config)
