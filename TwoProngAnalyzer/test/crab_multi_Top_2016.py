from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_twoprongntuplizer_Top'
config.section_('JobType')
config.JobType.psetName = '../cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['globalTag=mc2016', 'includeMCInfoBranches=True', 'includeBaseTwoProngs=False' ,'includeLooseTwoProngs=True']
config.JobType.pluginName = 'Analysis'
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/bchiari1/cms_area/twoprong/trees/Top/date'
config.Data.outputDatasetTag = "twoprongntuplizer"
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 200000
config.Data.totalUnits =  -1
config.Data.publication = False
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_Rutgers'


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

    config.General.requestName = 'TT'
    config.Data.outputDatasetTag = "twoprongntuplizer"
    config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcN=10139950','mcXS=509.4'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

#    config.General.requestName = 'ST_t_top'
#    config.Data.outputDatasetTag = "twoprongntuplizer"
#    config.Data.inputDataset = '/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    config.JobType.pyCfgParams.extend(['mcN=67240808','mcXS=123.3'])
#    submit(config)
#    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

#    config.General.requestName = 'ST_t_antitop'
#    config.Data.outputDatasetTag = "twoprongntuplizer"
#    config.Data.inputDataset = '/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#    config.JobType.pyCfgParams.extend(['mcN=38811017','mcXS=74.41'])
#    submit(config)
#    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

#    config.General.requestName = 'ST_tW_top'
#    config.Data.outputDatasetTag = "twoprongntuplizer"
#    config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#    config.JobType.pyCfgParams.extend(['mcN=6952830','mcXS=38.09'])
#    submit(config)
#    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

#    config.General.requestName = 'ST_tW_antitop'
#    config.Data.outputDatasetTag = "twoprongntuplizer"
#    config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#    config.JobType.pyCfgParams.extend(['mcN=6952830','mcXS=38.09'])
#    submit(config)
#    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]
