from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_twoprongntuplizer_DY_reg2prong'
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['isTauTau=True', 'globalTag=mc2016']
config.JobType.pluginName = 'Analysis'
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/bchiari1/cms_area/twoprong/ztagandprobe/nofiltering/trees/reg2prong/DY/'
config.Data.outputDatasetTag = ""
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

    config.General.requestName = 'DY_10to50'
    config.Data.outputDatasetTag = "DY_10to50"
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcN=35291566','mcXS=16260'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'DY_10to50_signal'
    config.Data.outputDatasetTag = "DY_10to50_signal_MuonHadronicFilter"
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['DYsignal=True'])
    config.JobType.pyCfgParams.extend(['mcN=35291566','mcXS=16260'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-3]

    config.General.requestName = 'DY_10to50_bkg'
    config.Data.outputDatasetTag = "DY_10to50_bkg_MuonHadronicExcludedFilter"
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['DYbkg=True'])
    config.JobType.pyCfgParams.extend(['mcN=35291566','mcXS=16260'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-3]

    config.General.requestName = 'DY_50_ext1'
    config.Data.outputDatasetTag = "DY_50_ext1"
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcN=49144274','mcXS=4956'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'DY_50_ext1_signal'
    config.Data.outputDatasetTag = "DY_50_ext1_signal_MuonHadronicFilter"
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['DYsignal=True'])
    config.JobType.pyCfgParams.extend(['mcN=49144274','mcXS=4956'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-3]

    config.General.requestName = 'DY_50_ext1_bkg'
    config.Data.outputDatasetTag = "DY_50_ext1_bkg_MuonHadronicExcludedFilter"
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['DYbkg=True'])
    config.JobType.pyCfgParams.extend(['mcN=49144274','mcXS=4956'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-3]

    config.General.requestName = 'DY_50_ext2'
    config.Data.outputDatasetTag = "DY_50_ext2"
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcN=96658943','mcXS=4956'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'DY_50_ext2_signal'
    config.Data.outputDatasetTag = "DY_50_ext2_signal_MuonHadronicFilter"
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['DYsignal=True'])
    config.JobType.pyCfgParams.extend(['mcN=96658943','mcXS=4956'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-3]

    config.General.requestName = 'DY_50_ext2_bkg'
    config.Data.outputDatasetTag = "DY_50_ext2_bkg_MuonHadronicExcludedFilter"
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['DYbkg=True'])
    config.JobType.pyCfgParams.extend(['mcN=96658943','mcXS=4956'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-3]
