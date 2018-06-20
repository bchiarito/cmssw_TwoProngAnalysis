from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'ztoll_truth3iprong_ex_con1_box01x04_dr15_taketwo_crab3jobs'
config.section_('JobType')
config.JobType.psetName = 'exodiphotonanalyzer_samples_cfg.py'
config.JobType.pyCfgParams = ['local=False','globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6','mcInfo=True','taupreselection=True','optionalExtraTrack=False','isTauTau=True','includeCands=True','constituentMinPt=1','photonBoxEta=0.1','photonBoxPhi=0.4','trackDR=0.15']
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2500 # default = 2000
config.JobType.maxJobRuntimeMin = 1315 # default = 1315 (22 hrs)
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/bchiari1/cms_area/twoprongstudies/trees/ztoll/truth3iprong_ex_con1_box01x04_dr15_taketwo'
config.Data.publication = False
config.Data.unitsPerJob = 500000
config.Data.totalUnits = -1
config.Data.splitting = 'EventAwareLumiBased'
config.Data.lumiMask = ''
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

    config.General.requestName = 'DYJetsToLL_M50_ext1'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=1'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'DYJetsToLL_M50_ext2'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=1'])
    submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'DYJetsToLL_M10to50'
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=1'])
    #submit(config)
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]
