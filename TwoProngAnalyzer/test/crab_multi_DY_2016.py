from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from multiprocessing import Process
###############
import crab_multi_helper
testfile = "/cms/chiarito/samples/miniaod/mc/DYJetsToLL_M50_RunIISumer16_MINIAOD_numEvent10000.root"
###############
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_DY_2016'
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['globalTag=mc2016', 'includeMCInfoBranches=True', 'oldData=True']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/bchiari1/cms_area/twoprong/trees/dy2016/'
config.Data.publication = False
config.Data.unitsPerJob = 250000
config.Data.totalUnits = -1
config.Data.splitting = 'EventAwareLumiBased'
config.Data.lumiMask = ''
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_Rutgers'
###############
crab_multi_helper.modify_config(config, testfile)
###############

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
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcN=35291566','mcXS=18610'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'DY_50_ext1'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcN=49144274','mcXS=5765.4'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'DY_50_ext2'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcN=96658943','mcXS=5765.4'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]
