import CRABClient
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
from multiprocessing import Process
###############
import crab_multi_helper
testfile = "/cms/chiarito/samples/miniaod/mc/GJets_40to100_RunIISummer16_MINIAOD_numEvent10000.root"
###############
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_GJets_2016'
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['globalTag=mc2016', 'includeMCInfoBranches=True', 'oldData=True']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/bchiari1/cms_area/twoprong/trees/gjets2016/'
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

    config.General.requestName = 'GJets_HT-40To100'
    config.Data.inputDataset = '/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=20810','mcN=4858154'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-100To200'
    config.Data.inputDataset = '/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=9223','mcN=5050534'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-200To400'
    config.Data.inputDataset = '/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=2303','mcN=10350849'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-400To600'
    config.Data.inputDataset = '/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=274.5','mcN=2530341'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-600ToInf'
    config.Data.inputDataset = '/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=93.52','mcN=2616911'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-40To100_ext1'
    config.Data.inputDataset = '/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=20780','mcN=4467985'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-100To200_ext1'
    config.Data.inputDataset = '/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=9241','mcN=5131873'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-200To400_ext1'
    config.Data.inputDataset = '/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=2302','mcN=10036487'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-400To600_ext1'
    config.Data.inputDataset = '/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=275.4','mcN=2529729'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'GJets_HT-600ToInf_ext1'
    config.Data.inputDataset = '/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=93.36','mcN=2463946'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]
