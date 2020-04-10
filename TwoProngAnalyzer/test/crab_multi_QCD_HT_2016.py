from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from multiprocessing import Process
###############
import crab_multi_helper
import sys
testfile = "/cms/chiarito/samples/miniaod/mc/QCD_HT100to200_RunIISummer16_MINIAOD_numEvent10000.root"
###############
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_QCD_HT_2016'
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['globalTag=mc2016', 'includeMCInfoBranches=True', 'includeLooseTwoProngs=True']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/bchiari1/cms_area/twoprong/trees/no_filter/qcd2016_HT/'
config.Data.publication = False
config.Data.unitsPerJob = 250000
config.Data.totalUnits = -1
config.Data.splitting = 'EventAwareLumiBased'
config.Data.lumiMask = ''
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_Rutgers'
###############
crab_multi_helper.modify_config(config)
if crab_multi_helper.options.command:
  crab_multi_helper.print_command(config.JobType.psetName, config.JobType.pyCfgParams, testfile)
  sys.exit()
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

    config.General.requestName = 'QCD_HT_50to100'
    config.Data.inputDataset = '/QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=4097049'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_100to200'
    config.Data.inputDataset = '/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=80684349'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_200to300'
    config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=18722416'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_300to500'
    config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=17035891'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_500to700'
    config.Data.inputDataset = '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=18929951'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_700to1000'
    config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=15629253'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_1000to1500'
    config.Data.inputDataset = '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=4767100'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_1500to2000'
    config.Data.inputDataset = '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=3970819'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_2000toInf'
    config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=1991645'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_200to300_ext1'
    config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=38857944'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_300to500_ext1'
    config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=37502012'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_500to700_ext1'
    config.Data.inputDataset = '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=43341392'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_700to1000_ext1'
    config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=29783527'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_1000to1500_ext1'
    config.Data.inputDataset = '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=10360193'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_1500to2000_ext1'
    config.Data.inputDataset = '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=7855883'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]

    config.General.requestName = 'QCD_HT_2000toInf_ext1'
    config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    config.JobType.pyCfgParams.extend(['mcXS=1','mcN=4047360'])
    #submit(config)
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    config.JobType.pyCfgParams = config.JobType.pyCfgParams[:-2]
