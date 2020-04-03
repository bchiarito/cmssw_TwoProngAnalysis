from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from multiprocessing import Process
###############
import crab_multi_helper
import sys
testfile = "/cms/chiarito/samples/miniaod/data/JetHT_2016_RunG_03Feb2017_MINIAOD_numEvent10000.root"
###############
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_JetHT_2016'
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['globalTag=data2016', 'includeLooseTwoProngs=True']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/%s/cms_area/twoprong/trees/no_filter/jet2016/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.unitsPerJob = 100
config.Data.totalUnits = -1
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
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
