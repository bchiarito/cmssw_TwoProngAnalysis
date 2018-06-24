from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_SinglePhoton_2015_50ns'
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_crab_cfg.py'
config.JobType.pyCfgParams = ['globalTag=GR_H_V58C', 'data2015=True']
config.JobType.pluginName = 'Analysis'
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/%s/cms_area/twoprong/trees/photon2015_50ns/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.unitsPerJob = 100
config.Data.totalUnits = 1
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_13TeV_16Dec2015ReReco_Collisions15_50ns_JSON_v2.txt'
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

    config.General.requestName = 'SinglePhoton_Run2015B'
    config.Data.inputDataset = '/SinglePhoton/Run2015B-16Dec2015-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_Run2015C'
    config.Data.inputDataset = '/SinglePhoton/Run2015C_50ns-16Dec2015-v1/MINIAOD'
    submit(config)
