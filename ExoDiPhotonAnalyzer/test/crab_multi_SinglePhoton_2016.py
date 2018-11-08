from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab_multi_SinglePhoton_2016'
config.section_('JobType')
config.JobType.psetName = 'cmssw_twoprongntuplizer_cfg.py'
config.JobType.pyCfgParams = ['globalTag=data2016', 'addConeHE=True', 'includeCands=False', 'filterOnPhoton=True']
config.JobType.pluginName = 'Analysis'
config.section_('Data')
config.Data.outLFNDirBase = '/store/user/%s/cms_area/twoprong/prelim/Nov5/photon2016/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.unitsPerJob = 100
config.Data.totalUnits = -1
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
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

    config.General.requestName = 'SinglePhoton_2016_RunH2'
    config.Data.inputDataset = '/SinglePhoton/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_2016_RunH3'
    config.Data.inputDataset = '/SinglePhoton/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_2016_RunG'
    config.Data.inputDataset = '/SinglePhoton/Run2016G-03Feb2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_2016_RunF'
    config.Data.inputDataset = '/SinglePhoton/Run2016F-18Apr2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_2016_RunE'
    config.Data.inputDataset = '/SinglePhoton/Run2016E-03Feb2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_2016_RunD'
    config.Data.inputDataset = '/SinglePhoton/Run2016D-03Feb2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_2016_RunC'
    config.Data.inputDataset = '/SinglePhoton/Run2016C-03Feb2017-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SinglePhoton_2016_RunB'
    config.Data.inputDataset = '/SinglePhoton/Run2016B-03Feb2017_ver2-v2/MINIAOD'
    submit(config)
