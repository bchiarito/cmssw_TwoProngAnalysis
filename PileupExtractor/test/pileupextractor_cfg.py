import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms/store/user/charaf/4A7B228E-B3ED-E111-98D3-00266CFFC544_RSGravGGkMpl01M3250Summer12AODSIM.root'
    )
)


# need to introduce the global tag now
# because the L1GtUtils method needs to fetch records...
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#use the right global tag!
process.GlobalTag.globaltag = 'START53_V7G::All'

# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('PileUpMC.root')
)


#load diphoton PU extractor
process.load("DiPhotonAnalysis.PileupExtractor.pileupextractor_cfi")
##process.diphotonPileUpExtractor.pileupCollection = cms.untracked.InputTag("addPileupInfo")


process.p = cms.Path(process.diphotonPileUpExtractor)
