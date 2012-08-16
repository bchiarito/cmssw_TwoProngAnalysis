import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M1000_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M1250_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M1500_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M1750_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M2000_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M2250_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M2500_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M3000_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    'file:/afs/cern.ch/work/c/chenders/private/storage/diphotons/signal_mc/RSGravGG_kMpl001_M750_TuneZ2star_8TeV_pythia6_cff_py_GEN.root',
    
    
    )
)


# file for all histograms for all modules
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('diphoton_tree.root')
)


process.demo = cms.EDAnalyzer('ExoDiPhotonGENOnlyAnalyzer'
)


process.p = cms.Path(process.demo)
