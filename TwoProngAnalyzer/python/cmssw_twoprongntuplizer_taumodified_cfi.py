import FWCore.ParameterSet.Config as cms

twoprongNtuplizer = cms.EDAnalyzer('TwoProngAnalyzer',
                                  # two-prong object options
                                  candidateMinPt = cms.untracked.double(20.0),
                                  candidateAbsMaxEta = cms.untracked.double(2.5),
                                  candidateTrackAsymmetryCut = cms.untracked.double(0.2),
                                  candidatePhotonAsymmetryCut = cms.untracked.double(0.15),
                                  candidateOptionalExtraTrack = cms.untracked.bool(True), # from False
                                  chargedHadronPairMinDR = cms.untracked.double(0.1), # from 0.05
                                  chargedHadronMinPt = cms.untracked.double(3.0),
                                  photonPtCut = cms.untracked.double(3.0),
                                  photonPhiBoxSize = cms.untracked.double(0.4), # from 0.2
                                  photonEtaBoxSize = cms.untracked.double(0.1), # from 0.05
                                  isolationConeR = cms.untracked.double(0.3),
                                  chargedIsoCut = cms.untracked.double(0.1),
                                  neutralIsoCut = cms.untracked.double(0.1),
                                  egammaIsoCut = cms.untracked.double(0.1),
                                  chargedIsoLooseMax = cms.untracked.double(0.3),
                                  neutralIsoLooseMax = cms.untracked.double(0.3),
                                  egammaIsoLooseMax = cms.untracked.double(0.3),
                                  generatorMatchDR = cms.untracked.double(0.1),
                                  # high-pt-photon-id options
                                  rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                  addPhotonCutDrConeHE = cms.untracked.bool(False),
                                  includeOldPhotons = cms.untracked.bool(True),
                                  # HLT paths
                                  bits = cms.InputTag("TriggerResults","","HLT"),
                                  prescales = cms.InputTag("patTrigger"),
                                  objects = cms.InputTag("selectedPatTrigger"),
                                  # MC options
                                  includeSignalGenParticles = cms.untracked.bool(True),
                                  includeMCInfo = cms.untracked.bool(False),
                                  mcXS = cms.untracked.double(1.0),
                                  mcN = cms.untracked.double(1.0),
                                  runningOnTauTauMC = cms.untracked.bool(False),
                                  # control output
                                  makeTrees = cms.untracked.bool(True),
                                  includeAllCandObjects = cms.untracked.bool(True),
                                  includeAllLooseObjects = cms.untracked.bool(False),
                                  debug = cms.untracked.bool(False),
                                  # histos
                                  fakeRateHistos = cms.untracked.bool(False),
                                  triggerEffHistos = cms.untracked.bool(False),
                                  twoprongYieldHistos = cms.untracked.bool(False),
                                  stackedDalitzHistos = cms.untracked.bool(False),
                                  # filters
                                  filterOnPhoton = cms.untracked.bool(False),
                                  filterOnTwoProng = cms.untracked.bool(False),
                                  )
