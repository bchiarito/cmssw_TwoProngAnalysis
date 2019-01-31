import FWCore.ParameterSet.Config as cms

twoprongNtuplizer = cms.EDAnalyzer('TwoProngAnalyzer',
                                  # two-prong object options
                                  dontIncludeTwoProngs = cms.untracked.bool(False),
                                  twoprong_MinPt = cms.untracked.double(20.0),
                                  twoprong_AbsMaxEta = cms.untracked.double(2.5),
                                  twoprong_MinTrackAsymmetry = cms.untracked.double(0.2),
                                  twoprong_MinPhotonAsymmetry = cms.untracked.double(0.15),
                                  twoprong_optionalExtraTrack = cms.untracked.bool(False),
                                  twoprong_chargedHadronPairMinDR = cms.untracked.double(0.05),
                                  twoprong_chargedHadronMinPt = cms.untracked.double(3.0),
                                  twoprong_photonPtCut = cms.untracked.double(3.0),
                                  twoprong_photonPhiBoxSize = cms.untracked.double(0.2),
                                  twoprong_photonEtaBoxSize = cms.untracked.double(0.05),
                                  twoprong_isolationConeR = cms.untracked.double(0.3),
                                  twoprong_chargedIsoCut = cms.untracked.double(0.1),
                                  twoprong_neutralIsoCut = cms.untracked.double(0.1),
                                  twoprong_egammaIsoCut = cms.untracked.double(0.1),
                                  twoprong_chargedIsoLooseMax = cms.untracked.double(0.3),
                                  twoprong_neutralIsoLooseMax = cms.untracked.double(0.3),
                                  twoprong_egammaIsoLooseMax = cms.untracked.double(0.3),
                                  twoprong_generatorMatchDR = cms.untracked.double(0.1),
                                  twoprong_includeDalitzVariables = cms.untracked.bool(False),
                                  # high-pt-photon-id options
                                  rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                  includeOldPhotons = cms.untracked.bool(False),
                                  includeBasePhotons = cms.untracked.bool(False),
                                  # MC options
                                  includeSignalGenParticles = cms.untracked.bool(False),
                                  includeZDecayGenParticles = cms.untracked.bool(False),
                                  includeMCInfo = cms.untracked.bool(False),
                                  mcXS = cms.untracked.double(1.0),
                                  mcN = cms.untracked.double(1.0),
                                  # optional branches
                                  includeCandTwoProngs = cms.untracked.bool(True),
                                  includeLooseTwoProngs = cms.untracked.bool(False),
                                  includeZTauHadBranches = cms.untracked.bool(False),
                                  includeZMuMuBranches = cms.untracked.bool(False),
                                  usePatTauForZPreBranches = cms.untracked.bool(False),
                                  # filters
                                  filterOnPhoton = cms.untracked.bool(False),
                                  filterOnTwoProng = cms.untracked.bool(False),
                                  filterForABCDStudy = cms.untracked.bool(False),
                                  # control output
                                  debug = cms.untracked.bool(False),
                                  makeTrees = cms.untracked.bool(True),
                                  includeDalitzHistos = cms.untracked.bool(False),
                                  )
