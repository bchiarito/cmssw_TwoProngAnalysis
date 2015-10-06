import FWCore.ParameterSet.Config as cms

diphotonSignalMCAnalyzer = cms.EDAnalyzer('ExoDiPhotonSignalMCAnalyzer',
                                          photonCollection = cms.untracked.InputTag("gedPhotons"),
                                          ptMin = cms.untracked.double(10),
                                          # careful with the HLT process name for MC samples!
                                          # it changes every time there is a re-reco done on the same RAW files!
                                          # eg for Spring10 35X samples, I believe the name is "REDIGI"
                                          hltResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                          rho25Correction = cms.InputTag("kt6PFJets25","rho"),
                                          pileupCorrection = cms.untracked.InputTag("addPileupInfo"),
                                          removeSpikes = cms.untracked.bool(False),
                                          requireTightPhotons = cms.untracked.bool(False),
                                          requireGenEventInfo = cms.untracked.bool(False),
                                          PUMCFileName = cms.untracked.string("PileUpMC.root"),
                                          PUDataFileName = cms.untracked.string("PileupDataAug10thHistogram.root"),
                                          PUMCHistName = cms.untracked.string("MCPileUpHistogram"),
                                          PUDataHistName = cms.untracked.string("pileup"),
                                          PFIDCategory = cms.untracked.string("Medium"),
                                          IDMethod = cms.untracked.string("Detector"),
                                          
                                          #Input taken from Ilya's code
                                          rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                          # Objects specific to AOD format
                                          photons = cms.InputTag("gedPhotons"),
                                          genParticles = cms.InputTag("genParticles"),
                                          # Objects specific to MiniAOD format
                                          photonsMiniAOD = cms.InputTag("slimmedPhotons"),
                                          genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                          # ValueMap names from the producer upstream
                                          full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                          phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                          phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                          phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                                          # Locations of files with the effective area constants.
                                          # The constants in these files below are derived for PHYS14 MC.
                                          effAreaChHadFile = cms.FileInPath
                                          #("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
                                          ("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfChargedHadrons_50ns.txt"),
                                          effAreaNeuHadFile= cms.FileInPath
                                          #("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
                                          ("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfNeutralHadrons_50ns.txt"),
                                          effAreaPhoFile   = cms.FileInPath
                                          #("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt"),
                                          ("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfPhotons_50ns.txt"),
                                          # ID decisions (common to all formats)
                                          #phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                                          #phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
                                          #phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight")
                                          phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
                                          phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"),
                                          phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight")
                                          
                                          )
