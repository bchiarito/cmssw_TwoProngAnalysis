import FWCore.ParameterSet.Config as cms

diphotonFilter = cms.EDFilter('DiPhotonFilter',
                              photonCollection = cms.untracked.InputTag("photons"),
                              ptMin_photon1 = cms.untracked.double(20),
                              ptMin_photon2 = cms.untracked.double(20)


)
