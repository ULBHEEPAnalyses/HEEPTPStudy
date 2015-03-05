#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  PileUpMerger
  ModifyBeamSpot
  ParticlePropagator
  StatusPid
  GenBeamSpotFilter

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  ModifyBeamSpotNoPU
  ParticlePropagatorNoPU
  ChargedHadronTrackingEfficiencyNoPU
  ElectronTrackingEfficiencyNoPU
  MuonTrackingEfficiencyNoPU
  ChargedHadronMomentumSmearingNoPU
  ElectronEnergySmearingNoPU
  MuonMomentumSmearingNoPU
  TrackMergerNoPU
  CalorimeterNoPU
  EFlowMergerNoPU
  FastJetFinderNoPU

  TrackMerger
  Calorimeter
  TrackPileUpSubtractor
  EFlowMerger

  GlobalRho
  Rho
  FastJetFinder
  GenJetFinder
  JetPileUpSubtractor

  NeutrinoFilter
  GenJetFinderNoNu
  GenMissingET

  EFlowChargedMerger
  RunPUPPI
  PuppiJetFinder
  PuppiRho
  PuppiJetPileUpSubtractor
  PuppiMissingET
  
  PhotonEfficiency
  PhotonIsolation

  ElectronEfficiency
  ElectronIsolation
  
  MuonEfficiency
  MuonIsolation
  
  MissingET

  BTagging
  BTaggingLoose
  TauTagging
  
  
  
  
  
  
  UniqueObjectFinderGJ
  UniqueObjectFinderEJ
  UniqueObjectFinderMJ
  
  ScalarHT
  
  PileUpJetID
  
  PileUpJetIDMissingET
  
  ConstituentFilter
  TreeWriter
}

module Merger PileUpJetIDMissingET {
  add InputArray TrackPileUpSubtractor/eflowTracks
  add InputArray MuonMomentumSmearing/muons
  add InputArray PileUpJetID/eflowTowers
  set MomentumOutputArray momentum
}

module Merger EFlowChargedMerger {
  add InputArray TrackPileUpSubtractor/eflowTracks
  add InputArray MuonMomentumSmearing/muons
  set OutputArray eflowTracks
}

module RunPUPPI RunPUPPI {
#  set TrackInputArray EFlowChargedMerger/eflowTracks
  set TrackInputArray Calorimeter/eflowTracks
  set NeutralInputArray Calorimeter/eflowTowers
  
  set TrackerEta 4.0

  set OutputArray weightedparticles
}

module FastJetFinder PuppiJetFinder {
  set InputArray RunPUPPI/weightedparticles
  set OutputArray jets

  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 0.

  # remove pileup again (using it for synchronization)
  set KeepPileUp 0
}

module FastJetFinder PuppiRho {
  set InputArray RunPUPPI/weightedparticles

  set ComputeRho true
  set RhoOutputArray rho
  
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 5

  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR 0.4
  set GhostEtaMax 5.0
  set RhoEtaMax 5.0
  
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0

  set JetPTMin 0.0
}

module JetPileUpSubtractor PuppiJetPileUpSubtractor {
  set JetInputArray PuppiJetFinder/jets
  set RhoInputArray PuppiRho/rho
  
  set OutputArray jets
  
  set JetPTMin 10.0
}


###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set NPUOutputArray NPU

  # Get rid of beam spot from http://red-gridftp11.unl.edu/Snowmass/MinBias100K_14TeV.pileup ...
  set InputBSX 2.44
  set InputBSY 3.39

  # ... and replace it with beam spot from CMSSW files
  set OutputBSX 0.24
  set OutputBSY 0.39

  # pre-generated minbias input file
  set PileUpFile MinBias.pileup

  # average expected pile up
  set MeanPileUp 50
  # spread in the beam direction in m (assumes gaussian)
  set ZVertexSpread 0.053
}

################
# ModifyBeamSpot
################

module ModifyBeamSpot ModifyBeamSpot {
  set ZVertexSpread 0.053
  set InputArray PileUpMerger/stableParticles
  set OutputArray stableParticles
  set PVOutputArray PV
}

module ModifyBeamSpot ModifyBeamSpotNoPU {
  set ZVertexSpread 0.053
  set InputArray Delphes/stableParticles
  set OutputArray stableParticles
  set PVOutputArray PV
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray PileUpMerger/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.00

  # magnetic field
  set Bz 3.8
}

module ParticlePropagator ParticlePropagatorNoPU {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.00

  # magnetic field
  set Bz 3.8



}


####################################
# StatusPidFilter
# This module removes all generated particles
# except electrons, muons, taus, and status==3
####################################

module StatusPidFilter StatusPid {
#    set InputArray Delphes/stableParticles
    set InputArray Delphes/allParticles
    set OutputArray filteredParticles

    set PTMin 1.0
}

#######################
# GenBeamSpotFilter
# Saves a particle intended to represent the beamspot
#######################

module GenBeamSpotFilter GenBeamSpotFilter {
    set InputArray ModifyBeamSpot/stableParticles
    set OutputArray beamSpotParticles

}



####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt} - Phase II
  set EfficiencyFormula { (pt <= 0.2) * (0.00) + \
(abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
(abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.87) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.82) + \
(abs(eta) > 4.0) * (0.00)}
}

module Efficiency ChargedHadronTrackingEfficiencyNoPU {
  set InputArray ParticlePropagatorNoPU/chargedHadrons
  set OutputArray chargedHadrons
  
    # add EfficiencyFormula {efficiency formula as a function of eta and pt} - Phase II
    set EfficiencyFormula { (pt <= 0.2) * (0.00) + \
(abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
(abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.87) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.82) + \
    (abs(eta) > 4.0) * (0.00)}
}


##############################
# Electron tracking efficiency - ID - Phase-II
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # tracking efficiency formula for electrons

  set EfficiencyFormula { (pt <= 0.2) * (0.00) + \
(abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
(abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0 && pt <= 10.0) * (0.82+pt*0.01) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 10.0) * (0.90) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt <= 10.0) * (0.8+pt*0.01) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0) * (0.85) + \
(abs(eta) > 4.0) * (0.00)}
}

module Efficiency ElectronTrackingEfficiencyNoPU {
  set InputArray ParticlePropagatorNoPU/electrons
  set OutputArray electrons

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # tracking efficiency formula for electrons

    set EfficiencyFormula { (pt <= 0.2) * (0.00) + \
(abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
(abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0 && pt <= 10.0) * (0.82+pt*0.01) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 10.0) * (0.90) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt <= 10.0) * (0.8+pt*0.01) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0) * (0.85) + \
    (abs(eta) > 4.0) * (0.00)}
}


##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons

  set EfficiencyFormula { (pt <= 0.2) * (0.00) + \
(abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.998) + \
(abs(eta) <= 1.2) * (pt > 1.0) * (0.998) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.99) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.99) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + \
(abs(eta) > 4.0) * (0.00)}
}
module Efficiency MuonTrackingEfficiencyNoPU {
  set InputArray ParticlePropagatorNoPU/muons
  set OutputArray muons

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons

    set EfficiencyFormula { (pt <= 0.2) * (0.00) + \
(abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.998) + \
(abs(eta) <= 1.2) * (pt > 1.0) * (0.998) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.99) + \
(abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.99) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + \
(abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + \
    (abs(eta) > 4.0) * (0.00)}
}



########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)}
}

module MomentumSmearing ChargedHadronMomentumSmearingNoPU {
  set InputArray ChargedHadronTrackingEfficiencyNoPU/chargedHadrons
  set OutputArray chargedHadrons

    # set ResolutionFormula {resolution formula as a function of eta and pt}
    set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
						 (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)}
} 



#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}
 set ResolutionFormula { (abs(eta) <= 1.5) * (energy > 0.1   && energy <= 2.5e1) * (energy*0.015) + \
                          (abs(eta) <= 1.5) * (energy > 2.5e1)                    * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
			  (abs(eta) > 1.5 && abs(eta) <= 4.0)                     * sqrt(energy^2*0.008^2 + energy*0.092^2 + 0.088^2)}
}

module EnergySmearing ElectronEnergySmearingNoPU {
  set InputArray ElectronTrackingEfficiencyNoPU/electrons
  set OutputArray electrons

    # set ResolutionFormula {resolution formula as a function of eta and energy}
    set ResolutionFormula { (abs(eta) <= 1.5) * (energy > 0.1   && energy <= 2.5e1) * (energy*0.015) + \
                          (abs(eta) <= 1.5) * (energy > 2.5e1)                    * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
				(abs(eta) > 1.5 && abs(eta) <= 4.0)                     * sqrt(energy^2*0.008^2 + energy*0.092^2 + 0.088^2)}
}


###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for muons

  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)}
}

module MomentumSmearing MuonMomentumSmearingNoPU {
  set InputArray MuonTrackingEfficiencyNoPU/muons
  set OutputArray muons

    # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for muons

    set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
						 (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)}
}


##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  set OutputArray tracks
}

module Merger TrackMergerNoPU {
# add InputArray InputArrxcay
  add InputArray ChargedHadronMomentumSmearingNoPU/chargedHadrons
  add InputArray ElectronEnergySmearingNoPU/electrons
  set OutputArray tracks
}



#############
# Calorimeter
#############

module Calorimeter Calorimeter {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray towers
  set PhotonOutputArray photons

  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowTowers

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 5 degrees towers
  set PhiBins {}
  for {set i -70} {$i <= 70} {incr i} {
    add PhiBins [expr {$i * $pi/70.0}]
  }

  foreach eta {-4 -3.825 -3.65 -3.475 -3.3 -3.125 -2.95 -2.868 -2.65 -2.5 -2.322 -2.172 -2.043 -1.93 -1.83 -1.74 -1.653 -1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.653 1.74 1.83 1.93 2.043 2.172 2.322 2.5 2.65 2.868 2.95 3.125 3.3 3.475 3.65 3.825 4} {
    add EtaPhiBins $eta $PhiBins
  }

  # 10 degrees towers
  set PhiBins {}
  for {set i -11} {$i <= 11} {incr i} {
    add PhiBins [expr {$i * $pi/11.0}]
  }
  foreach eta {-4.35 -4.175 -4 4.175 4.35 4.525} {
    add EtaPhiBins $eta $PhiBins
  }

  # 20 degrees towers
  set PhiBins {}
  for {set i -9} {$i <= 9} {incr i} {
    add PhiBins [expr {$i * $pi/9.0}]
  }
  foreach eta {-5 -4.7 -4.525 4.7 5} {
    add EtaPhiBins $eta $PhiBins
  }

  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {0.0 1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0 0.0}
  add EnergyFraction {22} {1.0 0.0}
  add EnergyFraction {111} {1.0 0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0 0.0}
  add EnergyFraction {13} {0.0 0.0}
  add EnergyFraction {14} {0.0 0.0}
  add EnergyFraction {16} {0.0 0.0}
  add EnergyFraction {1000022} {0.0 0.0}
  add EnergyFraction {1000023} {0.0 0.0}
  add EnergyFraction {1000025} {0.0 0.0}
  add EnergyFraction {1000035} {0.0 0.0}
  add EnergyFraction {1000045} {0.0 0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3 0.7}
  add EnergyFraction {3122} {0.3 0.7}


  # set ECalResolutionFormula {resolution formula as a function of eta and energy}
  set ECalResolutionFormula {(abs(eta) <= 1.5)                   * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
                             (abs(eta) > 1.5 && abs(eta) <= 4.0) * sqrt(energy^2*0.008^2 + energy*0.092^2 + 0.088^2) + \
                             (abs(eta) > 4.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.08^2 + energy*1.97^2)}



  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
 set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) + \
                             (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.050^2 + energy*0.706^2) + \
                             (abs(eta) > 3.0 && abs(eta) <= 4.9) * sqrt(energy^2*0.05^2 + energy*1.00^2)}
}

module Calorimeter CalorimeterNoPU {
  set ParticleInputArray ParticlePropagatorNoPU/stableParticles
  set TrackInputArray TrackMergerNoPU/tracks

  set TowerOutputArray towers
  set PhotonOutputArray photons

  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowTowers

    set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 5 degrees towers
    set PhiBins {}
    for {set i -70} {$i <= 70} {incr i} {
	add PhiBins [expr {$i * $pi/70.0}]
    }

 foreach eta {-4 -3.825 -3.65 -3.475 -3.3 -3.125 -2.95 -2.868 -2.65 -2.5 -2.322 -2.172 -2.043 -1.93 -1.83 -1.74 -1.653 -1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.653 1.74 1.83 1.93 2.043 2.172 2.322 2.5 2.65 2.868 2.95 3.125 3.3 3.475 3.65 3.825 4} {
    add EtaPhiBins $eta $PhiBins
    }
  
  # 10 degrees towers
    set PhiBins {}
    for {set i -11} {$i <= 11} {incr i} {
	add PhiBins [expr {$i * $pi/11.0}]
    }
    foreach eta {-4.35 -4.175 -4 4.175 4.35 4.525} {
    add EtaPhiBins $eta $PhiBins
    }

  # 20 degrees towers
    set PhiBins {}
    for {set i -9} {$i <= 9} {incr i} {
	add PhiBins [expr {$i * $pi/9.0}]
    }
    foreach eta {-5 -4.7 -4.525 4.7 5} {
    add EtaPhiBins $eta $PhiBins
    }

    # default energy fractions {abs(PDG code)} {Fecal Fhcal}
    add EnergyFraction {0} {0.0 1.0}
  # energy fractions for e, gamma and pi0
    add EnergyFraction {11} {1.0 0.0}
    add EnergyFraction {22} {1.0 0.0}
    add EnergyFraction {111} {1.0 0.0}
  # energy fractions for muon, neutrinos and neutralinos
    add EnergyFraction {12} {0.0 0.0}
    add EnergyFraction {13} {0.0 0.0}
    add EnergyFraction {14} {0.0 0.0}
    add EnergyFraction {16} {0.0 0.0}
    add EnergyFraction {1000022} {0.0 0.0}
    add EnergyFraction {1000023} {0.0 0.0}
    add EnergyFraction {1000025} {0.0 0.0}
    add EnergyFraction {1000035} {0.0 0.0}
    add EnergyFraction {1000045} {0.0 0.0}
  # energy fractions for K0short and Lambda
    add EnergyFraction {310} {0.3 0.7}
    add EnergyFraction {3122} {0.3 0.7}

    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    set ECalResolutionFormula {(abs(eta) <= 1.5)                   * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
                             (abs(eta) > 1.5 && abs(eta) <= 4.0) * sqrt(energy^2*0.008^2 + energy*0.092^2 + 0.088^2) + \
				   (abs(eta) > 4.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.08^2 + energy*1.97^2)}


  
    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) + \
                             (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.050^2 + energy*0.706^2) + \
						     (abs(eta) > 3.0 && abs(eta) <= 4.9) * sqrt(energy^2*0.05^2 + energy*1.00^2)}
}


##########################
# Track pile-up subtractor
##########################

module TrackPileUpSubtractor TrackPileUpSubtractor {
# add InputArray InputArray OutputArray
  add InputArray Calorimeter/eflowTracks eflowTracks
  add InputArray ElectronEnergySmearing/electrons electrons
  add InputArray MuonMomentumSmearing/muons muons

  set PVInputArray  ModifyBeamSpot/PV

  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution
  # Z vertex resolution in m
  set ZVertexResolution 0.0001
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray TrackPileUpSubtractor/eflowTracks
  add InputArray Calorimeter/eflowTowers
  add InputArray MuonMomentumSmearing/muons
  set OutputArray eflow
}

module Merger EFlowMergerNoPU {
# add InputArray InputArray
  add InputArray CalorimeterNoPU/eflowTracks
  add InputArray CalorimeterNoPU/eflowTowers
  add InputArray MuonMomentumSmearingNoPU/muons
  set OutputArray eflow
}


#############
# Rho pile-up
#############

module FastJetFinder Rho {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set ComputeRho true
  set RhoOutputArray rho

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 5

  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR 0.4
  set GhostEtaMax 5.0
  set RhoEtaMax 5.0

  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0

  set JetPTMin 0.0
}

module FastJetFinder GlobalRho {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set ComputeRho true
  set RhoOutputArray rho

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 5
  
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR 0.4
  set GhostEtaMax 5.0
  set RhoEtaMax 5.0
  
  add RhoEtaRange 0.0 5.0

  set JetPTMin 0.0
}


#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray Delphes/stableParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 5.0

}

module NeutrinoFilter NeutrinoFilter {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles  
}

module FastJetFinder GenJetFinderNoNu {
  set InputArray NeutrinoFilter/stableParticles
  
  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4
  
  set JetPTMin 5.0

}


############
# Jet finder
############

module FastJetFinder FastJetFinder {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 5

  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 5.0
}

module FastJetFinder FastJetFinderNoPU {
#  set InputArray CalorimeterNoPU/towers
  set InputArray EFlowMergerNoPU/eflow

  set OutputArray jets

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 5

  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 5.0
}



############
# Cambridge-Aachen Jet finder
############

module FastJetFinder CAJetFinder {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow
  set OutputArray jets
  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set AreaAlgorithm 5
  set JetAlgorithm 5
  set ParameterR 0.8
  # 200 GeV needed for boosted W bosons, 300 GeV is safe for boosted tops
  set JetPTMin 200.0
}

####################
# Constituent filter
####################

module ConstituentFilter ConstituentFilter {

  set ConEMin 1.

# # add JetInputArray InputArray
   add JetInputArray GenJetFinderNoNu/jets

# SZ changed this but it seems sensible
#   add JetInputArray FastJetFinder/jets
   add JetInputArray UniqueObjectFinderMJ/jets

#   add JetInputArray CAJetFinder/jets


# # add ConstituentInputArray InputArray OutputArray
   add ConstituentInputArray Delphes/stableParticles stableParticles
   add ConstituentInputArray TrackPileUpSubtractor/eflowTracks eflowTracks
   add ConstituentInputArray Calorimeter/eflowTowers eflowTowers
   add ConstituentInputArray MuonMomentumSmearing/muons muons
  # }



###########################
# Jet Pile-Up Subtraction
###########################

module JetPileUpSubtractor JetPileUpSubtractor {
  set JetInputArray FastJetFinder/jets
  set RhoInputArray Rho/rho

  set OutputArray jets

  set JetPTMin 5.0
}

module JetPileUpSubtractor CAJetPileUpSubtractor {
  set JetInputArray CAJetFinder/jets
  set RhoInputArray Rho/rho
  set OutputArray jets
  set JetPTMin 20.0
}


###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray Calorimeter/photons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.9635) + \
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 10.0)  * (0.9624) + \
                         (abs(eta) > 4.0)                                   * (0.00)}
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowMerger/eflow
  set RhoInputArray Rho/rho

  set OutputArray photons

  set DeltaRMax 0.3

  set PTMin 1.0

  set PTRatioMax 0.4
}

#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray TrackPileUpSubtractor/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
    set EfficiencyFormula {                                      (pt <= 4.0)  * (0.00) + \
                         (abs(eta) <= 1.45 ) * (pt >  4.0 && pt <= 6.0)   * (0.50) + \
                         (abs(eta) <= 1.45 ) * (pt >  6.0 && pt <= 8.0)   * (0.70) + \
                         (abs(eta) <= 1.45 ) * (pt >  8.0 && pt <= 10.0)  * (0.85) + \
                         (abs(eta) <= 1.45 ) * (pt > 10.0 && pt <= 30.0)  * (0.94) + \                                                      
                         (abs(eta) <= 1.45 ) * (pt > 30.0 && pt <= 50.0)  * (0.97) + \                          
                         (abs(eta) <= 1.45 ) * (pt > 50.0 && pt <= 70.0)  * (0.98) + \          
                         (abs(eta) <= 1.45 ) * (pt > 70.0 )  * (1.0) + \                                                                                                                 
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt >  4.0 && pt <= 10.0)   * (0.35) + \
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 10.0 && pt <= 30.0)   * (0.40) + \   
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 30.0 && pt <= 70.0)   * (0.45) + \                                 
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 70.0 )  * (0.55) + \    
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt >  4.0 && pt <= 10.0)  * (0.75) + \
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 10.0 && pt <= 30.0)  * (0.85) + \                                                      
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 30.0 && pt <= 50.0)  * (0.95) + \                          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 50.0 && pt <= 70.0)  * (0.95) + \          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 70.0 )  * (1.0) + \   
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt >  4.0 && pt <= 10.0)  * (0.65) + \
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 10.0 && pt <= 30.0)  * (0.75) + \                                                      
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 30.0 && pt <= 50.0)  * (0.90) + \                          
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 50.0 && pt <= 70.0)  * (0.90) + \          
                         (abs(eta) >= 2.0 && abs(eta) <= 4.0 ) * (pt > 70.0 )  * (0.90) + \                                                                                               
                          (abs(eta) > 4.0)                              * (0.00)}

}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowMerger/eflow
  set RhoInputArray Rho/rho

  set OutputArray electrons

  set DeltaRMax 0.3

  set PTMin 1.0

  set PTRatioMax 0.4
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray TrackPileUpSubtractor/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
    set EfficiencyFormula {                                    (pt <= 2.0)  * (0.00) + \  
                         (abs(eta) <= 4.00) * (pt >  2.0 && pt <= 3.0)  * (0.51) + \
                         (abs(eta) <= 4.00) * (pt >  3.0 && pt <= 4.0)  * (0.85) + \ 
                         (abs(eta) <= 4.00) * (pt >  4.0 && pt <= 11.0) * (0.93) + \               
                         (abs(eta) <= 4.00) * (pt >  11. && pt <= 50.)  * (0.96) + \   
                         (abs(eta) <= 4.00) * (pt >  50. && pt <= 70.)  * (0.98) + \                      
                         (abs(eta) <= 4.00) * (pt > 70.0 )  * (1.00) + \   
                         (abs(eta) > 4.00)  * (0.00)}

}

################
# Muon isolation
################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set IsolationInputArray EFlowMerger/eflow
  set RhoInputArray Rho/rho

  set OutputArray muons

  set DeltaRMax 0.3

  set PTMin 1.0

  set PTRatioMax 0.4
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

module Merger GenMissingET {
#  add InputArray Delphes/stableParticles
  add InputArray NeutrinoFilter/stableParticles
  set MomentumOutputArray momentum
}

module Merger PuppiMissingET {
  add InputArray RunPUPPI/weightedparticles
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################
module Merger ScalarHT {
# add InputArray InputArray
  add InputArray UniqueObjectFinderMJ/jets
  add InputArray UniqueObjectFinderEJ/electrons
  add InputArray UniqueObjectFinderGJ/photons
  add InputArray UniqueObjectFinderMJ/muons
  set EnergyOutputArray energy
}


###########
# b-tagging
###########

module BTagging BTagging {
  set PartonInputArray Delphes/partons
#  set JetInputArray FastJetFinder/jets
  set JetInputArray JetPileUpSubtractor/jets

  set BitNumber 0
  set DeltaR 0.4
  set PartonPTMin 1.0
  set PartonEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsBTV
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}

  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {                                      (pt <= 15.0) * (0.000) + \
                                                (abs(eta) <= 1.2) * (pt > 15.0) * (0.1873*tanh(pt*0.0183 - 0.2196)) + \
                              (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.1898*tanh(pt*0.00997 - 0.143)) + \
                              (abs(eta) > 4.0)                                  * (0.000)}

  # efficiency formula for b-jets
    add EfficiencyFormula {5} {                                      (pt <= 15.0) * (0.000) + \
                                                (abs(eta) <= 1.2) * (pt > 15.0) * (0.629858*tanh(pt*0.0166188 + 0.300119)) + \
                              (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.584522*tanh(pt*0.0144387 + 0.397034)) + \
									 (abs(eta) > 4.0)                                  * (0.000)}
}

module BTagging BTaggingLoose {
  set PartonInputArray Delphes/partons
#  set JetInputArray FastJetFinder/jets
  set JetInputArray JetPileUpSubtractor/jets

  set BitNumber 1
  set DeltaR 0.4
  set PartonPTMin 1.0
  set PartonEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsBTV
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.02}

  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {                                      (pt <= 15.0) * (0.000) + \
                                                (abs(eta) <= 1.2) * (pt > 15.0) * (0.29*tanh(pt*0.0183 - 0.2196)) + \
                              (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.29*tanh(pt*0.00997 - 0.143)) + \
                              (abs(eta) > 4.0)                                  * (0.000)}

  # efficiency formula for b-jets

    add EfficiencyFormula {5} {                                      (pt <= 15.0) * (0.000) + \
                                                (abs(eta) <= 1.2) * (pt > 15.0) * (0.629858*tanh(pt*0.0166188 + 0.300119)) + \
                              (abs(eta) > 1.2 && abs(eta) <= 4.0) * (pt > 15.0) * (0.584522*tanh(pt*0.0144387 + 0.397034)) + \
									 (abs(eta) > 4.0)                                  * (0.000)}

}


module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
#  set JetInputArray FastJetFinder/jets
  set JetInputArray JetPileUpSubtractor/jets

  set DeltaR 0.4

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.004}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.65}
}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

#module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
#  add InputArray PhotonIsolation/photons photons
#  add InputArray ElectronIsolation/electrons electrons
#  add InputArray JetPileUpSubtractor/jets jets
#}

module UniqueObjectFinder UniqueObjectFinderGJ {
   add InputArray PhotonIsolation/photons photons
   add InputArray JetPileUpSubtractor/jets jets
}

module UniqueObjectFinder UniqueObjectFinderEJ {
   add InputArray ElectronIsolation/electrons electrons
   add InputArray UniqueObjectFinderGJ/jets jets
}

module UniqueObjectFinder UniqueObjectFinderMJ {
   add InputArray MuonIsolation/muons muons
   add InputArray UniqueObjectFinderEJ/jets jets
}

### 
#Pileup jet id
###

module PileUpJetID PileUpJetID {
  set JetInputArray JetPileUpSubtractor/jets
  set OutputArray jets
  set NeutralsInPassingJets eflowTowers

  # Using constituents does not make sense with Charged hadron subtraction                                                                                                           
  # In 0 mode, dR cut used instead                                                                                                                                                   
  set UseConstituents 0

  set TrackInputArray Calorimeter/eflowTracks
  set NeutralInputArray Calorimeter/eflowTowers
  set ParameterR 0.4

  set JetPTMin 5.0

  set MeanSqDeltaRMaxBarrel 0.13
  set BetaMinBarrel 0.16
  set MeanSqDeltaRMaxEndcap 0.07
  set BetaMinEndcap 0.06
  set MeanSqDeltaRMaxForward 0.01

}



##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
  add Branch StatusPid/filteredParticles Particle GenParticle
  add Branch GenBeamSpotFilter/beamSpotParticles BeamSpotParticle GenParticle

  add Branch FastJetFinder/jets RawJet Jet
  add Branch FastJetFinderNoPU/jets RawJetNoPU Jet

  add Branch GenJetFinder/jets GenJetWithNu Jet
  add Branch GenJetFinderNoNu/jets GenJet Jet
  add Branch UniqueObjectFinderMJ/jets Jet Jet
  add Branch UniqueObjectFinderEJ/electrons Electron Electron
  add Branch UniqueObjectFinderGJ/photons Photon Photon
  add Branch UniqueObjectFinderMJ/muons Muon Muon

  add Branch PileUpJetIDMissingET/momentum PileUpJetIDMissingET MissingET
  add Branch GenMissingET/momentum GenMissingET MissingET
  add Branch PuppiMissingET/momentum PuppiMissingET MissingET


  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
  add Branch Rho/rho Rho Rho
  add Branch GlobalRho/rho GlobalRho Rho
  add Branch PileUpMerger/NPU NPU ScalarHT

  add Branch PuppiJetFinder/jets PuppiJet Jet

  set OffsetFromModifyBeamSpot 0

#  add Branch RunPUPPI/weightedparticles PuppiWeightedParticles GenParticle
#  add Branch Delphes/allParticles Particle GenParticle
#  add Branch Calorimeter/eflowTracks EFlowTrack Track
#  add Branch Calorimeter/eflowTowers EFlowTower Tower
#  add Branch MuonMomentumSmearing/muons EFlowMuon Muon
#  add Branch PuppiJetPileUpSubtractor/jets SubtractedPuppiJet Jet
#  add Branch PuppiRho/rho PuppiRho Rho
}

# # add Branch InputArray BranchName BranchClass
#  # add Branch Delphes/allParticles Particle GenParticle
  # add Branch StatusPid/filteredParticles Particle GenParticle
#  # add Branch TrackMerger/tracks Track Track
#  # add Branch Calorimeter/towers Tower Tower
#  # add Branch ConstituentFilter/eflowTracks EFlowTrack Track
#  # add Branch ConstituentFilter/eflowTowers EFlowTower Tower
#  # add Branch ConstituentFilter/muons EFlowMuon Muon
  # add Branch GenJetFinder/jets GenJet Jet
  # add Branch CAJetPileUpSubtractor/jets CAJet Jet
  # add Branch UniqueObjectFinderMJ/jets Jet Jet
  # add Branch UniqueObjectFinderEJ/electrons Electron Electron
  # add Branch UniqueObjectFinderGJ/photons Photon Photon
  # add Branch UniqueObjectFinderMJ/muons Muon Muon
  # add Branch MissingET/momentum MissingET MissingET
  # add Branch ScalarHT/energy ScalarHT ScalarHT
  # add Branch Rho/rho Rho ScalarHT



