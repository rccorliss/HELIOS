HELIOS 

Axel Drees 11/16/2021

Standalone fast simulation package that run in ROOT framework (version 6). 
Desigined to aid systematic studies of measurements of photons and lepton 
pairs in high energy p+p, p+A and A+A collisions, in particular with PHENIX. 

HELIOS package components in /HELIOSLibrary/ 

include file
- HELIOSLibrary.h             use #include "...path../HELIOSLibrary.h" to include all necessary 
                              components of HELIOS in your code
Particle generator
- PDG.h                       Collection of particle properties and constants
- KinematicDistributions.h    Collection of momentum spectra of various particles
- Particle.C, Particle.h      Particle Class - inherets from TLorentzVector, 
                              defines, particles and decays, includes 3-momentum generator 
                              and generates decay particles using decay class 
- Decay.C, Decay.h            Decay Class, used by Particle to generate decays
- DecayAuxilliaryFunctions.h  collection of functions used in Decay Class 

Interaction with matter
- InteractionWithMaterial.h   simulation of physics processes like Photon conversions, bremmsstrahlung etc. 

Experiment specific simulation
- PHENIXSetup.h               collection of constants defining PHENIX
- MyPHENIX.h                  PHENIX detector class

////// HELIOSLibrary.h /////////////////////////////////////////////////////////////////////////////////////
//
//  includes all necessary code for the HELIOS packackage and 
//  defines string with full pathlength to local directory of HELIOS package
//  used for some data files
//  
//  this string need to be updated for installation on new computer 
//  e.g.  C:/root_v6.22.06/macros/HELIOS/
//
//  Axel Drees 11/15/2021
//
//
///// PDG.h /////////////////////////////////////////////////////////////////////////////////////////////////
//
// Collection Particle Properties used in HELIOS from PDG: https://pdg.lbl.gov/2021
//
//
// Axel Drees 8/30/2021
//            9/21/2021 updated
//  
//// Particle.h //////////////////////////////////////////////////////////////////////////////////////////////
//
// HELIOS Particle Class used for stand alone simulation in ROOT 6
//
// inherits from TLorentzVector
//
// Particle(name)  - Defines particles as Lorentz vectors, with parameters define in PDG.h. 
//                 - Defines decay branches as particle property if particle is not stable. 
//                 - Decays are handeled in Class Decay.C 
//
// currently implemented particles and decay branches by name used in HELIOS
// 
//     hadrons       decays 
//     pi0           pi0->gg, pi0->gee
//     pi+           stable
//     pi-           stable
//     eta           eta->gg, eta->gee
//     etap          etap->gg, etap->gee
//     rho0          rho0->ee
//     omega         omega->ee, omega->pi0g, omega->pi0ee     
//     phi           phi->ee
//
//     leptons     
//     photon        stable
//     electron      stable
//     positron      stable
//
// generate random 3 vector 
// GenerateP()     - there are three options implemented through over loading 
//                 - flat, from pt spectrum, or pt, eta and phi spectra (see below)
//
// decays particle, three different options 
// Decay()         - generate one random decay with probability equal to branching ratio
// DecaySingleBranch(name) - forces a particular decay
// DeacyFlat        - generates one random decay with equal probability for all defined branches
//
// member variables
//
// Double_t weight - allows to assign a weight to a given particle, e.g. dN/dpt for p generated flat
//                 - Weight() returns value
// Int_t PID       - PDG particle id, ID() returns value 
// TString Name    - name used by HELIOS - defined in PDG, Name() returns value 
// Int_t Charge    - charge of particle - Charge() returns value
//
// member functions:
//
// SetWeight(w)       - sets value of weight
// AddWeight(w)       - multiplies w to existing weight
// ResetP()           - set Lorentz vector to (0,0,0,m)
// UpdateParticle(p)  - sets TLorentz vector of particle to p
// 
// access to decay particle information:
//
//  Int_t GetNumberOfDaughters()    - returns number of daugthers generated for this particle
//  Int_t GetDaughterID(i)          - returns PDG ID of daughter i 
//  Double_t GetDaughterWeight(i)   - returns get weight for daughter, depending on method used for decay generator  
//  TLorentzVector GetDecayDaughter(i) - returns 4 vector of daughter
// 
//  
// Axel Drees 11/19/2019
// integrated to HELIOS 9/21/2021
//
////// Decay.h //////////////////////////////////////////////////////////////////////////////////////////////////
//
// Header for Class Decay that handels all decay branches defined for particles in HELIOS
//
// currently defined decay names (selfexplanatory) used to identify the correct parent specific decay branch 
//
//     "pi0->gg"       "pi0->gee"
//     "eta->gg"       "eta->gee"
//     "etap->gg"      "etap->gee"        "etap->rho0g"
//     "omega->pi0g"   "omega->pi0ee"     "omega->ee"
//     "rho0->ee"
//     "phi->ee" 
// 
// currently defined decay types of decay branches used internaly to select correct decay function                                   
//     "TwoBody"  two body decay 
//     "Dalitz"   Dalitz decay
//
// Particle class defines each decay branch for non stable particle with Decay constructor and sets 
// properties using:
//   SetName("name")                   as defined above
//   SetType("type")                   as defined above
//   SetBR(BR)                         branching ratio of decay 
//   SetNumberOfDecayParticles(n)      # of decay particles
//   SetMassDistributions()            of lepton pairs from Dalitz or vector meson decays 
//   DefineDaughters(*4vector)         preset array of TLorentz vectors of decay particles  
//
// random decay particle generation for defined branch:
//   GenerateDecay()                   produces random decay particles for this branch boosted in to parent frame
//   SetParent(pt,eta,phi,m)           used to update parent particle prior to decay simulation
//
// decay particles and properties can be accessed via:
//  GetNumberDecayParticles()    returns number of decay particles for this decay branch
//  GetDaughter(i)               returns TLorentzVector of daughter i    
//  GetBR()                      returns branching ratio
//  GetName()                    returns name of decay branch
//  GetType()                    returns type of decay branch
//
// defined operators:    " = "   set one decay branch equal to another
//
//
// originally developed Axel Drees 8/23/2021
// addpted to HELIOS               9/21/2021
//
////// DecayAuxilliaryFunctions.h ///////////////////////////////////////////////////////////////////////////////////////
//
// Auxiliary functions used in Class Decays  
//
// mass distributions functions use to define TF1's for Dalitz Decays: 
//
// f_pi0Dalitz 
// f_etaDalitz 
// f_etapDalitz 
// f_omegaDalitz 
//
// mass distributions functions used to define TF1's for 2 body decays with finite width:
//
// f_omegaee
// f_rho0ee
// f_phiee
//
// Support functions:
// 
// KrollWada:        is a generalized Kroll Wada Function for p -> o e+e-, where p is the parent 
//                   and o is the "other" particle. For o being a photon this will be the well 
//                   known Kroll Wada distribution. 
//
// BreitWigner:      Breit-Wigner parameterisation of resonances.
//
// GounarisSakurai:  Generalized Breit-Wigner function, commonly used for all 2 body decays 
//                   with finite width. 
//
// EMTFormFactor:    Lepton G parameterization of Electro Magnetic Tranition Form Factor in 
//                   Dalitz decays.
//
// Axel Drees 8/30/2021
// expanded   9/9/2021    
//
//// InteractionWithMaterial.h //////////////////////////////////////////////////////////////////////////
//
// collection of functions to simulate interaction of particles with material
//
// currently implemented:
//
//   PhotonConversion
//   
//
// Axel Drees  9/22/2021 
//
/////// PHENIXSetup.h //////////////////////////////////////////////////////////////////////////////////////////
//
// Collection of constants defining PHENIX Setup 
// used in Class PHENIX 
//
// Axel Drees 9/14/2021
//  