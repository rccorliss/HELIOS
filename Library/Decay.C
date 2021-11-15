//// Decay.C /////////////////////////////////////////////////////////////////////////////////////////////////
//
// see Decay.h for detailed description
//
// handels decay branches defined as particle properties in HELIOS 
// 
//
// Axel Drees 11/19/2019
// addaped to HELIOS 9/21/2021
//  
//
//

#include "Decay.h"

Decay::Decay(){             // Default constructor
  BR = 0;
  DecayType = "unkown";
  DecayName = "unkown";
  NumberDecayParticles = 0; 
}
Decay::~Decay(void){}      // destructor

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::GenerateDecay(){
//
// generate random decay particles depending on the defined decay type
//
  TLorentzVector parent;
  TLorentzVector daughter1;
  TLorentzVector daughter2;
  TLorentzVector daughter3;
  Double_t pt, eta, phi, eemass; 

  if (DecayType == "TwoBody") {
//    std::cout << "Decay name " << DecayName << std::endl; 
    if (DecayName == "omega->ee" or DecayName == "rho0->ee" 
        or DecayName == "phi->ee") {
       TF1 *h = (TF1 *)gROOT->FindObject(DecayName);   // get eepair mass
       eemass = h->GetRandom(); 
       pt  = Parent.Pt();
       eta = Parent.Eta();
       phi = Parent.Phi();
//       std::cout << " ----- " << eemass  << std:: endl; 
       parent.SetPtEtaPhiM(pt,eta,phi,eemass);
    } else {
      parent = Parent;
    }
    daughter1.SetPtEtaPhiM(0.,0.,0.,Daughter[0].M());
    daughter2.SetPtEtaPhiM(0.,0.,0.,Daughter[1].M());
    TwoBodyDecay(parent,daughter1,daughter2);
    Daughter[0] = daughter1;
    Daughter[1] = daughter2;
  } else if (DecayType == "Dalitz" ){
    parent = Parent;
    daughter1.SetPtEtaPhiM(0.,0.,0.,Daughter[0].M());
    daughter2.SetPtEtaPhiM(0.,0.,0.,Daughter[1].M());
    daughter3.SetPtEtaPhiM(0.,0.,0.,Daughter[2].M());
    DalitzDecay(parent,daughter1,daughter2,daughter3);
    Daughter[0] = daughter1;
    Daughter[1] = daughter2;
    Daughter[2] = daughter3;
  } else {
    BR = 0;
    NumberDecayParticles=0;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::TwoBodyDecay(TLorentzVector &parent, TLorentzVector &daughter1, TLorentzVector &daughter2){
//
// Isotropic two body decay of parent into two particles
//
// input:  parent 4 vector
// output: daughter1 &2 4 vector 
//
// Axel Drees 10/17/2018 - checked to work for 2 photon decays
//            8/29/2021 - adapted from standalone function to this class
//
  TRandom3 randy = TRandom3(0);                // Random Generator
  Double_t pi = 3.14159;                              // define pi

  Double_t px,py,pz,E;                                // define generic 4 vector component 
  Double_t mass = parent.M();                         // mass of parent particle
  Double_t m1 = daughter1.M();                      // mass of decay particle 1 
  Double_t m2 = daughter2.M();                      // mass of decay particle 2
  TVector3 boost = parent.BoostVector();              // boost vector of parent particle

//
// decay in rest frame 
//
  Double_t p     = sqrt((mass*mass-(m1+m2)*(m1+m2))*(mass*mass-(m1-m2)*(m1-m2)))/2/mass; 
//  std::cout << "----- p --------- " << p << std::endl;
  Double_t phi = randy.Uniform(0.,2*pi);            // phi random from 0. to 2pi
  Double_t z   = randy.Uniform(-1.,1.);             // z used to calculate theta 
//  std::cout << "----- phi,z --------- " << phi << " " << z << std::endl;
  Double_t theta = acos(z);                           //     so that each angle in 4pi is equaly likely
  pz = p*cos(theta);                                  // longitudinal momentum
  px = p*sin(theta)*cos(phi);                         // transverse momentum x component
  py = p*sin(theta)*sin(phi);                         // transverse momentum y component
  E = sqrt(px*px+py*py+pz*pz+m1*m1);                  // energy of decay particle 1
  daughter1.SetPxPyPzE(px,py,pz,E);                 // set 4 vector of 1 decay particle
  E = sqrt(px*px+py*py+pz*pz+m2*m2);                  // energy of decay particle 1
  daughter2.SetPxPyPzE(-px,-py,-pz,E);              // set 4 vector of 2 decay particle in opposit direction
//  std::cout << "-----------------" << E << std::endl;
//
// boost to parent frame
//
  daughter1.Boost(boost);                           // boost with parent momentum
  daughter2.Boost(boost);                           // boost with parent momentum
//  std::cout << "---------------- " <<Daughter[0].E() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::DalitzDecay(TLorentzVector &parent, TLorentzVector &decay1, TLorentzVector &decay2, TLorentzVector &decay3){
//
// Isotropic two body decay of parent into photon and virtual photon, followed by two body decay of virtual photon 
//
// input:  parent              4 vector
// output: decay particle 3    4 vector must have E = mass p=0
//         decay particle 1&2  4 vector of electron and positron, must have E = mass electron and p=0
//         MassDistribution    - function with mass distribution for Dalitz ee pair
// 
// output: parent and decay    4 vectors
//
// The decay is calculated in several seps
//
//   1. determine mass of ee pair using MassDistribution
//
// Axel Drees 11/16/2019 - 
//
  Double_t pi = 3.14159;                              // define pi
  Double_t px,py,pz,E;                                // define generic 4 vector component 
  Double_t eemass=0;
  Double_t mass = parent.M();                         // mass of parent particle
  Double_t m1 = decay1.M();                           // mass of decay particle 1 
  Double_t m2 = decay2.M();                           // mass of decay particle 2
  Double_t m3 = decay3.M();                           // mass of decay particle 2
  TLorentzVector eepair;
  TVector3 parentBoost = parent.BoostVector();        // boost vector of parent particle
//
// mass of ee pair
  TF1 *h = (TF1 *)gROOT->FindObject(DecayName);           // get eepair mass
  eemass = h->GetRandom(); 
  eepair.SetPxPyPzE(0.,0.,0.,eemass);                 // set eepair lorentz vector

  // std::cout << std::endl;
  // std::cout << "ee par mass: " << eemass << std::endl;      
//
// first 2 body decay parent->decay3 eepair
  TwoBodyDecay(parent,decay3,eepair);
//
// second 2 body decay of eepair to decay1 and decay2
  TwoBodyDecay(eepair,decay1,decay2);

  
  // parent.Print();
  // decay3.Print();                   
  // decay1.Print();
  // decay2.Print();
  // eepair.Print();                   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Decay::SetMassDistributions(){
//
// define TF1's used to sample mass for various decays including
//
//     - virtual photon mass from Dalitz Decays
//     - short lived resonances omega, phi J/psi etc
//     - longer lived rho meson
//
//
//
//   
  Double_t c1 = 2*alpha/3/pi;                         // constant factor
  
  if (DecayName == "pi0->gee") {
      TF1 *h_Mass = new TF1(DecayName,f_pi0Dalitz,2.*eMass,pi0Mass);
  } else if (DecayName == "eta->gee") {
      TF1 *h_Mass = new TF1(DecayName,f_etaDalitz,2.*eMass,etaMass);
  } else if (DecayName == "etap->gee") {
      TF1 *h_Mass = new TF1(DecayName,f_etapDalitz,2.*eMass,etapMass);
  } else if (DecayName == "etap->rho0g") {
      TF1 *h_Mass = new TF1(DecayName,f_etapRho0g,2.*eMass,etapMass);
  } else if (DecayName == "omega->pi0ee") {
      TF1 *h_Mass = new TF1(DecayName,f_omegaDalitz,2.*eMass,omegaMass-pi0Mass);
  } else if (DecayName == "omega->ee") {
      TF1 *h_Mass = new TF1(DecayName,f_omegaee,2.*pi0Mass,2*omegaMass); 
  } else if (DecayName == "rho0->ee") {
      TF1 *h_Mass = new TF1(DecayName,f_rho0ee,2.*pi0Mass,2*rho0Mass); 
  } else if (DecayName == "phi->ee") {
      TF1 *h_Mass = new TF1(DecayName,f_phiee,2.*pi0Mass,2*phiMass); 
  }


}

