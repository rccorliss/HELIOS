//////   PHENIXDetector.h //////////////////////////////////////////////////////////////////
//
// This Class containing a collection of PHENIX specific functions 
// describing the EMCal and charged particle tracking
//
// changed from collection of functions to class structure 9/13/2021
//  
//   InAcceptance(particle,charge)                        Ideal PHENIX acceptance - checks if particle is in acceptance 
//   EMCalSector(phi)                                     EMCal sector number 1 - 8 (phi only)
//   SmearEnergy(randomGenerator,energy,option,c1,c2)     EMCal energy resoltion function - parameters c1,c2 are optional
//   NonLinearEnergy(energy,constant1,constant2)          parameterized non linear energy respose of EMCal
//
//
// Axel Drees
// updated 11/16/2021
//
///////////////////////////////////////////////////////////////////////////////////

#ifndef PHENIX_h
#define PHENIX_h

#include <TLorentzVector.h>                                    // MyParticle inherits from TLorentz Vector
#include <TRandom3.h>
#include <TMath.h>
#include "PDG.h"                                             // 
#include "PHENIXSetup.h"
#include "InteractionWithMaterial.h"

class PHENIXDetector
{

public:
  PHENIXDetector();                                                // default constructor
  virtual ~PHENIXDetector();                                       // destructor

// member functions - see definitions below for details 
// general central arm functions  
  Int_t InAcceptance(TLorentzVector particle, Int_t q); 

// calculate inportant variables for given particle
  void CharacterizeTrack(TLorentzVector VT, Int_t q, Int_t id);

// EMCal 
  TLorentzVector ReconstructShower(TLorentzVector VT, Int_t id);
  Int_t    EMCalSector(Double_t phi);
  Int_t    EMCalSectorCoordinates(Double_t phi, Double_t theta, Double_t& y, Double_t& z);
  Bool_t   EMCalLive(Int_t sector);
  Double_t NonLinearEnergy(Double_t energy, Double_t c0=1.003, Double_t c1=0.05, Double_t c2=1.77);
  Double_t SmearEnergy(Double_t energy, Int_t opt=0, Double_t c1=0.081, Double_t c2=0.021 );
  Double_t SmearEMCalPosition(Double_t energy, Double_t x, Int_t opt=0 );
  Double_t ShowerEnergyCorrection(Double_t energy, Double_t sinT);
  Double_t ShowerEnergyCorrection(Double_t energy, Double_t y, Double_t z);
  Double_t EMCalImpactAngle(TLorentzVector VT, Int_t id);
  Double_t EMCalPhi(TLorentzVector VT, Double_t q); 

// Tracking
  TLorentzVector ReconstructTrack(TLorentzVector VT, Int_t id);
  Double_t SmearPt(Double_t pt, Int_t &q);
  Double_t SmearDCphi0(Double_t pt, Double_t phi0);
  Double_t SmearDCeta(Double_t pt, Double_t eta0);
  Double_t DCPhi(TLorentzVector VT, Double_t q); 
  Double_t RICHPhi(TLorentzVector VT, Double_t q); 

  Double_t ElectronMomentum(Double_t E, Double_t p);

// physics processes 
  Int_t VTXConversion();

// get Private variables
  Int_t Arm(){
    return arm;
  }
  Int_t Sector(){
    return sector;
  }
  Double_t SectorY(){
    return sectorY;
  }
  Double_t SectorZ(){
    return sectorZ;
  }
  Double_t SectorSinT(){
//   if (sector != 0) std::cout << " --+-- " << sectorSinT << std::endl;
    return sectorSinT;
  }
  Double_t Phi_EMCal(){
    return phi_EMCal;
  }
  Double_t Phi_DC(){
    return phi_DC;
  }
  Double_t Phi_RICH(){
    return phi_RICH;
  }

private:

// internal variables
  TRandom3 randy = TRandom3(0);                                // Random Generator
  Double_t SectorX0[8],SectorXmax[8],SectorXmin[8];
  Double_t SectorY0[8],SectorYmax[8],SectorYmin[8];
  Double_t P_R_EMCal_max[8];
  Bool_t Intersect(Double_t &Px,Double_t &Py,Double_t *x, Double_t *y);

  Double_t  phi_EMCal;
  Double_t  phi_RICH;
  Double_t  phi_DC;
  Double_t  alpha_DC;
  Double_t  alpha_EMCal;
  Double_t  sectorZ;
  Double_t  sectorY;
  Double_t  sectorSinT;
  Int_t     sector;
  Int_t     arm; 

};

#endif

