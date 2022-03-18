////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Particle Properties used in HELIOS from PDG: https://pdg.lbl.gov/2021
//
//
// Axel Drees 8/30/2021
//            9/21/2021 updated
//  
#ifndef PDG_h
#define PDG_h

// constants
  Double_t pi = 3.14159;                              // define pi
  Double_t alpha = 1/137.;                            // fine structure constant

// masses of defined particles
  Double_t pi0Mass  = 0.134977;                        // pi0
  Double_t piMass   = 0.1395702;                       // pi+/pi-
  Double_t etaMass  = 0.54786;                         // eta
  Double_t eMass    = 0.000511;                        // e+/e-
  Double_t muMass   = 0.10566;                         // mu+/mu-
  Double_t rho0Mass = 0.77526;                         // rho meson
  Double_t omegaMass= 0.78266;                         // omega meson
  Double_t etapMass = 0.95778;                         // eta' meson
  Double_t phiMass  = 1.0195;                          // phi meson

  Double_t DeltaMass = 1.232;                          // Delta baryon average
  Double_t NucleonMass = 0.939;                         // Nucleon (p+n)/2

// width of particles
  Double_t omegaWidth = 0.00849;  
  Double_t rho0Width  = 0.1491;   
  Double_t phiWidth   = 0.0042;

//  Double_t DeltaWidth = 0.117;
  
// PDG group particle IDs
  Int_t pi0ID      =  111;
  Int_t pipID      =  211;
  Int_t pimID      = -211;
  Int_t etaID      =  221;
  Int_t photonID   =  22;
  Int_t virtualg   = -22;
  Int_t electronID =  11;
  Int_t positronID = -11;
  Int_t mupID      =  13; 
  Int_t mumID      = -13; 
  Int_t omegaID    = 223;
  Int_t rho0ID     = 113;
  Int_t etapID     = 331;
  Int_t phiID      = 333;

  Int_t DeltaID    = 2114;                            // this is the ID of the Delta_0
  Int_t NucleonID  = 2112;                            // neutron ID  

// Branching ratios of defined decays
  Double_t BR_pi0_gg        = 0.98823; 
  Double_t BR_pi0_Dalitz    = 1-0.98823;
  Double_t BR_eta_gg        = 0.3941;
  Double_t BR_eta_Dalitz    = 0.0069;
  Double_t BR_eta_Dalitz2   = 0.0031;
  Double_t BR_eta_mm        = 0.000058;
  Double_t BR_omega_pi0g    = 0.0834;
  Double_t BR_omega_pi0ee   = 0.00077;
  Double_t BR_omega_ee      = 0.0000739;
  Double_t BR_omega_pi0mm   = 0.000134;
  Double_t BR_omega_mm      = 0.000074;
  Double_t BR_etap_gg       = 0.0231;
  Double_t BR_etap_Dalitz   = 0.000491;
  Double_t BR_etap_Dalitz2  = 0.000113;
  Double_t BR_etap_wg       = 0.0252;
  Double_t BR_etap_rho0g    = 0.295;
  Double_t BR_rho0_ee       = 0.0000472;
  Double_t BR_rho0_mm       = 0.000045;
  Double_t BR_phi_ee        = 0.000297;
  Double_t BR_phi_mm        = 0.000286;
  Double_t BR_Delta_Ng      = 0.006;

// transition formfactor slope values -> NA60 and Lepton G
  Double_t b_pi0_Dalitz = 5.5;
  Double_t b_eta_Dalitz = 1.934;
  Double_t b_etap_Dalitz;
  Double_t b_omega_pi0ll = 2.223; 

// effective temperature defining the line shape of rho -> ll (NA60)
  Double_t T_rho0_ll = 0.161; 

// particle ratios at pt->inft. (can be used as weights for mt-scaled pt distributions)
  Double_t Eta_to_Pi0 = 0.487; // 0.487±0.024 from Yuanjee Ren Thesis
  Double_t Etap_to_Pi0 =  0.25; // 0.25±0.075 from Wenqing Fan Thesis table 4.4
  Double_t Omega_to_Pi0 =  0.81; // 0.81+/-0.02(stat)+/-0.09(sys) from ppg118 
//  Double_t Rho_to_Pi0 = ; 
//  Double_t Phi_to_Pi0 = ; 

#endif 