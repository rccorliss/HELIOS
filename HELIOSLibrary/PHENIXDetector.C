////// PHENIXDetector.C //////////////////////////////////////////////////////////////////
//
// This Class containing a collection of PHENIX specific functions 
// describing the EMCal and charged particle tracking
//
// see PHENIXDetector.h for overview and below for more specific details
//  
// Axel Drees 
// most recent update 11/16/2021
//
///////////////////////////////////////////////////////////////////////////////////

#include "PHENIXDetector.h"

////////////////////////////////////////////////////////////////////////////////////
//
// Constructor of PEHNIXDetector, calculates a number of internal members variables
//  
  PHENIXDetector::PHENIXDetector(){                // constructor

    for (int i=0; i<8; i++){                       // calculate sector edges in x-y plane for internal use
       SectorX0[i] = sin(P_centerPhi[i]*pi)*P_R_EMCal;
       SectorY0[i] = cos(P_centerPhi[i]*pi)*P_R_EMCal;
       P_R_EMCal_max[i] = P_R_EMCal*sqrt(1+pow(tan(P_dPhi[i]*pi),2));
       SectorXmin[i] = sin((P_centerPhi[i]-P_dPhi[i])*pi)*P_R_EMCal_max[i];
       SectorXmax[i] = sin((P_centerPhi[i]+P_dPhi[i])*pi)*P_R_EMCal_max[i];
       SectorYmin[i] = cos((P_centerPhi[i]-P_dPhi[i])*pi)*P_R_EMCal_max[i];
       SectorYmax[i] = cos((P_centerPhi[i]+P_dPhi[i])*pi)*P_R_EMCal_max[i];
    }
  }

  PHENIXDetector::~PHENIXDetector(void){}

/////////////////////////////////////////////////////////////////////////////////////////////
//
// PHENIX Butsyk acceptance for central arms in ++ field configuration
//
// input   TLorentzVector Particle     4 vector of input particle
//         Int_t   q                   charge
//
// returns 0 not in acceptance  
//         1 in arm 0
//         2 in arm 1
//
// for q = +1 or -1 particle will be in DC and RICH acceptance considering bent in B-field
// for q = 0 only ECAL acceptance is considered, no magnetic field
//
// assumes acceptance for ECAL and RICH are the same
//
// note orientation of phi using Lorentz Vectors is such that phi = 0 
// is along the PHENIX y-axis, i.e. straight up.  
//
//  
// Axel Drees  10/20/2018 - updated 9/30/2019
//
Int_t PHENIXDetector::InAcceptance(TLorentzVector particle, Int_t q)
{
  Double_t phi_rich, phi_DC, phi_ECAL;
  bool dep = 0;

  Double_t phi = particle.Phi();
  Double_t pT  = particle.Pt();
  Double_t eta   = particle.Eta();    

  if (eta > P_eta or eta < -P_eta)  return 0;
  if (pT<0.1)     return 0;

  phi_DC   = DCPhi(particle,q); 
  phi_rich = RICHPhi(particle,q); 
  phi_ECAL = EMCalPhi(particle,q);  

  arm = 0;
  if ( abs(q) == 1 &&
      ( phi_DC >= P_phi_east_min && phi_DC <= P_phi_east_max) && 
      ( phi_rich >= P_phi_east_min && phi_rich <= P_phi_east_max)) arm = arm_east;

  if ( abs (q) == 1 &&
      ( phi_DC >= P_phi_west_min && phi_DC <= P_phi_west_max) &&
      ( phi_rich >= P_phi_west_min && phi_rich <= P_phi_west_max))  arm = arm_west;
    
  if ( q == 0 && ( phi_ECAL >= P_phi_east_min && phi_ECAL <= P_phi_east_max)) arm =arm_east;
  if ( q == 0 && ( phi_ECAL >= P_phi_west_min && phi_ECAL <= P_phi_west_max)) arm = arm_west;

  return arm;
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// PHENIX EMCal Sector
//
// input   Double_t             - phi angle at EMCal
//
// returns sector number
//
//               E3   /5   4\   W3
//               E2  |6     3|  W2
//               E1  |7     2|  W1
//               E0   \8   1/   W0
//
// based on numbers from mEmcGeometryModule.C 
// 
// phi = 0 needs to be understood as pointing downward
//
// Axel Drees 10/25/2018
// 
Int_t PHENIXDetector::EMCalSector(Double_t phi){
   
  sector = 0;
  for(int i=0;i<8;i++){
//    cout << (P_centerPhi[i]+fdPhi[i])*pi << " " <<  phi << "  " << (P_centerPhi[i]-fdPhi[i])*pi << endl;
    if ( phi >= (P_centerPhi[i]-P_fdPhi[i])*pi && 
         phi <= (P_centerPhi[i]+P_fdPhi[i])*pi) {sector=i+1;} 
  }
  return sector; 
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// PHENIX EMCal Sector Coordinates
//
// input   Double_t             - phi angle at EMCal
//         Double_t             - theta (eqivalent to eta angle at mid rapidity)          
//
// return  Double_t             - y is coordinate in phi direction in cm from center
//         Double_t             - z is coordinate in eta direction in cm from center
//
//
// based on numbers from mEmcGeometryModule.C 
//
// Axel Drees 8/20/21
// 
Int_t PHENIXDetector::EMCalSectorCoordinates(Double_t phi, Double_t theta, Double_t& y, Double_t& z){

  sector=0;
  Double_t L = P_R_EMCal; // distance from origin in cm

  for(int i=0;i<8;i++){
    if ( phi >= (P_centerPhi[i]-P_dPhi[i])*pi && 
         phi <= (P_centerPhi[i]+P_dPhi[i])*pi) {
        sector=i+1;
        z = tan(pi/2.-theta)*L;
        y = tan(phi-P_centerPhi[i]*pi)*L;
//        std::cout << " -----------  y=" << y << "    z=" << z << " phi=" << tan(phi-P_centerPhi[i]*pi) << std::endl;
    } 
  }
  return sector; 
}
/////////////////////////////////////////////////////////////////////////////////
//
// PHENIX EMCAL Energy resolution 
//
// input   Double_t E           - input energy
//         Int_t    opt         - if 0 use PbSc resolution (default, 1 use PlGl, else external parameters
//
// returns smeared energy  
//
// uses values from ppg176
//          PbSc 8.1%/sqrt(E) + 2.1%
//          PbGl 5.9%/sqrt(E) + 0.8%
//
// lower energy cut at 100 MeV to avoid negative energies
// impliment poisson distribution to avoid negative energies
//
//          sigmaE/E = c1/sqrt(E) = 1/sqrt(Neff)         
//
// here Neff is the effective number of counts that th epoisson statistics is based on. 
//
// Axel Drees        10/23/2018  
// Axel Drees updated 7/26/2021
// dito               8/03/2021
// 
Double_t PHENIXDetector::SmearEnergy(Double_t energy, Int_t opt, Double_t c1, Double_t c2 ){

  Double_t E,E1,E2,Neff,sigmaE;

  if (energy < 0.001) return 0.;  

  if (opt == 0){                                         // set resolution constants 
    c1 = sigma_E_PbSc_c1;
    c2 = sigma_E_PbSc_c2; 
  }
  else if (opt == 1) {
    c1 = sigma_E_PbGl_c1;
    c2 = sigma_E_PbGl_c2;
  }

  Neff   = energy/c1/c1;                                 // calculate effective number of counts for 1/sqrt(e) term
  E1     = energy-energy/Neff * randy.PoissonD(Neff);  // get random variable and calculate energy shift
  sigmaE = energy*c2;                                    // calculate constant term   
  E2     = energy - randy.Gaus(energy,sigmaE);         // get random variable and calculate energy shift
  E      = energy + E1 + E2;                             // add energy shifts to energy

  if (E<0) E=0.;

//  std::cout << " PHENIX EMCal resolution " << energy << "   " << E << "   " << "  " << E1 << "  " << E2 << "  " << sigmaE/energy << "   " << Neff << endl;
  return E;
}
/////////////////////////////////////////////////////////////////////////////////
//
// PHENIX EMCAL position resolution 
//
// input   Double_t E           - input energy
//         Double_t x           - coordinate in y or z from sector center in cm
//         Int_t    opt         - if 0 use PbSc resolution (default, 1 use PlGl, else external parameters
//
// returns random smearing angle in radian to be added to phi/theta   
//
// uses values from Wenging Fan's thesis (page 55)
//          PbSc 5.9mm/sqrt(E) + 1.6mm 
//          PbGl 5.9%/sqrt(E) + 0.2mm
//
// extended to acount for impact angle term in PbSc
//          PbSc  0.8*21mm sin(theta)
//
// devide by 5m to convert coefficients to mrad
//
// Axel Drees        8/13/2021  
// 
Double_t PHENIXDetector::SmearEMCalPosition(Double_t energy, Double_t x, Int_t opt ){

  Double_t X = 0,sigmaX = 0;
  Double_t c1,c2,d;
  Double_t sintheta;               

  if (energy < 0.001) return 0.;  

  if (opt == 0){                                         // set resolution constants 
    c1 = sigma_x_PbSc_c1/P_R_EMCal;
    c2 = sigma_x_PbSc_c2/P_R_EMCal; 
    d  = sigma_x_PbSc_d/P_R_EMCal; 
  }
  else if (opt == 1) {
    c1 = sigma_x_PbGl_c1/P_R_EMCal;
    c2 = sigma_x_PbGl_c2/P_R_EMCal; 
    d  = sigma_x_PbGl_d/P_R_EMCal; 
  } else {
    return X;
  }

  sintheta = sin(abs(atan(x/P_R_EMCal)));
  sigmaX = sqrt(c1*c1/energy + c2*c2 + d*d*sintheta*sintheta);  // calculate energy dependent sigma 
  X = randy.Gaus(0,sigmaX);                             // get random variable and calculate energy shift

//  std::cout << " PHENIX EMCal resolution " << sintheta << "   " 
//  << x/P_R_EMCal << "   " << atan(x/P_R_EMCal) << std::endl;
  return X;
}
////////////////////////////////////////////////////////////////////////////////
//
//  PHENIX EMCal live area
//
//  input  Int_t                - sector number as deterimed by EMCalSector(phi)
//
//  returns .true. if in live area
//
// Live area is determined statistically based on dead maps from AuAU Run14
//
// Axel Drees  25/10/2018
//
Bool_t PHENIXDetector::EMCalLive(Int_t sector){
  Bool_t live=0;

  live =  (randy.Uniform(0.,1.) >= P_SectorLive[sector-1]);
  return live;
}
///////////////////////////////////////////////////////////////////////
//
//  PHENIX EMCal non linearity parameterization
//
//         E* = E (1 - c1*(1- exp(-c2/E)))
//
//  input  Double_t E     true energy
//         Double_t c1    constant 1       
//         Double_t c2    constant 2
//
//  returns nuon liner energy response
//
//  Axel Drees 10/24/2021
//
Double_t PHENIXDetector::NonLinearEnergy(Double_t energy, Double_t c0, Double_t c1, Double_t c2){

  Double_t E = c0*energy / (1 + c1*exp(-c2*energy));
  return E; 
}
///////////////////////////////////////////////////////////////////////
//
//  PHENIX VTX conversion probability 
//         
//
//  input  TRandom3       random generator
//
//         VTX layer 1&2  X/X0 ~1.3% each        
//         VTX layer 1&2  X/X0 ~5.2% each        
//
// use approximation 7/9 X/X0 as conversion probability 
//
//  returns layer of conversion 
//
//  Axel Drees 11/1/2018
//
Int_t PHENIXDetector::VTXConversion(){

  if (randy.Uniform(0.,1.) < 7./9.*VTX_X0[0]) return 1;
  if (randy.Uniform(0.,1.) < 7./9.*VTX_X0[1]) return 2;
  if (randy.Uniform(0.,1.) < 7./9.*VTX_X0[2]) return 3;
  if (randy.Uniform(0.,1.) < 7./9.*VTX_X0[3]) return 4;

  return 0; 
}

///////////////////////////////////////////////////////////////////////
//
// PHENIX DC momentum resolution
//
// Input: pt   in GeV 
//        q    charge
// Return:     smeared pt
//
// calculates angular resolution on alpha according to 
//   
//   sigma_alpha/alpha = 1/delta/kDC (sigma_ms (+) sigma_det*pt )  
//
// valid for charged particles from interaction point
//
// Axel Drees 9/15/2020
//
 Double_t PHENIXDetector::SmearPt(Double_t pt, Int_t &q){
   
    Double_t a;   
    Double_t alpha = alpha_DC; 
    Double_t sigma_ms = sigma_alpha_ms;    // units rad GeV
    Double_t sigma_det = sigma_alpha_det;   // units rad 

    Double_t sigma_alpha = sqrt(pow(sigma_ms/pt,2)+sigma_det*sigma_det);

//    std::cout << " pt " << pt << " alpha " << alpha << " sigma_alpha " << sigma_alpha;

    a = randy.Gaus(alpha,sigma_alpha);

//    std::cout << " reco alpha " << a << std::endl; 

    return P_delta*P_kDC/a;
 }
   
///////////////////////////////////////////////////////////////////////
//
// PENIX DC/PC1 angular resolution in phi
//
// input:  pt     in GeV
//         phi0   phi angle of particle at interaction point
// 
// return: smeared phi0
//
// calculates angular resolution from error propagation on:
//
//   phi0 = (1-delta)/delta alpha + phiDC  
//
// Axel Drees 9/16/2021
//
Double_t PHENIXDetector::SmearDCphi0(Double_t pt, Double_t phi0)
{
  Double_t phi;
  Double_t ff=(1-P_delta)/P_delta;

  Double_t sigma_phi0 = sqrt ( (pow(ff*sigma_alpha_ms,2)+pow(sigma_phi_ms,2))/pt/pt 
                              + pow(ff*sigma_alpha_det,2) );
  phi = randy.Gaus(phi0,sigma_phi0);

  // if (abs(phi-phi0)>0.1) {
  //   std::cout << " phi0 " << phi0 << " phi " << phi << " pt " << pt << " sigma " << sigma_phi0 << std::endl;
  // }
  return phi;
  
}

///////////////////////////////////////////////////////////////////////
//
// PENIX DC/PC1 angular resolution in eta
//
// input:  pt     in GeV
//         eta   pseudo rapidity of particle at interaction point
// 
// return: smeared eta
//
// calculates angular resolution from error propagation on:
//
//   phi0 = (1-delta)/delta alpha + phiDC  
//
// Axel Drees 9/16/2021
//
Double_t PHENIXDetector::SmearDCeta(Double_t pt, Double_t eta0)
{
  Double_t eta;

  Double_t sigma_eta0 = sqrt ( (pow(sigma_phi_ms,2))/pt/pt 
                              + pow(sigma_eta_det,2) );
  eta = randy.Gaus(eta0,sigma_eta0);

//  std::cout << " eta " << eta0 << " eta " << eta << " pt " << pt << " sigma " << sigma_eta0 << std::endl;
  return eta;
  
}

////////////////////////////////////////////////////////////////////////////////////////
//
// characterize track first to avoid multiple calculations of the same quantities
//
// Axel Drees 10/2/2021
//
void PHENIXDetector::CharacterizeTrack(TLorentzVector VT, Int_t q, Int_t id){
  
  Double_t y=0,z=0;
   
   alpha_EMCal = q*(P_kBfield-P_kEMCal_phi)/VT.Pt(); 
   alpha_DC    = P_delta*P_kDC / VT.Pt();
   arm         = InAcceptance(VT,q);
   phi_DC      = DCPhi(VT,q);
   phi_RICH    = RICHPhi(VT,q);
   phi_EMCal   = EMCalPhi(VT,q);
   sector      = EMCalSectorCoordinates(phi_EMCal, VT.Theta(), y, z);
   sectorY     = y;
   sectorZ     = z;
   sectorSinT  =  EMCalImpactAngle(VT, id);

//   std::cout << std::endl;
//   std::cout << "characterize track " << sector << " " << y << " " << z << std::endl;
 
}  

////////////////////////////////////////////////////////////////////////////////////////
//
// simulate track reconstruction
// 
// Axel Drees  9/17/2021
//

TLorentzVector PHENIXDetector::ReconstructTrack(TLorentzVector VT, Int_t q)
{
  TLorentzVector Reco;
  Double_t pt,eta,phi,m; 

  if (q != 0) {
    if (arm > 0){
      m    = VT.M();
      pt   = SmearPt(VT.Pt(), q);
      phi  = SmearDCphi0(VT.Pt(), VT.Phi());
      eta  = SmearDCeta(VT.Pt(), VT.Eta());
      Reco.SetPtEtaPhiM(pt,eta,phi,m);
//      std::cout << " --- track reco ----- " << VT.Pt() << "  " << q << "  " << pt << std::endl;
      return Reco;
    } else {
//      std::cout << "     --------- not in acceptance " << std::endl;      
      return Reco;
    }
  } else {
    std::cout << " not charged particle " << std::endl;
    return Reco;
  } 
}

Double_t PHENIXDetector::ElectronMomentum(Double_t E, Double_t p)
{
  Double_t value;
  float we,wp;
  if (sector>0 && sector<7){
    we = 1./(sigma_E_PbSc_c1*sigma_E_PbSc_c1/E+sigma_E_PbSc_c2*sigma_E_PbSc_c2);
  } else {
    we = 1./(sigma_E_PbSc_c1*sigma_E_PbSc_c1/E+sigma_E_PbSc_c2*sigma_E_PbSc_c2);
  }
  wp = 1/(sigma_p_c1*sigma_p_c1*p*p + sigma_p_c2*sigma_p_c2);

  value = (we*E+wp*p)/(we+wp);
  return value;

}
////////////////////////////////////////////////////////////////////////////////////////
//
// simulate shower reconstruction
// 
// currently only setup for electromagnetic showers for photons and electrons
// position resolution only implemented for photons
//
// Axel Drees  9/17/2021
//

TLorentzVector PHENIXDetector::ReconstructShower(TLorentzVector VT, Int_t id)
{
  TLorentzVector Reco;
  Int_t is,ia,iopt,q;
  Double_t a1,a2;
  Double_t E=0,px=0,py=0,pz=0; 
  Double_t y=0,z=0,theta=0;
  Double_t pt=0, phi=0, eta=0, m=0;

//  std::cout << abs(id) << "  " << pipID << std::endl;

  if ( id == photonID or  abs(id) == electronID) {
    q = -id/abs(id);
    if (id == photonID) q=0;
    is =  sector; 
    ia =  arm; 
    y = sectorY;
    z = sectorZ;

    if ( ia > 0 && is > 0){
      iopt = 0;
      if (is == 7 or is ==8) iopt = 1;
//    std::cout << " reconstructing photon in arm " << ia <<  " sector " << is << " and y = " << y << " and z =  " << z <<  std::endl;

      E = SmearEnergy(VT.E(),iopt);                 // apply default energy resolution
      if (abs(id) == photonID ) {  
        a1 = SmearEMCalPosition(VT.E(),y,iopt);  // apply position resolution to phi
        a2 = SmearEMCalPosition(VT.E(),z,iopt);  // apply position resolution to eta ~ theta
      } else {
        a1 = 0; 
        a2 = 0;
      }
      phi = VT.Phi() + a1;                         
      eta = VT.Eta() + a2;
      theta = pi/2. - 2.* atan(exp(-eta));
 //     theta = eta;
      pt = E*cos(theta);                         // recalculate pt
      m = VT.M();

//    std::cout << " eta " << eta << "  theta " << theta << std::endl;  
      Reco.SetPtEtaPhiM(pt,eta,phi,m);
//      Reco.Print();

//    std::cout << std::endl;



      return Reco;
    } else {
//      std::cout << " not in acceptance " << std::endl;      
      return Reco;
    }
  } else {
    std::cout << " not photon or electron " << std::endl;
    return Reco;
  } 
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// These routines calculate the energy correction depending on impact angle 
//
// addapted from EmcScSectorRec::CorrectECore(float Ecore, float x, float y, float* Ecorr)
// PHENIX CVS offline -> packages -> emc (see also Tadaaki Isobe's thesis equation 5.21)
//
// Axel Drees 9/26/2021
//
Double_t  PHENIXDetector::ShowerEnergyCorrection(Double_t energy, Double_t sinT)
// use sin(theta) calculated externally 
// must be used for charged particles
{
  const float par1 = 1; // 0.918;
  const float par2 = 1.35;
  const float par3 = 0.003;
  float corr;

   if( energy < 0.01 ) return energy;
  corr = par1 * ( 1 - par2*sinT*sinT*sinT*sinT*(1 - par3*log(energy)) );
//   std::cout << " correction with sinT " << sinT << "  corr=" << corr <<  std::endl;
  return corr; 
}

Double_t  PHENIXDetector::ShowerEnergyCorrection(Double_t energy, Double_t y, Double_t z)
// calculates sin(theta) from y,z position in local coordinates
// this corresponds to what is assumed in the EMCal reconstruction 
// but is valid only for photons
{
  const float par1 = 1.; // 0.918;
  const float par2 = 1.35;
  const float par3 = 0.003;
  float corr;
  float sinT; 

  sinT = sqrt((y*y+z*z)/(P_R_EMCal*P_R_EMCal+y*y+z*z));
  if( energy < 0.01 ) return energy;
  corr = par1 * ( 1 - par2*sinT*sinT*sinT*sinT*(1 - par3*log(energy)) );
//   std::cout << " correction with y,z " << sinT << "  corr=" << corr << std::endl;
  return corr; 
}

/////////////////////////////////////////////////////////////////////////////////
//
// calculates impact angle on calorimeter for electrons and photons 
// assumes particle is in PHENIX acceptance 
//
// Axel Drees 9/27/2021
//
Double_t PHENIXDetector::EMCalImpactAngle(TLorentzVector VT, Int_t id){

  Int_t is;
  Double_t y=0,z=0,pT;
  Double_t q, phi, tanphi, tantheta;
  Double_t sinT=-999;

  // std::cout << " --- id=" << id << "  phi=" << VT.Phi() << "  eta=" << VT.Eta() << "   pt=" << VT.Pt() << std::endl;

  y = sectorY;
  z = sectorZ;
  if (id == photonID)
  {
    if (sector>0) sinT = sqrt((y*y+z*z)/(P_R_EMCal*P_R_EMCal+y*y+z*z));
  } 
  else if (abs(id) == electronID) 
  {
//    q = -id/abs(id);
//    pT=VT.Pt();
//    is =  EMCalSectorCoordinates(phi_EMCal,VT.Theta(),y,z);
    if (sector>0) {
      tantheta = z/P_R_EMCal;
      phi  = atan(y/P_R_EMCal)+alpha_EMCal;
      tanphi = tan(phi);
      sinT = sin(atan(sqrt(tantheta*tantheta+tanphi*tanphi)));
    //  std::cout << " -+-+-+-+-+- sinT= " << sinT << "  or y z = " 
   //           << sqrt((y*y+z*z)/(P_R_EMCal*P_R_EMCal+y*y+z*z)) 
    //          << " alpha_EM= " << alpha_EMCal 
   //            << "  z= " << z/P_R_EMCal << "  tan(theta)=" << tan(theta) 
   //            << " phi " << y/P_R_EMCal << "  phi_EMCal = " << phi_EMCal 
   //           <<  << " y'= " << phi <<std::endl;
    //          << std::endl; 
    }
  }
  return sinT;

}

/////////////////////////////////////////////////////////////////////////////////
//
// calculates phi angle of track at DC,Rich, and calorimeter for electrons and photons 
// assumes particle is in PHENIX acceptance 
//
// Axel Drees 10/2/2021
//
Double_t PHENIXDetector::EMCalPhi(TLorentzVector VT, Double_t q){

  Double_t phi;
  Double_t phiDC,xDC,yDC,xE,yE;
  Double_t alpha_EMCal = q*(P_kBfield-P_kEMCal_phi)/VT.Pt();
  Double_t x[4],y[4], Px, Py; 
  Bool_t test = false , inter = false;

  if (q==0) {
    phi=VT.Phi();
  } else {
    phi=0;
    phiDC = DCPhi(VT,q); 
    x[0]  = sin(phiDC)*P_R_DC;
    y[0]  = cos(phiDC)*P_R_DC;
    x[1]  = x[0]+sin(phiDC-alpha_EMCal)*1.5*P_R_DC;
    y[1]  = y[0]+cos(phiDC-alpha_EMCal)*1.5*P_R_DC;

    while (!test) {
   
      for (int i=0; i<8; i++){
        x[2] = SectorXmin[i];
        y[2] = SectorYmin[i];
        x[3] = SectorXmax[i];
        y[3] = SectorYmax[i];
        inter = Intersect(Px,Py,x,y);
        if (inter) break;
      }
      test = true;
    } 
    phi = atan2(Px,Py);
  }

//  std::cout << std::endl;


  return phi;

}
Double_t PHENIXDetector::DCPhi(TLorentzVector VT, Double_t q){

  Double_t phi; 
  if(q==0) {
    phi = VT.Phi();
  } else {
    phi = VT.Phi() + q* P_kDC_phi/VT.Pt();
    if (phi>pi)    phi =  - (phi - 2.*pi);
    if (phi<-pi)   phi = - (phi + 2.*pi);
  }
  return phi;

}
Double_t PHENIXDetector::RICHPhi(TLorentzVector VT, Double_t q){

  Double_t phi;
  if(q==0) {
    phi = VT.Phi();
  } else {
    phi = VT.Phi() + q* P_kRICH_phi/VT.Pt();
    if (phi>pi)  phi = - (phi - 2.*pi);
    if (phi<-pi) phi = - (phi + 2.*pi);
  }
  return phi;

}

//////////////////////////////////////////////////////////////////////////////////////
//
// intersection of two line segments (x[0],y[0])(x[1],y[1]) and (x[2],y[2])(x[3],y[3])
// 
// if lines intersect it will return the coordinates Px,Py of intersection

Bool_t PHENIXDetector::Intersect(Double_t &Px,Double_t &Py,Double_t *x, Double_t *y){
  Bool_t inter = false;

  Double_t A1,B1,C1;
  Double_t A2,B2,C2;
  Px = 0;
  Py = 0;

  A1 = y[1]-y[0];
  B1 = x[0]-x[1];
  C1 = A1*x[0]+B1*y[0];

  A2 = y[3]-y[2];
  B2 = x[2]-x[3];
  C2 = A2*x[2]+B2*y[2];

  Double_t d = A1*B2-A2*B1; 

  if (d==0) {
    return inter;
  } else {
    Px = (B2*C1 - B1*C2)/d;
    Py = (A1*C2 - A2*C1)/d;
  }
  inter =    (std::min(x[0],x[1]) <= Px && Px <= std::max(x[0],x[1])) 
          && (std::min(x[2],x[3]) <= Px && Px <= std::max(x[2],x[3])) 
          && (std::min(y[0],y[1]) <= Py && Py <= std::max(y[0],y[1])) 
          && (std::min(y[2],y[3]) <= Py && Py <= std::max(y[2],y[3]));
  if (!inter) {
    Px=0;
    Py=0;
  }
  return inter;
}