//// KinematicDistributions.h //////////////////////////////////////////////////////////////////////////
//
// collection of momentum spectra of various particles
//
// currently implemented:
//
//   TF1 Hagedorn (name, mass, upperlim, lowerlim, A=504.5, a=0.5169, b=0.1626, p0=0.7366, n=-8.274)
//   TF1 HagedornYield (name, mass, upperlim, lowerlim, A=504.5, a=0.5169, b=0.1626, p0=0.7366, n=-8.274)
//       default parameters for Au+Au 200 (see Alan Dion thesis for details)
//       pp 200 parameters: 377., 0.356, 0.068, 0.7, -8.25
//
//   TF1 PromptPhotonYield(name, upperlim, lowerlim, a=0.0066, b=6.4, c=0.4, n=17.6, s=200.)
//       Prompt photon function fit from PPG140 fit to pp 200
//
//   double  Weight_GPR_pion_pp200(pt, opt=0)
//           - uses GPR for pi0 from ppg202, and pi+/pi-  ppg030 & ppg101  - Roli Esha 11/10/2021  
//           - return value+opt*error
//   double  Weight_GPR_etapi_uni(pt, opt=0)
//           - universal eta/pi ratio fro Yuanjie Ren thesis  
//
// Axel Drees  11/11/2021 
//

#ifndef KinDist_h
#define KinDist_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>

#include "/GPRfiles/GPRFileNames.h"

//*******************************************************************
//
// Hagedorn Function for invariant yield mesons, 1/pt dN/dpt
//
// see thesis of Alan Dion for papameters for pi0
// all other particle distributions are derived using mt scaling
// assumes mt scaling 
//
Double_t Hagedorn_Func(Double_t *x, Double_t *p) {
//
// Hagedorn function in 1/pt dN/pt from ppg088 
//
  Double_t pt   = x[0];
  Double_t mass = p[0];
  Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.13498*0.13498);
  Double_t A    = p[1];
  Double_t a    = p[2];
  Double_t b    = p[3];
  Double_t p0   = p[4];
  Double_t n    = p[5];
 
  Double_t value = A* pow( exp(-a*mt-b*mt*mt)+mt/p0 ,n);
  return value;
}
//*******************************************************************
//
// Hagedorn Function for mesons yield, dN/dpt
//
// see thesis of Alan Dion for papameters for pi0
// all other particle distributions are derived using mt scaling
// assumes mt scaling 
//
Double_t Hagedorn_Yield_Func(Double_t *x, Double_t *p) {
//
// Hagedorn function in 1/pt dN/pt from ppg088 
//
  Double_t pt   = x[0];
  Double_t mass = p[0];
  Double_t mt   = TMath::Sqrt(pt * pt + mass * mass - 0.13498*0.13498);
  Double_t A    = p[1];
  Double_t a    = p[2];
  Double_t b    = p[3];
  Double_t p0   = p[4];
  Double_t n    = p[5];
 
  Double_t value = pt*A/2/3.141* pow( exp(-a*mt-b*mt*mt)+mt/p0 ,n);
  return value;
}

///////////////////////////////////////////////////////////////////////////////
TF1 *HagedornYield(const Char_t *name, Double_t mass, Double_t upperlim, Double_t lowerlim, Double_t A=504.5, Double_t a=0.5169, Double_t b=0.1626, Double_t p0=0.7366, Double_t n=-8.274){
//
// modified Hagedorn function used in ppg088, default parameters for min bias AuAu as documented in Alan Dion's PhD thesis
//
//
  if (lowerlim <= 0.) lowerlim = 0.001;
  TF1 *fHagedorn = new TF1(name, Hagedorn_Yield_Func, lowerlim, upperlim, 6);
  fHagedorn->FixParameter(0,mass); 
  fHagedorn->SetParameters(mass, A, a, b, p0, n); 
  fHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return fHagedorn;
    }

///////////////////////////////////////////////////////////////////////////////
TF1 *Hagedorn(const Char_t *name, Double_t mass, Double_t upperlim, Double_t lowerlim, Double_t A=504.5, Double_t a=0.5169, Double_t b=0.1626, Double_t p0=0.7366, Double_t n=-8.274){
//
// modified Hagedorn function used in ppg088, default parameters for min bias AuAu as documented in Alan Dion's PhD thesis
//
//
  if (lowerlim <= 0.) lowerlim = 0.001;
  TF1 *fHagedorn = new TF1(name, Hagedorn_Func, lowerlim, upperlim, 6);
  fHagedorn->FixParameter(0,mass); 
  fHagedorn->SetParameters(mass, A, a, b, p0, n); 
  fHagedorn->SetParNames("mass", "A", "a", "b", "p0", "n");
  return fHagedorn;
    }

//*******************************************************************
//
Double_t PromptPhoton_Func(Double_t *x, Double_t *p) {
//
// Prompt photon function fit from PPG140 
//
// the units are in mb/GeV^2
//
  Double_t pt    = x[0];
  Double_t a     = p[0];
  Double_t b     = p[1];
  Double_t c     = p[2];
  Double_t n     = p[3];
  Double_t s     = p[4];
  Double_t xt    = 2.*pt/s;
  Double_t value = pt* a / pow(pt, b+c*log(xt)) * pow(1-xt*xt,n);
//  std::cout << "... fit ... pt " << pt << " xt " << xt << " a " << a << " value " << value << std::endl;
  return value;
}

TF1 *PromptPhotonYield(const Char_t *name, Double_t upperlim, Double_t lowerlim, Double_t a=0.0066, Double_t b=6.4, Double_t c=0.4, Double_t n=17.6, Double_t s=200.){
//
//
  if (lowerlim <= 0.) lowerlim = 0.001;
  TF1 *fHagedorn = new TF1(name, PromptPhoton_Func, lowerlim, upperlim, 5);
  fHagedorn->FixParameter(5,s); 
  fHagedorn->SetParameters(a, b, c, n, s); 
  fHagedorn->SetParNames("a", "b", "c", "n", "s");
  return fHagedorn;
    }

////////////////////////////////////////////////////////////////////////////////////////////
//
// distributions created with Gaussian Process Regression
//    
// txt files located in /GPRfiles/... 
//
//    GPR_pion_pp200.txt - uses pi0 from ppg202, and pi+/pi-  ppg030 & ppg101  - Roli Esha 11/10/2021
//    GPR_etapi_universal.txt - for pp  and pA collisons at all energies - from Yuanjie Ren MA thesis 2020
//
const int GPR_pi_pp200_n=1000;
const double GPR_pi_pp200_ptmax = 25.;   // GeV/c
const double GPR_pi_pp200_ptmin = 0.3;   // GeV/c
const double GPR_pi_pp200_gap = (GPR_pi_pp200_ptmax-GPR_pi_pp200_ptmin)/GPR_pi_pp200_n;

const int GPR_etapi_uni_n=500;
const double GPR_etapi_uni_ptmax = 1.477;  // log10(pt in GeV/c)
const double GPR_etapi_uni_ptmin = -1.8;    // log10(pt in GeV/c)
const double GPR_etapi_uni_gap = (GPR_etapi_uni_ptmax-GPR_etapi_uni_ptmin)/GPR_etapi_uni_n;

double GPR_pi_pp200_pt[GPR_pi_pp200_n], GPR_pi_pp200_yield[GPR_pi_pp200_n], GPR_pi_pp200_yield_error[GPR_pi_pp200_n];
double GPR_etapi_uni_pt[GPR_etapi_uni_n], GPR_etapi_uni[GPR_etapi_uni_n], GPR_etapi_uni_error[GPR_etapi_uni_n];

//////// pp 200 pion pt distribution ////////////////////////////////////////////////////////////////
// read textfile
void Read_GPR_pion_pp200(){
  std::ifstream infile;
  infile.open(GPR_pion_pp200);
  for(int i = 0; i < GPR_pi_pp200_n; i++){
    infile >> GPR_pi_pp200_pt[i] >> GPR_pi_pp200_yield[i] >> GPR_pi_pp200_yield_error[i];
  }
  infile.close();
}
// calculate weight 
double Weight_GPR_pion_pp200(double pt, double opt=0){
  
  double gap, value=0, error=0, target_index, alpha;
  int index;

  if (pt<=GPR_pi_pp200_ptmin or pt>=GPR_pi_pp200_ptmax) return 0;
  gap = GPR_pi_pp200_gap;
  target_index = (pt-GPR_pi_pp200_ptmin)/gap;
  index = floor(target_index);
  alpha = target_index - index;

  if (index < GPR_pi_pp200_n-1) {
     value = GPR_pi_pp200_yield[index]*(1-alpha) + GPR_pi_pp200_yield[index+1]*alpha; 
     error = GPR_pi_pp200_yield_error[index]*(1-alpha) + GPR_pi_pp200_yield_error[index+1]*alpha;
  }
//   std::cout << target_index << "\t" << alpha << "\t" << value << "\t" << error << std::endl;
  return value+opt*error;
}
////////// eta to pi0 ration univesal function valid for all pp and pA collisions //////////////////////
void Read_GPR_etapi_uni(){
  std::string trashline;
  std::ifstream infile;
  infile.open(GPR_etapi_universal);
  getline(infile,trashline); 
//  std::cout<<"trash: "<<trashline<<std::endl;  
  getline(infile,trashline); 
//  std::cout<<"trash: "<<trashline<<std::endl;  
  for(int i = 0; i < GPR_etapi_uni_n; i++){
    infile >> GPR_etapi_uni_pt[i] >> GPR_etapi_uni[i] >> GPR_etapi_uni_error[i];
//    std::cout<<"i="<<i<<", (x,y,ey)=("<<GPR_etapi_uni_pt[i]<<","<<GPR_etapi_uni[i]<<","<<GPR_etapi_uni_error[i]<<")"<<std::endl;
  }
  infile.close();
}
double Weight_GPR_etapi_uni(double pt, double opt=0){
  double gap, lpt, value=0, error=0, target_index, alpha;
  int index;

  lpt = log10(pt);

  if (lpt>=GPR_etapi_uni_ptmax) return 0;
  if (lpt<=GPR_etapi_uni_ptmin) return GPR_etapi_uni[0]+opt*GPR_etapi_uni_error[0];
  gap = GPR_etapi_uni_gap; 
  target_index = (lpt-GPR_etapi_uni_ptmin)/gap;
  index = floor(target_index);
  alpha = target_index - index;

  if (index < GPR_etapi_uni_n-1) {
     value = GPR_etapi_uni[index]*(1-alpha) + GPR_etapi_uni[index+1]*alpha; 
     error = GPR_etapi_uni_error[index]*(1-alpha) + GPR_etapi_uni_error[index+1]*alpha;
  }  
  return value+opt*error;
}

double Pt_Mt_Scaled(double pt, double m, double m_ref){

  double pt_ref = sqrt(pt*pt + m*m - m_ref*m_ref);
  return pt_ref;
}

#endif
