#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
// 
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include "HELIOSLibrary/HELIOSLibrary.h"
#include "MyPlotting/MyPlot.C"


using namespace std;

///////////////////////////////////////////////////////////////////////////////
void TestDecay(){
//
  TRandom3 MyRandy = TRandom3(0);                // Random Generator
  TLorentzVector temp;
  Double_t ww = 0; 
  Int_t id = 0;

  Int_t iopt = 2;                                   

  cout << endl;

  MyPlot plot;

  Particle pi0("pi0");
  Particle eta("eta");
  Particle omega("omega");
  Particle rho0("rho0");
  Particle etap("etap");
  Particle phi_meson("phi");
//  pi0.Print();
  Particle gamma1("photon");
  Particle gamma2("photon");
  Particle electron("electron");
  Particle positron("positron");
  TLorentzVector eepair;


  Double_t E1, E2, OpeningAngle;
  Double_t e1, e2;
  Double_t pt_min = 0.1;
  Double_t pt_max = 25.;
  Double_t pt_low = 0.1;
  Double_t pt_high = 15.;
  Double_t ptcut = 0.2;
  Double_t ecut = 0.5;
  Int_t nevt = 1000000;
  Bool_t ek  = 1;
  Bool_t out = 0;
  Bool_t elecFound = false;
  Bool_t posiFound = false;
  Double_t weight = 1.;
  Int_t ndecay =0;

  TF1 *piHagedorn      = HagedornYield("piHagedorn", pi0Mass, pt_max, pt_min);
  TF1 *etaHagedorn     = HagedornYield("etaHagedorn", etaMass, pt_max, pt_min);
  TF1 *rho0Hagedorn    = HagedornYield("omegaHagedorn", rho0Mass, pt_max, pt_min);
  TF1 *omegaHagedorn   = HagedornYield("omegaHagedorn", omegaMass, pt_max, pt_min);
  TF1 *etapHagedorn    = HagedornYield("etapHagedorn", etapMass, pt_max, pt_min);
  TF1 *phiHagedorn     = HagedornYield("phiHagedorn", phiMass, pt_max, pt_min);
//  TF1 *pi0Phi         = new TF1("pi0Phi","1+[0]*x",-2*pi,2*pi);
//  TF1 *pi0Rapidity    = new TF1("pi0Rapidity","1+[0]*x",-.5,.5);

  TH1D *h_ptpi             = new TH1D("h_ptpi","pt",100,pt_low,pt_high);
  TH1D *h_pt_gpi           = new TH1D("h_pt_gpi","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_epi           = new TH1D("h_pt_epi","pt_e",100,pt_low,pt_high);
  TH1D *h_masspi           = new TH1D("hmasspi","mass",300,0.,1.5);

// oid plot.StyleMe(TH1D *tE, int marker = 20, int color = 1, double msize = 1.2, int lstyle = 1, float lwidth = 2)
  plot.StyleMe(h_ptpi,   20, kBlue, 1., 1, 2); 
  plot.StyleMe(h_pt_gpi, 20, kBlue, 1., 2, 2); 
  plot.StyleMe(h_pt_epi, 20, kBlue, 1., 3, 2); 
  plot.StyleMe(h_masspi, 20, kBlue, 1., 1, 2); 
  
  TH1D *h_pteta             = new TH1D("h_pteta","pt",100,pt_low,pt_high);
  TH1D *h_pt_geta           = new TH1D("h_pt_geta","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_eeta           = new TH1D("h_pt_eeta","pt_e",100,pt_low,pt_high);
  TH1D *h_masseta           = new TH1D("hmasseta","mass",300,0.,1.5);

  plot.StyleMe(h_pteta,   20, kRed, 1., 1, 2); 
  plot.StyleMe(h_pt_geta, 20, kRed, 1., 2, 2); 
  plot.StyleMe(h_pt_eeta, 20, kRed, 1., 3, 2); 
  plot.StyleMe(h_masseta, 20, kRed, 1., 1, 2); 

  TH1D *h_ptrho0             = new TH1D("h_ptrho0","pt",100,pt_low,pt_high);
  TH1D *h_pt_grho0           = new TH1D("h_pt_grho0","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_erho0           = new TH1D("h_pt_erho0","pt_e",100,pt_low,pt_high);
  TH1D *h_massrho0           = new TH1D("hmassrho0","mass",300,0.,1.5);

  plot.StyleMe(h_ptrho0,   20, kMagenta, 1., 1, 2); 
  plot.StyleMe(h_pt_grho0, 20, kMagenta, 1., 2, 2); 
  plot.StyleMe(h_pt_erho0, 20, kMagenta, 1., 3, 2); 
  plot.StyleMe(h_massrho0, 20, kMagenta, 1., 1, 2); 

  TH1D *h_ptomega             = new TH1D("h_ptomega","pt",100,pt_low,pt_high);
  TH1D *h_pt_gomega           = new TH1D("h_pt_gomega","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_eomega           = new TH1D("h_pt_eomega","pt_e",100,pt_low,pt_high);
  TH1D *h_massomega           = new TH1D("hmassomega","mass",300,0.,1.5);

  plot.StyleMe(h_ptomega,   20, kGreen+2, 1., 1, 2); 
  plot.StyleMe(h_pt_gomega, 20, kGreen+2, 1., 2, 2); 
  plot.StyleMe(h_pt_eomega, 20, kGreen+2, 1., 3, 2); 
  plot.StyleMe(h_massomega, 20, kGreen+2, 1., 1, 2); 

  TH1D *h_ptetap             = new TH1D("h_ptetap","pt",100,pt_low,pt_high);
  TH1D *h_pt_getap           = new TH1D("h_pt_getap","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_eetap           = new TH1D("h_pt_eetap","pt_e",100,pt_low,pt_high);
  TH1D *h_massetap           = new TH1D("hmassetap","mass",300,0.,1.5);

  plot.StyleMe(h_ptetap,   20, kOrange+2, 1., 1, 2); 
  plot.StyleMe(h_pt_getap, 20, kOrange+2, 1., 2, 2); 
  plot.StyleMe(h_pt_eetap, 20, kOrange+2, 1., 3, 2); 
  plot.StyleMe(h_massetap, 20, kOrange+2, 1., 1, 2); 

  TH1D *h_ptphi             = new TH1D("h_ptphi","pt",100,pt_low,pt_high);
  TH1D *h_pt_gphi           = new TH1D("h_pt_gphi","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_ephi           = new TH1D("h_pt_ephi","pt_e",100,pt_low,pt_high);
  TH1D *h_massphi           = new TH1D("hmassphi","mass",300,0.,1.5);

  plot.StyleMe(h_ptphi,   20, kTeal, 1., 1, 2); 
  plot.StyleMe(h_pt_gphi, 20, kTeal, 1., 2, 2); 
  plot.StyleMe(h_pt_ephi, 20, kTeal, 1., 3, 2); 
  plot.StyleMe(h_massphi, 20, kTeal, 1., 1, 2); 



  TH1D *h_pt_g           = new TH1D("h_pt_g","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_e           = new TH1D("h_pt_e","pt_e",100,pt_low,pt_high);
  TH1D *h_mass           = new TH1D("hmass","mass",300,0.,1.5);

  plot.StyleMe(h_pt_g,   20, kBlack, 1., 1, 2); 
  plot.StyleMe(h_pt_e,   20, kBlack, 1., 1, 2); 
  plot.StyleMe(h_mass,   20, kBlack, 1., 1, 2); 


  for (int i=0; i<nevt; i++){                                     // event loop
    if(i%100000==0) cout << "event loop at " << i << endl;         // event counter printed 
/////////////////////////////////////////////// pi decays ////////////////////////////////////////////////
    pi0.ResetP();
    pi0.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi0 with flat distribution
    pi0.SetWeight(piHagedorn->Eval(pi0.Pt())/float(nevt));


//    pi0.GenerateP(piHagedorn,pi0Phi,pi0Rapidity);               // generate 4 vector for pi0 
//    pi0.GenerateP(piHagedorn);                                  // generate 4 vector for pi0 
    h_ptpi->Fill(pi0.Pt(),pi0.Weight());                         // histogram pion proprety

    if (iopt == 2) pi0.DecayFlat();
    if (iopt == 1) pi0.Decay();
    if (iopt == 4) pi0.DecaySingleBranch("pi0->Dalitz");
    if (iopt == 3) pi0.DecaySingleBranch("pi0->gg");


    ndecay = pi0.GetNumberOfDaughters();
    if (i<10)  cout << "----- Pi0Deacy --------------------------------------------------------" << endl; 

    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 
    if (i<10) cout << " weight " << pi0.GetDaughterWeight(0) << endl;

    elecFound = false;
    posiFound = false;    

    for (Int_t j=0; j< ndecay; j++) {
      temp = pi0.GetDecayDaughter(j);
      if (i<10) temp.Print();
      ww   = pi0.GetDaughterWeight(j)*pi0.Weight();
      id   = pi0.GetDaughterID(j);
      if (i<10) cout << " weight " << ww << " id: " << id << endl;
      if (id == photonID) {
         h_pt_gpi->Fill(temp.Pt(),ww);    
      } else if (id == electronID) {
        electron.UpdateParticle(temp);
        electron.SetWeight(ww);
        elecFound = true;
        h_pt_epi->Fill(electron.Pt(),electron.Weight());          
//          if (i<10) electron.Print();
      } else if (id == positronID) {
        positron.UpdateParticle(temp);
        positron.SetWeight(ww);
        posiFound = true;
        h_pt_epi->Fill(positron.Pt(),positron.Weight());         
//          if (i<10) positron.Print();
      }
    }
    if (elecFound && posiFound) {
      eepair = electron+positron;
//      eepair.Print();
      if (i<10) eepair.Print();
      h_masspi->Fill(eepair.M(),ww);     
    }
/////////////////////////////////////////////// eta decays ////////////////////////////////////////////////
    eta.ResetP();
    eta.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi0 with flat distribution
    eta.SetWeight(etaHagedorn->Eval(eta.Pt())/float(nevt)*0.5);
    h_pteta->Fill(eta.Pt(),eta.Weight());                         // histogram pion proprety
      
    if (iopt == 1) eta.Decay();
    if (iopt == 2) eta.DecayFlat();
    if (iopt == 3) eta.DecaySingleBranch("eta->gg");
    if (iopt == 4) eta.DecaySingleBranch("eta->Dalitz");

    ndecay = eta.GetNumberOfDaughters();
    if (i<10)  cout << "----- EtaDeacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 
     
    elecFound = false;
    posiFound = false;    

    for (Int_t j=0; j< ndecay; j++) {
      temp = eta.GetDecayDaughter(j);
      if (i<10) temp.Print();
      ww   = eta.GetDaughterWeight(j)*eta.Weight();
      id   = eta.GetDaughterID(j);
      if (i<10) cout << " weight " << ww << " id: " << id << endl;
      if (id == photonID) {
         h_pt_geta->Fill(temp.Pt(),ww);    
      } else if (id == electronID) {
        electron.UpdateParticle(temp);
        electron.SetWeight(ww);
        elecFound = true;
        h_pt_eeta->Fill(electron.Pt(),electron.Weight());          
//          if (i<10) electron.Print();
      } else if (id == positronID) {
        positron.UpdateParticle(temp);
        positron.SetWeight(ww);
        posiFound = true;
        h_pt_eeta->Fill(positron.Pt(),positron.Weight());         
//          if (i<10) positron.Print();
      }
    }
    if (elecFound && posiFound) {
      eepair = electron+positron;
//        eepair.Print();
      if (i<10) eepair.Print();
      h_masseta->Fill(eepair.M(),ww);     
    }
/////////////////////////////////////////////// rho0 decays ////////////////////////////////////////////////
    rho0.ResetP();
    electron.ResetP();
    positron.ResetP();
    rho0.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi0 with flat distribution
    rho0.SetWeight(rho0Hagedorn->Eval(rho0.Pt())/float(nevt)*1.15);
    h_ptrho0->Fill(rho0.Pt(),rho0.Weight());                         // histogram pion proprety

    if (iopt == 1) rho0.Decay();
    if (iopt == 2) rho0.DecayFlat();
//      if (iopt == 5) rho0.DecaySingleBranch("rho0->ee");

    ndecay = rho0.GetNumberOfDaughters();
    if (i<10)  cout << "----- rho0Deacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 

    elecFound = false;
    posiFound = false;    
      
    for (Int_t j=0; j< ndecay; j++) {
       temp = rho0.GetDecayDaughter(j);
       if (i<10) temp.Print();
       ww   = rho0.GetDaughterWeight(j)*rho0.Weight();
       id   = rho0.GetDaughterID(j);
       if (i<10) cout << " weight " << ww << " id: " << id << endl;
       if (id == photonID) {
          h_pt_grho0->Fill(temp.Pt(),ww);    
       } else if (id == electronID) {
         electron.UpdateParticle(temp);
         electron.SetWeight(ww);
         elecFound = true;
//           h_pt_eetap->Fill(electron.Pt(),electron.Weight());          
           if (i<10) electron.Print();
       } else if (id == positronID) {
         positron.UpdateParticle(temp);
         positron.SetWeight(ww);
         posiFound = true;
//           h_pt_eetap->Fill(positron.Pt(),positron.Weight());         
           if (i<10) positron.Print();
       }
    }   
    if (elecFound && posiFound) {
       eepair = electron+positron;
       if (i<10) cout << "                            " << eepair.M() << endl;
       if (i<10) eepair.Print();
       h_massrho0->Fill(eepair.M(),ww);
    }    
/////////////////////////////////////////////// omega decays ////////////////////////////////////////////////
    omega.ResetP();
    electron.ResetP();
    positron.ResetP();
    omega.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi0 with flat distribution
    omega.SetWeight(omegaHagedorn->Eval(omega.Pt())/float(nevt)*0.9);
    h_ptomega->Fill(omega.Pt(),omega.Weight());                         // histogram pion proprety

    
    if (iopt == 1) omega.Decay();
    if (iopt == 2) omega.DecayFlat();
    if (iopt == 3) omega.DecaySingleBranch("omega->gg");
    if (iopt == 4) omega.DecaySingleBranch("omega->pi0ee");
//      if (iopt == 5) omega.DecaySingleBranch("omega->ee");

    ndecay = omega.GetNumberOfDaughters();
    if (i<10)  cout << "----- omegaDeacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 

    elecFound = false;
    posiFound = false;    
      
    for (Int_t j=0; j< ndecay; j++) {
       temp = omega.GetDecayDaughter(j);
       if (i<10) temp.Print();
       ww   = omega.GetDaughterWeight(j)*omega.Weight();
       id   = omega.GetDaughterID(j);
       if (i<10) cout << " weight " << ww << " id: " << id << endl;
       if (id == photonID) {
          h_pt_gomega->Fill(temp.Pt(),ww);    
       } else if (id == electronID) {
         electron.UpdateParticle(temp);
         electron.SetWeight(ww);
         elecFound = true;
//           h_pt_eetap->Fill(electron.Pt(),electron.Weight());          
           if (i<10) electron.Print();
       } else if (id == positronID) {
         positron.UpdateParticle(temp);
         positron.SetWeight(ww);
         posiFound = true;
//           h_pt_eetap->Fill(positron.Pt(),positron.Weight());         
           if (i<10) positron.Print();
       }
    }   
    if (elecFound && posiFound) {
       eepair = electron+positron;
       if (i<10) cout << "                            " << eepair.M() << endl;
       if (i<10) eepair.Print();
       h_massomega->Fill(eepair.M(),ww);
    }    
/////////////////////////////////////////////// eta' decays ////////////////////////////////////////////////
    etap.ResetP();
    electron.ResetP();
    positron.ResetP();
    etap.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi0 with flat distribution
    etap.SetWeight(etapHagedorn->Eval(etap.Pt())/float(nevt)*0.25);
    h_ptetap->Fill(etap.Pt(),etap.Weight());                         // histogram pion proprety
      
    if (iopt == 1) etap.Decay();
    if (iopt == 2) etap.DecayFlat();
    if (iopt == 3) etap.DecaySingleBranch("etap->gg");
    if (iopt == 4) etap.DecaySingleBranch("etap->Dalitz");

    ndecay = etap.GetNumberOfDaughters();
    if (i<10)  cout << "----- etapDeacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 
     
    elecFound = false;
    posiFound = false;    

    for (Int_t j=0; j< ndecay; j++) {
      temp = etap.GetDecayDaughter(j);
      if (i<10) temp.Print();
      ww   = etap.GetDaughterWeight(j)*etap.Weight();
      id   = etap.GetDaughterID(j);
      if (i<10) cout << " weight " << ww << " id: " << id << endl;
      if (id == photonID) {
         h_pt_getap->Fill(temp.Pt(),ww);    
      } else if (id == electronID) {
        electron.UpdateParticle(temp);
        electron.SetWeight(ww);
        elecFound = true;
        h_pt_eetap->Fill(electron.Pt(),electron.Weight());          
//          if (i<10) electron.Print();
      } else if (id == positronID) {
        positron.UpdateParticle(temp);
        positron.SetWeight(ww);
        posiFound = true;
        h_pt_eetap->Fill(positron.Pt(),positron.Weight());         
//          if (i<10) positron.Print();
      }
    }
    if (elecFound && posiFound) {
      eepair = electron+positron;
//        eepair.Print();
      if (i<10) eepair.Print();
      h_massetap->Fill(eepair.M(),ww);     
    }
/////////////////////////////////////////////// phi decays ////////////////////////////////////////////////
    phi_meson.ResetP();
    electron.ResetP();
    positron.ResetP();
    phi_meson.GenerateP(pt_min,pt_max);                                // generate 4 vector for pi0 with flat distribution
    phi_meson.SetWeight(phiHagedorn->Eval(phi_meson.Pt())/float(nevt)*0.25);
    h_ptphi->Fill(phi_meson.Pt(),phi_meson.Weight());                         // histogram pion proprety
    
    if (iopt == 1) phi_meson.Decay();
    if (iopt == 2) phi_meson.DecayFlat();
    if (iopt == 5) phi_meson.DecaySingleBranch("phi->ee");

    ndecay = phi_meson.GetNumberOfDaughters();
    if (i<10)  cout << "----- phiDeacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 
     
    elecFound = false;
    posiFound = false;    

   for (Int_t j=0; j< ndecay; j++) {
     temp = phi_meson.GetDecayDaughter(j);
     if (i<10) temp.Print();
     ww   = phi_meson.GetDaughterWeight(j)*phi_meson.Weight();
     id   = phi_meson.GetDaughterID(j);
     if (i<10) cout << " weight " << ww << " id: " << id << endl;
     if (id == photonID) {
        h_pt_gphi->Fill(temp.Pt(),ww);    
     } else if (id == electronID) {
       electron.UpdateParticle(temp);
       electron.SetWeight(ww);
       elecFound = true;
       h_pt_ephi->Fill(electron.Pt(),electron.Weight());          
//          if (i<10) electron.Print();
     } else if (id == positronID) {
       positron.UpdateParticle(temp);
       positron.SetWeight(ww);
       posiFound = true;
       h_pt_ephi->Fill(positron.Pt(),positron.Weight());         
//          if (i<10) positron.Print();
     }
   }
  if (elecFound && posiFound) {
    eepair = electron+positron;
//        eepair.Print();
    if (i<10) eepair.Print();
    h_massphi->Fill(eepair.M(),ww);     
  }
  

 ////////////////////////////////////////////////////////////////////////////
  }

  TH1D *h_etapi       = (TH1D*)h_pteta->Clone("h_etapi");
  TH1D *h_omegapi     = (TH1D*)h_ptomega->Clone("h_omegapi");
  TH1D *h_etappi      = (TH1D*)h_ptetap->Clone("h_etappi");
  
  TH1D *h_g_pih       = (TH1D*)h_pt_gpi->Clone("h_g_pih");
  TH1D *h_g_etah      = (TH1D*)h_pt_geta->Clone("h_g_etah");
  TH1D *h_g_omegah    = (TH1D*)h_pt_gomega->Clone("h_g_omegah");
  TH1D *h_g_etaph     = (TH1D*)h_pt_getap->Clone("h_g_etaph");
 
  h_pt_g->Add(h_pt_gpi);
  h_pt_g->Add(h_pt_geta);
  h_pt_g->Add(h_pt_gomega);
  h_pt_g->Add(h_pt_getap);

  h_g_pih->Divide(h_pt_g);
  h_g_etah->Divide(h_pt_g);
  h_g_omegah->Divide(h_pt_g);
  h_g_etaph->Divide(h_pt_g);

  h_etapi->Divide(h_ptpi);
  h_omegapi->Divide(h_ptpi);
  h_etappi->Divide(h_ptpi);
  
  h_mass->Add(h_masspi);
  h_mass->Add(h_masseta);
  h_mass->Add(h_massomega);
  h_mass->Add(h_massetap);
  h_mass->Add(h_massrho0);
  h_mass->Add(h_massphi);


  TCanvas *c1 = plot.Canvas ("c1",400,600,10,10,1);
  TH1D *pt1 = plot.Frame("pt1","pt","counts",pt_low,pt_high*.999,1e-12,1.);
  pt1->Draw();
   h_ptpi->Draw("sameL Chist");
   h_pteta->Draw("sameL Chist");
   h_ptomega->Draw("sameL Chist");
   h_ptrho0->Draw("sameL Chist");
   h_ptetap->Draw("sameL Chist");
  h_ptphi->Draw("sameL Chist");

  TCanvas *c11 = plot.Canvas ("c11",400,600,410,10,1);
  TH1D *pt2 = plot.Frame("pt2","pt","counts",pt_low,pt_high*.999,1e-12,1.);
  pt2->Draw();
  h_pt_gpi->Draw("same Chist" );
  h_pt_geta->Draw("same Chist");
  h_pt_gomega->Draw("same Chist");
  h_pt_getap->Draw("same Chist");

//  h_pt_eomega->Draw("same");

  TCanvas *c2 = plot.Canvas ("c2",400,300,10,610);
  TH1D *ptr = plot.Frame("ptr","pt","ratio",pt_low,pt_high*.999,0.0,0.99);
  ptr->Draw();
  h_etapi->Draw("same Chist" );
  h_omegapi->Draw("same Chist" );
  h_etappi->Draw("same Chist" );

  TCanvas *c21 = plot.Canvas ("c21",400,300,410,610,1);
  TH1D *ptr1 = plot.Frame("ptr1","pt","ratio",pt_low,pt_high*.999,0.002,1.2);
  ptr1->Draw();
  h_g_pih->Draw("same Chist");
  h_g_etah->Draw("same Chist");
  h_g_omegah->Draw("same Chist");
  h_g_etaph->Draw("same Chist");


  TCanvas *c3 = plot.Canvas ("c3",600,600,810,10,1);
  TH1D *mass = plot.Frame("mass","mass","counts",0.,1.5,3e-12,3e-3);
  mass->Draw();
  h_masspi->Draw("same Chist");
  h_masseta->Draw("same Chist");
  h_massetap->Draw("same Chist");
  h_massomega->Draw("same Chist");
  h_massrho0->Draw("same Chist");
  h_massphi->Draw("same Chist");

  h_mass->Draw("same Chist");
}




