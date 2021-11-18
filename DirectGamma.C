#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
// 
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TLorentzVector.h>

using namespace std;

#include "HELIOSLibrary/HELIOSLibrary.h"
#include "MyPlotting/MyPlot.C"

PHENIXDetector MyPHENIX;                       // define PHENIX detector 

TLorentzVector RecoPhoton(TLorentzVector photon, Int_t id = photonID){        // use fast simulation to reconstruct photon 
    MyPHENIX.CharacterizeTrack(photon, 0., id); 
    TLorentzVector reco = MyPHENIX.ReconstructShower(photon, id);  
//    reco.Print();
    return reco;
}

Double_t fitf(Double_t *x,Double_t *par) {
      Double_t fitval = par[0];
      return fitval;
}

///////////////////////////////////////////////////////////////////////////////
void DirectGamma(){
//
  TRandom3 MyRandy = TRandom3(0);                // Random Generator
                                
  MyPlot plot;

  Particle pi0("pi0");
  Particle eta("eta");
  Particle omega("omega");
  Particle etap("etap");
  Particle photon("photon");
  TLorentzVector gamma;
  TLorentzVector electron,positron;   
  TLorentzVector temp, temp1, temp2;

  Double_t pt_min = 0.3;     // used for all simulation
  Double_t pt_max = 25.;
  Double_t pt_low = 0.0;     // used for plotting
  Double_t pt_high = 25.;
  Double_t ecut = 0.5;
  Double_t opacut = 0;
  Int_t nevt = 10000000;
  Int_t ndecay =0;
  Double_t ww = 0, ww1, ww2; 
  Int_t id = 0;

  Bool_t il = true;          // include leptons
  Bool_t reco = false;       // reconstruct as photon
  Bool_t VTXconv = false;    // flag that photon converted

//// momentum distributions dN/dpt for different particles
// ppg088 parameters for pi0 pp, all other mt scaled to pi0
  TF1 *piHagedorn     = Hagedorn("piHagedorn", pi0Mass, pt_max, pt_min, 377., 0.356, 0.068, 0.7, -8.25);
  TF1 *etaHagedorn    = Hagedorn("etaHagedorn", etaMass, pt_max, pt_min,  377., 0.356, 0.068, 0.7, -8.25);
  TF1 *omegaHagedorn  = Hagedorn("omegaHagedorn", omegaMass, pt_max, pt_min, 377., 0.356, 0.068, 0.7, -8.25);
  TF1 *etapHagedorn   = Hagedorn("etapHagedorn", etapMass, pt_max, pt_min,  377., 0.356, 0.068, 0.7, -8.25);
// direct photon fit fom PPG140
  TF1 *directPhoton   = PromptPhotonYield("directPhoton", pt_max, pt_min);
  Read_GPR_pion_pp200();
  Read_GPR_etapi_uni();


// histograms for parent particle dN/dpt
  TH1D *h_ptpi        = new TH1D("h_ptpi","pt",100,pt_low,pt_high);
  TH1D *h_pteta       = new TH1D("h_pteta","",100,pt_low,pt_high);
  TH1D *h_ptomega     = new TH1D("h_ptomega","",100,pt_low,pt_high);
  TH1D *h_ptetap      = new TH1D("h_ptetap","",100,pt_low,pt_high);
  TH1D *h_ptdirect    = new TH1D("h_ptdirect","pt",100,pt_low,pt_high);
  plot.StyleMe(h_ptdirect,   20, kMagenta, 1., 1, 2); 
  plot.StyleMe(h_ptpi,   20, kBlue, 1., 1, 2); 
  plot.StyleMe(h_pteta,   20, kRed, 1., 1, 2); 
  plot.StyleMe(h_ptomega,   20, kGreen+2, 1., 1, 2); 
  plot.StyleMe(h_ptetap,   20, kOrange+2, 1., 1, 2);   

// dito but in in variant cross section
  TH1D *h_ptdirectx   = new TH1D("h_ptdirectx","pt",100,pt_low,pt_high);
  TH1D *h_ptpix       = new TH1D("h_ptpix","pt",100,pt_low,pt_high);
  plot.StyleMe(h_ptdirectx,   20, kMagenta, 1., 1, 2); 
  plot.StyleMe(h_ptpix,   20, kBlue, 1., 1, 2); 

// decay photon spectra from parent particles
  TH1D *h_pt_gpi      = new TH1D("h_pt_gpi","pt_g",100,pt_low,pt_high);
  TH1D *h_pt_geta     = new TH1D("h_pt_geta","",100,pt_low,pt_high);
  TH1D *h_pt_gomega   = new TH1D("h_pt_gomega","",100,pt_low,pt_high);
  TH1D *h_pt_getap    = new TH1D("h_pt_getap","",100,pt_low,pt_high);
  plot.StyleMe(h_pt_gpi, 20, kBlue, 1., 2, 2); 
  plot.StyleMe(h_pt_geta, 20, kRed, 1., 2, 2); 
  plot.StyleMe(h_pt_gomega, 20, kGreen+2, 1., 2, 2); 
  plot.StyleMe(h_pt_getap, 20, kOrange+2, 1., 2, 2); 

// check opening angle distribution for merging effects
  TH2D *h_opa_pigg   = new TH2D("h_opa_pigg","",100,pt_low,pt_high,100,0.,50.);
  TH2D *h_opa_etagg   = new TH2D("h_opa_etagg","",100,pt_low,pt_high,100,0.,50.);
  TH2D *h_ptpt_pigg   = new TH2D("h_ptpt_pigg","",100,pt_low,pt_high,100,pt_low,pt_high);
  TH2D *h_ptpt_etagg   = new TH2D("h_ptpt_etagg","",100,pt_low,pt_high,100,pt_low,pt_high);
  TH3D *h_opaptpt_pigg = new TH3D("h_opaptpt_pigg","",100,pt_low,pt_high,100,pt_low,pt_high,100,0.,50.);


// decay photons in PHENIX 
  TH1D *h_pt_gpi_r         = new TH1D("h_pt_gpi_r","pt_g_r",100,pt_low,pt_high);
  TH1D *h_pt_geta_r         = new TH1D("h_pt_geta_r","",100,pt_low,pt_high);
  TH1D *h_pt_gomega_r         = new TH1D("h_pt_gomega_r","",100,pt_low,pt_high);
  TH1D *h_pt_getap_r         = new TH1D("h_pt_getap_r","",100,pt_low,pt_high);
  plot.StyleMe(h_pt_gpi_r, 20, kBlue, 1., 2, 2); 
  plot.StyleMe(h_pt_geta_r, 20, kRed, 1., 2, 2); 
  plot.StyleMe(h_pt_gomega_r, 20, kGreen+2, 1., 2, 2); 
  plot.StyleMe(h_pt_getap_r, 20, kOrange+2, 1., 2, 2); 

// all phtons from hadron decays
  TH1D *h_pt_ghadron           = new TH1D("h_pt_ghadron","pt_g",100,pt_low,pt_high);
  plot.StyleMe(h_pt_ghadron,   20, kBlack, 1., 1, 2); 

// systematic study
// eta/pi uncertainty
  TH1D *h_gh_sys_10         = new TH1D("h_gh_sys_10","pt_g",100,pt_low,pt_high);
  TH1D *h_gh_sys_11         = new TH1D("h_gh_sys_10","pt_g",100,pt_low,pt_high);
  plot.StyleMe(h_gh_sys_10,   20, kBlack, 1., 2, 2); 
  plot.StyleMe(h_gh_sys_11,   20, kBlack, 1., 2, 2); 
// omega & eta'
  TH1D *h_gh_sys_20         = new TH1D("h_gh_sys_20","pt_g",100,0.5,pt_high);
  TH1D *h_gh_sys_21         = new TH1D("h_gh_sys_20","pt_g",100,0.5,pt_high);
  plot.StyleMe(h_gh_sys_20,   20, kBlack, 1., 2, 2); 
  plot.StyleMe(h_gh_sys_21,   20, kBlack, 1., 2, 2); 

// dito reconstructed
  TH1D *h_pt_ghadron_r           = new TH1D("h_pt_ghadron_r","pt_g",100,pt_low,pt_high);
  plot.StyleMe(h_pt_ghadron_r,   20, kBlack, 1., 1, 2); 

//all photoons
  TH1D *h_pt_gincl             = new TH1D("h_pt_gincl","pt_g",100,pt_low,pt_high);
  plot.StyleMe(h_pt_gincl,   20, kCyan, 1., 1, 2); 


///////////// event loop /////////////////////////////////////////////////////////////////////////////////
  cout << endl;
  for (int i=0; i<nevt; i++){                                      // event loop
    if(i%100000==0) cout << "event loop at " << i << endl;         // event counter printed 

///////////////////////////////////////// direct photons //////////////////////////////////////////////////
    photon.ResetP();                                                // reset monemtum   
    photon.GenerateP(pt_min,pt_max);                                // generate new momentum flat in pt
    photon.SetWeight(directPhoton->Eval(photon.Pt())/float(nevt));  // set weight 
    h_ptdirect->Fill(photon.Pt(),photon.Weight());                  
    h_ptdirectx->Fill(photon.Pt(),photon.Weight()/photon.Pt());     
    h_pt_gincl->Fill(photon.Pt(),photon.Weight());     

    if (i<10)  cout << "----- Direct Photon --------------------------------------------------------" << endl; 
    if (i<10) photon.Print(); 
    if (i<10) cout << " pt " << photon.Pt() << " weight " << photon.Weight() << endl;

/////////////////////////////////////////////// pi decays ////////////////////////////////////////////////
    pi0.ResetP();                                                 // reset
    pi0.GenerateP(pt_min,pt_max);                                 // generate
//    pi0.SetWeight(piHagedorn->Eval(pi0.Pt())/float(nevt));        // set weight
    pi0.SetWeight(pi0.Pt()*Weight_GPR_pion_pp200(pi0.Pt())/float(nevt));
    h_ptpix->Fill(pi0.Pt(),pi0.Weight()/pi0.Pt());                 // histogram 
    h_ptpi->Fill(pi0.Pt(),pi0.Weight());                         // 

    pi0.DecayFlat();                                              // generate decay photons
    ndecay = pi0.GetNumberOfDaughters();                          // 2 for pi0->gg and 3 for pi0-e+e-g

    if (i<10)  cout << "----- Pi0Deacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 
    if (i<100) cout << "pt pi0 " << pi0.Pt() << " gamma pt " << pi0.GetDecayDaughter(0).Pt() << " weight " << pi0.GetDaughterWeight(0) << " " << pi0.Weight() << " " << piHagedorn->Eval(pi0.Pt())/float(nevt) << endl;
    if (i<100) cout << "pt pi0 " << pi0.Pt() << " gamma pt " << pi0.GetDecayDaughter(1).Pt() << " weight " << pi0.GetDaughterWeight(1) << " " << pi0.Weight() << " " << piHagedorn->Eval(pi0.Pt())/float(nevt) << endl;

    

    Double_t opa=100.;
    if(ndecay ==2 && pi0.GetDaughterID(0)==pi0.GetDaughterID(1) && pi0.GetDaughterID(0)==photonID){
        ww   = pi0.GetDaughterWeight(0)*pi0.Weight(); 
        temp1 = pi0.GetDecayDaughter(0); 
        temp2 = pi0.GetDecayDaughter(1); 
        opa = temp1.Angle(temp2.Vect())*P_R_EMCal*100.; 
        h_opa_pigg->Fill(temp1.Pt(),opa,ww);
        h_opa_pigg->Fill(temp2.Pt(),opa,ww);
        h_ptpt_pigg->Fill(temp1.Pt(),temp2.Pt(),ww);
        h_opaptpt_pigg->Fill(temp1.Pt(),temp2.Pt(),opa,ww);
    }

    for (Int_t j=0; j< ndecay; j++) {                             // loop over decay particles
      temp = pi0.GetDecayDaughter(j);                             // temporary TLorentzVector of daughter
      ww   = pi0.GetDaughterWeight(j)*pi0.Weight();               // update weight with BR 

      id   = pi0.GetDaughterID(j);                                // get daugther id
      if (i<10) temp.Print();
     
      VTXconv = false; 
      if (id == photonID) {                                       // if its a photon 
        h_pt_gpi->Fill(temp.Pt(),ww);                             // histogram "true" information
        h_pt_ghadron->Fill(temp.Pt(),ww);
        h_gh_sys_10->Fill(temp.Pt(),ww);
        h_gh_sys_11->Fill(temp.Pt(),ww);
        h_gh_sys_20->Fill(temp.Pt(),ww);
        h_gh_sys_21->Fill(temp.Pt(),ww);
         h_pt_gincl->Fill(temp.Pt(),ww);
        VTXconv = (MyPHENIX.VTXConversion()>0);                     // check if photon will convert 
      }  
//                                                                 check if decay daughter needs to be reconstrcuted
      reco = (!VTXconv && id==photonID) or (il && (id == electronID or id == positronID));
      if  (reco && opa>opacut) {  
        gamma = RecoPhoton(temp,id);                              // fast simulation for reconstruction          
        if (i<10) cout << " reconstructing decay daughter " << endl;
        if (i<10) gamma.Print(); 
        if (gamma.E()>ecut) {                                     // pass basic energy cut
            h_pt_gpi_r->Fill(gamma.Pt(),ww);     
            h_pt_ghadron_r->Fill(gamma.Pt(),ww);
        } 
      } else if (reco && opa<=opacut) {
        temp = temp1+temp2;
        gamma = RecoPhoton(temp,id);
        if (gamma.E()>ecut) {                                     // pass basic energy cut
            h_pt_gpi_r->Fill(gamma.Pt(),ww/2);     
            h_pt_ghadron_r->Fill(gamma.Pt(),ww/2);
        }   
      } 
      if (VTXconv){                                               // if photon converts reconstruct e+e-
        PhotonConvert(temp, temp1, temp2);                        // as photons 
        electron = RecoPhoton(temp1,electronID); 
        positron = RecoPhoton(temp2,positronID);
        if (i<10) cout << " reconstructing photon conversion " << endl;
        if (electron.E()>ecut){
            if (i<100) electron.Print(); 
            h_pt_gpi_r->Fill(electron.Pt(),ww);     
            h_pt_ghadron_r->Fill(electron.Pt(),ww);            
        } 
        if (positron.E()>ecut){
            if (i<100) positron.Print(); 
            h_pt_gpi_r->Fill(positron.Pt(),ww);     
            h_pt_ghadron_r->Fill(positron.Pt(),ww);            
        }      
      }
    }
/////////////////////////////////////////////// eta decays ////////////////////////////////////////////////
    eta.ResetP();                                                // dito for eta
    eta.GenerateP(pt_min,pt_max);                                // 
//    eta.SetWeight(etaHagedorn->Eval(eta.Pt())/float(nevt)*Eta_to_Pi0);
    eta.SetWeight(eta.Pt()*Weight_GPR_pion_pp200(eta.Pt())/float(nevt) * Weight_GPR_etapi_uni(eta.Pt(),0));
    h_pteta->Fill(eta.Pt(),eta.Weight()/eta.Pt());                         // histogram pion proprety
      
    eta.DecayFlat();
    ndecay = eta.GetNumberOfDaughters();

    if (i<10)  cout << "----- EtaDeacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 

    if(ndecay ==2 && eta.GetDaughterID(0)==eta.GetDaughterID(1) && eta.GetDaughterID(0)==photonID){
        ww   = eta.GetDaughterWeight(0)*eta.Weight(); 
        temp1 = eta.GetDecayDaughter(0); 
        temp2 = eta.GetDecayDaughter(1); 
        Double_t a = temp1.Angle(temp2.Vect())*P_R_EMCal*100.; 
        h_opa_etagg->Fill(temp1.Pt(),a,ww);
        h_opa_etagg->Fill(temp2.Pt(),a,ww);
        h_ptpt_etagg->Fill(temp1.Pt(),temp2.Pt(),ww);
    }
     
    for (Int_t j=0; j< ndecay; j++) {
      temp = eta.GetDecayDaughter(j);
      if (i<10) temp.Print();
      ww   = eta.GetDaughterWeight(j)*eta.Weight();
      ww1  = eta.GetDaughterWeight(j)*eta.Pt()*Weight_GPR_pion_pp200(eta.Pt())/float(nevt) 
                  * Weight_GPR_etapi_uni(eta.Pt(),-1.);
      ww2  = eta.GetDaughterWeight(j)*eta.Pt()*Weight_GPR_pion_pp200(eta.Pt())/float(nevt) 
                  * Weight_GPR_etapi_uni(eta.Pt(),1.);
      id   = eta.GetDaughterID(j);

      VTXconv = false; 
      if (id == photonID) {                                       // if its a photon 
        h_pt_geta->Fill(temp.Pt(),ww);                             // histogram "true" information
        h_pt_ghadron->Fill(temp.Pt(),ww);
        h_gh_sys_10->Fill(temp.Pt(),ww1);
        h_gh_sys_11->Fill(temp.Pt(),ww2);
        h_gh_sys_20->Fill(temp.Pt(),ww);
        h_gh_sys_21->Fill(temp.Pt(),ww);
         h_pt_gincl->Fill(temp.Pt(),ww);
        VTXconv = (MyPHENIX.VTXConversion()>0);                     // check if photon will convert 
      }  
                                                                  // check if decay daughter needs to be reconstrcuted
      reco = (!VTXconv && id==photonID) or (il && (id == electronID or id == positronID));
      if  (reco) {  
        gamma = RecoPhoton(temp,id);                              // fast simulation for reconstruction          
        if (i<10) cout << " reconstructing decay daughter " << endl;
        if (i<10) gamma.Print(); 
        if (gamma.E()>ecut) {                                     // pass basic energy cut
            h_pt_geta_r->Fill(gamma.Pt(),ww);     
            h_pt_ghadron_r->Fill(gamma.Pt(),ww);
        }
      } 
      if (VTXconv){                                               // if photon converts reconstruct e+e-
        PhotonConvert(temp, temp1, temp2);                        // as photons 
        electron = RecoPhoton(temp1,electronID); 
        positron = RecoPhoton(temp2,positronID);
        if (i<10) cout << " reconstructing photon conversion " << endl;
        if (electron.E()>ecut){
            if (i<100) electron.Print(); 
            h_pt_geta_r->Fill(electron.Pt(),ww);     
            h_pt_ghadron_r->Fill(electron.Pt(),ww);            
        } 
        if (positron.E()>ecut){
            if (i<100) positron.Print(); 
            h_pt_geta_r->Fill(positron.Pt(),ww);     
            h_pt_ghadron_r->Fill(positron.Pt(),ww);            
        }      
      } 

    }
/////////////////////////////////////////////// omega decays ////////////////////////////////////////////////
    omega.ResetP();                                                // dito for omega 
    omega.GenerateP(pt_min,pt_max);                                // 
//    omega.SetWeight(omegaHagedorn->Eval(omega.Pt())/float(nevt)*Omega_to_Pi0);
    Double_t pt_scaled = Pt_Mt_Scaled(omega.Pt(), omegaMass, etaMass);
    omega.SetWeight( pt_scaled*Weight_GPR_pion_pp200(pt_scaled)/float(nevt) 
          * Weight_GPR_etapi_uni(pt_scaled,0)*Omega_to_Pi0/Eta_to_Pi0);
     h_ptomega->Fill(omega.Pt(),omega.Weight()/omega.Pt());         // 
    
    omega.DecayFlat();                        // only one decay branch with gamma
    ndecay = omega.GetNumberOfDaughters();

    if (i<10)  cout << "----- omegaDeacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 
   
    for (Int_t j=0; j< ndecay; j++) {
       temp = omega.GetDecayDaughter(j);
       if (i<10) temp.Print();
       ww   = omega.GetDaughterWeight(j)*omega.Weight();
       ww1 = ww*(1+sqrt(0.111*0.111+0.025*0.025+0.026*0.026));   // based on omega/pi = 0.81+/-0.02(stat)+/-0.09(sys)
       ww2 = ww/(1+sqrt(0.111*0.111+0.025*0.025+0.026*0.026));   // and BR 8.4+/-0.22%
       id   = omega.GetDaughterID(j);

      VTXconv = false; 
      if (id == photonID) {                                       // if its a photon 
        h_pt_gomega->Fill(temp.Pt(),ww);                             // histogram "true" information
        h_pt_ghadron->Fill(temp.Pt(),ww);
        h_gh_sys_10->Fill(temp.Pt(),ww);
        h_gh_sys_11->Fill(temp.Pt(),ww);
        h_gh_sys_20->Fill(temp.Pt(),ww1);
        h_gh_sys_21->Fill(temp.Pt(),ww2);
        h_pt_gincl->Fill(temp.Pt(),ww);
        VTXconv = (MyPHENIX.VTXConversion()>0);                     // check if photon will convert 
      }  
                                                                  // check if decay daughter needs to be reconstrcuted
      reco = (!VTXconv && id==photonID) or (il && (id == electronID or id == positronID));
      if  (reco) {  
        gamma = RecoPhoton(temp,id);                              // fast simulation for reconstruction          
        if (i<10) cout << " reconstructing decay daughter " << endl;
        if (i<10) gamma.Print(); 
        if (gamma.E()>ecut) {                                     // pass basic energy cut
            h_pt_gomega_r->Fill(gamma.Pt(),ww);     
            h_pt_ghadron_r->Fill(gamma.Pt(),ww);
        }
      } 
      if (VTXconv){                                               // if photon converts reconstruct e+e-
        PhotonConvert(temp, temp1, temp2);                        // as photons 
        electron = RecoPhoton(temp1,electronID); 
        positron = RecoPhoton(temp2,positronID);
        if (i<10) cout << " reconstructing photon conversion " << endl;
        if (electron.E()>ecut){
            if (i<100) electron.Print(); 
            h_pt_gomega_r->Fill(electron.Pt(),ww);     
            h_pt_ghadron_r->Fill(electron.Pt(),ww);            
        } 
        if (positron.E()>ecut){
            if (i<100) positron.Print(); 
            h_pt_gomega_r->Fill(positron.Pt(),ww);     
            h_pt_ghadron_r->Fill(positron.Pt(),ww);            
        }      
      }
    }   
/////////////////////////////////////////////// eta' decays ////////////////////////////////////////////////
    etap.ResetP();                                                // dito for eta'
    etap.GenerateP(pt_min,pt_max);                                // 
//    etap.SetWeight(etapHagedorn->Eval(etap.Pt())/float(nevt)*Etap_to_Pi0);
    pt_scaled = Pt_Mt_Scaled(etap.Pt(), etapMass, etaMass);
    etap.SetWeight( pt_scaled*Weight_GPR_pion_pp200(pt_scaled)/float(nevt) 
          * Weight_GPR_etapi_uni(pt_scaled,0)*Etap_to_Pi0/Eta_to_Pi0);
    h_ptetap->Fill(etap.Pt(),etap.Weight()/etap.Pt());            // 
      
    etap.DecayFlat();
    ndecay = etap.GetNumberOfDaughters();

    if (i<10)  cout << "----- etapDeacy --------------------------------------------------------" << endl; 
    if (i<10) cout << " found decay with " << ndecay << " daughters " << endl; 
     
    for (Int_t j=0; j< ndecay; j++) {
      temp = etap.GetDecayDaughter(j);
      if (i<10) temp.Print();
      ww   = etap.GetDaughterWeight(j)*etap.Weight();
       ww1 = ww*(1+sqrt(0.3*0.3+0.014*0.014));   // based on etap/pi = 0.25+/-0.075(sys)
       ww2 = ww/(1+sqrt(0.3*0.3+0.014*0.014));   // and BR 2.37+/-0.033%
      id   = etap.GetDaughterID(j);

      VTXconv = false; 
      if (id == photonID) {                                       // if its a photon 
        h_pt_getap->Fill(temp.Pt(),ww);                             // histogram "true" information
        h_pt_ghadron->Fill(temp.Pt(),ww);
        h_gh_sys_10->Fill(temp.Pt(),ww);
        h_gh_sys_11->Fill(temp.Pt(),ww);
        h_gh_sys_20->Fill(temp.Pt(),ww1);
        h_gh_sys_21->Fill(temp.Pt(),ww2);
         h_pt_gincl->Fill(temp.Pt(),ww);
        VTXconv = (MyPHENIX.VTXConversion()>0);                     // check if photon will convert 
      }  
                                                                  // check if decay daughter needs to be reconstrcuted
      reco = (!VTXconv && id==photonID) or (il && (id == electronID or id == positronID));
      if  (reco) {  
        gamma = RecoPhoton(temp,id);                              // fast simulation for reconstruction          
        if (i<10) cout << " reconstructing decay daughter " << endl;
        if (i<10) gamma.Print(); 
        if (gamma.E()>ecut) {                                     // pass basic energy cut
            h_pt_getap_r->Fill(gamma.Pt(),ww);     
            h_pt_ghadron_r->Fill(gamma.Pt(),ww);
        }
      } 
      if (VTXconv){                                               // if photon converts reconstruct e+e-
        PhotonConvert(temp, temp1, temp2);                        // as photons 
        electron = RecoPhoton(temp1,electronID); 
        positron = RecoPhoton(temp2,positronID);
        if (i<10) cout << " reconstructing photon conversion " << endl;
        if (electron.E()>ecut){
            if (i<100) electron.Print(); 
            h_pt_getap_r->Fill(electron.Pt(),ww);     
            h_pt_ghadron_r->Fill(electron.Pt(),ww);            
        } 
        if (positron.E()>ecut){
            if (i<100) positron.Print(); 
            h_pt_getap_r->Fill(positron.Pt(),ww);     
            h_pt_ghadron_r->Fill(positron.Pt(),ww);            
        }      
      }
    } 
}

///////////////////// end of event loop ///////////////////////////////////////////////////////
// Calculate some ratios of histograms for plotting

// parent to pi0 ratios
  TH1D *h_etapi       = (TH1D*)h_pteta->Clone("h_etapi");
  TH1D *h_omegapi     = (TH1D*)h_ptomega->Clone("h_omegapi");
  TH1D *h_etappi      = (TH1D*)h_ptetap->Clone("h_etappi");
  h_etapi->Divide(h_ptpix);
  h_omegapi->Divide(h_ptpix);
  h_etappi->Divide(h_ptpix);

// ratio for decay photons to all photons from hadron decays
  TH1D *h_g_pih       = (TH1D*)h_pt_gpi->Clone("h_g_pih");
  TH1D *h_g_etah      = (TH1D*)h_pt_geta->Clone("h_g_etah");
  TH1D *h_g_omegah    = (TH1D*)h_pt_gomega->Clone("h_g_omegah");
  TH1D *h_g_etaph     = (TH1D*)h_pt_getap->Clone("h_g_etaph");
  h_g_pih->Divide(h_pt_ghadron);
  h_g_etah->Divide(h_pt_ghadron);
  h_g_omegah->Divide(h_pt_ghadron);
  h_g_etaph->Divide(h_pt_ghadron);
  plot.StyleMe(h_g_pih, 20, kBlue, 1.2, 1., 2.);
  plot.StyleMe(h_g_etah, 20, kRed, 1.2, 1., 2.);
  plot.StyleMe(h_g_omegah, 20, kGreen+2, 1.2, 1., 2.);
  plot.StyleMe(h_g_etaph, 20, kOrange+2, 1.2, 1., 2.);

// dito but reconstructed in PHENIX
  TH1D *h_g_pih_r       = (TH1D*)h_pt_gpi_r->Clone("h_g_pih_r");
  TH1D *h_g_etah_r      = (TH1D*)h_pt_geta_r->Clone("h_g_etah_r");
  TH1D *h_g_omegah_r    = (TH1D*)h_pt_gomega_r->Clone("h_g_omegah_r");
  TH1D *h_g_etaph_r     = (TH1D*)h_pt_getap_r->Clone("h_g_etaph_r");
  h_g_pih_r->Divide(h_pt_ghadron_r);
  h_g_etah_r->Divide(h_pt_ghadron_r);
  h_g_omegah_r->Divide(h_pt_ghadron_r);
  h_g_etaph_r->Divide(h_pt_ghadron_r);

// gamma hadron / gamma pi0
  TH1D *h_g_htopi     = (TH1D*) h_pt_ghadron->Clone("h_g_htopi");
  h_g_htopi->Divide(h_pt_gpi);
  TH1D *h_g_htopi_10  = (TH1D*) h_gh_sys_10->Clone("h_g_htopi_10");
  h_g_htopi_10->Divide(h_pt_gpi);
  TH1D *h_g_htopi_11  = (TH1D*) h_gh_sys_11->Clone("h_g_htopi_11");
  h_g_htopi_11->Divide(h_pt_gpi);

  TH1D *h_g_htopi_r     = (TH1D*) h_pt_ghadron_r->Clone("h_g_htopi_r");
  h_g_htopi_r->Divide(h_pt_gpi_r);

// gamma eta+omega+etap / gamma eta
  TH1D *h_g_htoeta = (TH1D*) h_pt_ghadron->Clone("h_g_htoeta");
  h_g_htoeta->Add(h_pt_gpi,-1.);
  h_g_htoeta->Divide(h_pt_geta);
  TH1D *h_g_htoeta_20 = (TH1D*) h_gh_sys_20->Clone("h_g_htoeta_20");
  h_g_htoeta_20->Add(h_pt_gpi,-1.);
  h_g_htoeta_20->Divide(h_pt_geta);
  TH1D *h_g_htoeta_21 = (TH1D*) h_gh_sys_21->Clone("h_g_htoeta_21");
  h_g_htoeta_21->Add(h_pt_gpi,-1.);
  h_g_htoeta_21->Divide(h_pt_geta);

// direct gamma / pi0
  TH1D *h_gtopi       = (TH1D*) h_ptdirectx->Clone("h_gtopi");
  h_gtopi->Divide(h_ptpix);

  TF1 *fhtopi = new TF1("htopi", fitf, 5., 25., 1); 
  h_g_htopi->Fit("htopi","R");
  TF1 *fhtoeta = new TF1("htoeta", fitf, 8., 23., 1); 
  h_g_htoeta->Fit("htoeta","R");
  TF1 *fhtoeta1 = new TF1("htoeta1", fitf, 8., 23., 1); 
  h_g_htoeta_20->Fit("htoeta1","R");
  TF1 *fhtoeta2 = new TF1("htoeta2", fitf, 8., 23., 1); 
  h_g_htoeta_21->Fit("htoeta2","R");

//////////////////// plot histograms //////////////////////////////////////////////////////
// parent pt spectra
  plot.SetLeftMargin(.2);
  TCanvas *c1 = plot.Canvas ("c1",400,500,10,10,1);
  plot.SetyTitleOffset(1.5);
  TH1D *f1 = plot.Frame("f1","p_{T} [GeV/c]","Ed^{3}#sigma/dp^{3} [mb GeV^{-2} c^{3}]",pt_low,25,1e-14,1.);
//  TH1D *f1 = plot.Frame("f1","p_{T} [GeV/c]","Ed^{3}#sigma/dp^{3} [mb GeV^{-2} c^{3}]",pt_low,pt_high*.999,1e-14,1.);
  f1->Draw();
  h_ptdirectx->Draw("sameL Chist");
  h_ptpix->Draw("sameL Chist");
  h_pteta->Draw("sameL Chist");
  h_ptomega->Draw("sameL Chist");
  h_ptetap->Draw("sameL Chist");

  TLegend *L1 = plot.Legend("p+p #sqrt{s} 200 GeV ",0.5,.67,.9,.95);
  L1->AddEntry(h_ptpix,"#pi^{0}","l");
  L1->AddEntry(h_pteta,"#eta","l");
  L1->AddEntry(h_ptomega,"#omega","l");
  L1->AddEntry(h_ptetap,"#eta'","l");
  L1->AddEntry(h_ptdirect,"direct #gamma","l");
  L1->Draw("same");
  plot.Reset();

// parent to pi0 ratios
  TCanvas *c2 = plot.Canvas ("c2",500,400,10,510);
  TH1D *f2 = plot.Frame("f2","p_{T} [GeV/c]"," hadron / #pi^{0}",pt_low,pt_high*.999,0.0,1.509);
  f2->Draw();
  h_etapi->Draw("same Chist" );
  h_omegapi->Draw("same Chist" );
  h_etappi->Draw("same Chist" );

  plot.SetLegendSize(0.06);
  plot.SetLegendColor(kGreen+2);
  TLegend *L2a = plot.Legend("#omega",0.85,.76,.9,.81);
  L2a->Draw("same");
  plot.SetLegendColor(kRed);
  TLegend *L2b = plot.Legend("#eta",0.85,.45,.9,.50);
  L2b->Draw("same");
  plot.SetLegendColor(kOrange+2);
  TLegend *L2c = plot.Legend("#eta'",0.85,.25,.9,.30);
  L2c->Draw("same");
  plot.Reset();

// photon pt spectra
  plot.SetLeftMargin(.2);
  TCanvas *c3 = plot.Canvas ("c3",400,500,410,10,1);
  plot.SetyTitleOffset(1.5);
  TH1D *f3 = plot.Frame("f3","p_{T} [GeV/c]","dN_{#gamma}/dp_{T} ",pt_low,pt_high*.999,1e-12,1.);
  f3->Draw();
//  h_pt_gincl->Draw("sameL Chist");
  h_ptdirect->Draw("same Chist");
//  h_pt_ghadron->Draw("same Chist");
  h_pt_gpi->Draw("same Chist");

  TLegend *L3 = plot.Legend("#gamma from p+p #sqrt{s} 200 GeV ",0.5,.77,.9,.95);
  L3->AddEntry(h_pt_gincl,"inclusive","l");
  L3->AddEntry(h_pt_ghadron,"hadron decays","l");
  L3->AddEntry(h_ptdirect,"direct","l");
  L3->Draw("same");
  plot.Reset();

// gamma decay to gamma hadron 
  TCanvas *c4 = plot.Canvas ("c4",500,400,510,510,1);
  TH1D *f4 = plot.Frame("f4","p_{T} [GeV/c]"," #gamma_{decay} / #gamma_{hadron}",pt_low,pt_high*.999,0.002,5.);
  f4->Draw();
  h_g_pih->Draw("same Chist");
  h_g_etah->Draw("same Chist");
  h_g_omegah->Draw("same Chist");
  h_g_etaph->Draw("same Chist");
  h_g_pih_r->Draw("same Chist");
  h_g_etah_r->Draw("same Chist");
  h_g_omegah_r->Draw("same Chist");
  h_g_etaph_r->Draw("same Chist");

  TLegend *L4 = plot.Legend("", .2, .82, .7, .96);
  L4->AddEntry(h_g_pih,"true distributions", "l");
  L4->AddEntry(h_g_pih_r,"reconstructed, with e^{+}e^{-} and conversions", "l");
  L4->Draw("same");
  plot.SetLegendColor(kBlue);
  plot.SetLegendSize(0.06);
  TLegend *L4a = plot.Legend("#pi^{0}",0.85,.7,.9,.75);
  L4a->Draw("same");
  plot.SetLegendColor(kRed);
  TLegend *L4b = plot.Legend("#eta",0.825,.53,.9,.58);
  L4b->Draw("same");
  plot.SetLegendColor(kGreen+2);
  TLegend *L4c = plot.Legend("#omega",0.875,.46,.9,.51);
  L4c->Draw("same");
  plot.SetLegendColor(kOrange+2);
  TLegend *L4d = plot.Legend("#eta'",0.85,.29,.9,.34);
  L4d->Draw("same");
  plot.Reset();

// gamma direct to pi0  
  TCanvas *c5 = plot.Canvas ("c5",500,400,1010,10);
  TH1D *f5 = plot.Frame("f5","p_{T} [GeV/c]","#gamma_{direct} / #pi^{0}",pt_low,pt_high*.999,0.0,0.375);
  f5->Draw();
  h_gtopi->Draw("same Chist");

// gamma hadron / gamma pi0
  TCanvas *c6 = plot.Canvas ("c6",500,400,1010,410);
  TH1D *f6 = plot.Frame("f6","p_{T} [GeV/c]","#gamma_{h} / #gamma_{#pi^{0}}",pt_low,21.99,0.95,1.359);
  f6->Draw();
  h_g_htopi->Draw("same Chist");
  h_g_htopi_10->Draw("same Chist");
  h_g_htopi_11->Draw("same Chist");
//  fhtopi->Draw("same");

  plot.SetLegendSize(0.06);
  TLegend *L6 = plot.Legend("#gamma_{hadr} / #gamma_{#pi^{0}} (p_{T}>8 GeV) = 1.23",0.35,.33,.9,.38);
//  L6->Draw("same");

// gamma hadron / gamma pi0
  TCanvas *c7 = plot.Canvas ("c7",500,400,510,10);
  TH1D *f7 = plot.Frame("f7","p_{T} [GeV/c]","#gamma_{#eta,#omega,#eta'} / #gamma_{#eta}",pt_low,22.99,0.9601,1.319);
  f7->Draw();
  h_g_htoeta->Draw("same Chist");
  h_g_htoeta_20->Draw("same Chist");
  h_g_htoeta_21->Draw("same Chist");
  fhtoeta->Draw("same");

  plot.SetLegendSize(0.06);
  TLegend *L7 = plot.Legend("#gamma_{#eta,#omega,#eta'} / #gamma_{#eta} (p_{T}>8 GeV) = 1.19#pm0.03",0.33,.33,.9,.38);
  L7->Draw("same");

}




