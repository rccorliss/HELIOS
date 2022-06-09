// Writes out ROOT output file 
// This example generates single pi0 events with the pi0 decay handled by HELIOS


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TString.h>

using namespace std;

#include "HELIOSLibrary/HELIOSLibrary.h"
PHENIXDetector MyPHENIX; 

void setTrack(WriteTrack& newTrack, int isFinal, int num, int id, int ist, float px, float py, float pz, float en, float mass, float xpos, float ypos, float zpos, int br, int parent_index, float weight);
void WriteROOT2Oscar(TString infile = "input.root", TString output = "oscar.dat", Int_t nstart = 0, Int_t nend = 1);

void WriteROOTOutput(TString OutputName){

	TRandom3 MyRandy = TRandom3(0);

	TFile* f = new TFile(OutputName,"RECREATE");

	WriteEvent MyEvent;
	WriteTrack MyTrack;


	TH2F* h_input = new TH2F("h_input", "PID vs decay", 2500, -0.5, 2499.5, 10, -0.5, 9.5);

	TH1F* h_input_pt = new TH1F("h_input_pt", "pt profile of input", 200, 0., 30.0);
	TH1F* h_input_eta = new TH1F("h_input_eta", "eta profile of input", 200, -1.0, 1.0);
	TH1F* h_input_phi = new TH1F("h_input_phi", "phi profile of input", 200, -10.0, 10.0);
	

    TTree* T = new TTree("T","Tree with my events");
    T->Branch("MyEvent",&MyEvent);


	TLorentzVector temp, temp1, temp2;
	Int_t ndecay = 0;
	TString DecayBranch;
	Int_t DecayNumber, mother_index, daughter_index;
	Bool_t VTXconv = false;

	Particle pi0("pi0");
	                              
	Double_t pt_min = 0.5;     											// used for all simulation
	Double_t pt_max = 15.;

	Int_t nevt = 100000000;
	Int_t nentry = 0;

	for (int i = 0; i < nevt; i++){   									// event loop

	  	if(i%10000000 == 0) cout << "event loop at " << i << endl;     	// event counter printed 

	  	MyEvent.ClearEvent();
	  	
	  	pi0.ResetP();                                                 	// reset
	 	pi0.GenerateP(pt_min,pt_max);                                 	// generate

	 	float trk_wt = 1;

	  	pi0.DecayFlat();                                              	// generate decay photon

    	ndecay = pi0.GetNumberOfDaughters();                          	// 2 for pi0->gg and 3 for pi0-e+e-g

    	DecayNumber = pi0.GetDecayNumber();								// 0 for pi0->gg and 1 for pi0-e+e-g

    	setTrack(MyTrack, 0, nentry, pi0.ID(), 0, pi0.Px(), pi0.Py(), pi0.Pz(), pi0.E(), pi0.M(), 0, 0, 0, DecayNumber, -999, trk_wt);
    	MyEvent.AddEntry(MyTrack); 
    	mother_index = nentry;
    	nentry++;

    	h_input->Fill(pi0.ID(), DecayNumber);

    	h_input_pt->Fill(pi0.Pt()); //pi0
    	h_input_eta->Fill(pi0.PseudoRapidity()); //pi0
    	h_input_phi->Fill(pi0.Phi()); //pi0
    	
    	for (Int_t j=0; j< ndecay; j++) {                             // loop over decay particles

    	  	temp = pi0.GetDecayDaughter(j);                             
    	  	int id = pi0.GetDaughterID(j);
    	  	trk_wt = pi0.GetDaughterWeight(j)*pi0.Weight();

    		setTrack(MyTrack, 1, nentry, id, 0, temp.Px(), temp.Py(), temp.Pz(), temp.E(), temp.M(), 0, 0, 0, 0, mother_index, trk_wt);
    		MyEvent.AddEntry(MyTrack);
    		nentry++;

    	}

    	MyEvent.SetNStable(ndecay);
    	T->Fill();

    	MyEvent.ClearEvent();
    	nentry = 0;
    }

    f->cd();

    h_input->Write();
    h_input_pt->Write();
    h_input_eta->Write();
    h_input_phi->Write();

    T->Write();
    f->Close();
}


void setTrack(WriteTrack& newTrack, int isFinal, int num, int id, int ist, float px, float py, float pz, float en, float mass, float xpos, float ypos, float zpos, int br, int parent_index, float weight){

	newTrack.SetFinal(isFinal);
	newTrack.SetNum(num);
	newTrack.SetID(id);
	newTrack.SetIst(ist);
	newTrack.SetPx(px);
	newTrack.SetPy(py);
	newTrack.SetPz(pz);
	newTrack.SetEnergy(en);
	newTrack.SetMass(mass);
	newTrack.SetXpos(xpos);
	newTrack.SetYpos(ypos);
	newTrack.SetZpos(zpos);
	newTrack.SetBranch(br);
	newTrack.SetParentIndex(parent_index);
	newTrack.SetWeight(weight);

}


void WriteROOT2Oscar(TString infile = "input.root", TString output = "oscar.txt", Int_t nstart = 0, Int_t nend = 1){

	TFile* input = new TFile(infile,"READ");
	if(!(input))
	{
	  cout << "no input file" << endl;
	  exit(1);
	}

	ofstream file(output);

	//Read in the TTrees 

	TTree* T = (TTree*)input->Get("T");
	TBranch* br = T->GetBranch("MyEvent");

	WriteEvent* event = 0;
	br->SetAddress(&event);

	for(int ievt = nstart; ievt < nend; ievt++){

		event->ClearEvent();
    	br->GetEntry(ievt);

		file << 0 << "\t" << event->GetNStable() << endl;

		for(int i = 0; i < event->GetNEntries(); i++){

			WriteTrack Track = event->GetWriteTrack(i);

			if(Track.GetFinal() == 1){

				file << Track.GetNum() << "\t" << Track.GetID() << "\t" << 0 << "\t" << Track.GetPx() << "\t" << Track.GetPy() << "\t" << Track.GetPz() << "\t" << Track.GetEnergy() << "\t" << Track.GetMass() << "\t" << Track.GetXpos() << "\t" << Track.GetYpos() << "\t" << Track.GetZpos() << "\t" << "0" << endl;

			}
		}

		file << "0" << "\t" << "0" << endl;
	}

}
