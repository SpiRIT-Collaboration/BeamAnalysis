#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TChain.h"
#include "TCutG.h"
#include <TFile.h>
#include "TCanvas.h"
#include "TH2D.h"
#include "TH2.h"
#include "cmath"

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "BeamRaw.h"
#include "BeamBeam.h"

using namespace std;


void drawPID(){
    gROOT->ForceStyle();

    BeamBeam *beam = new BeamBeam();
    BeamRaw *raw = new BeamRaw();

    int first_run=2225;
    int last_run=2255;
    //add runs to the chain
    for (int ii=first_run; ii<=last_run; ii++){
        beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",ii),0,"beam");
        raw->fChain->AddFile(Form("data/run%i.ridf.root",ii),0,"raw");
    }

    beam->Init();
    raw->Init();

    Long64_t nentries = beam->fChainBeam->GetEntriesFast();
    cout << " Entries = " << nentries << endl;

    TCanvas *cvs2 = new TCanvas("cvspid2", "", 0, 0, 800, 580);
    auto *histBeamPID = new TH2D("histBeamPID", "", 500, 2.1, 2.3, 500, 40, 60);

    for(int ientry=0; ientry < nentries; ientry++){
      if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
      beam->fChainBeam->GetEvent(ientry);
      raw->fChain->GetEvent(ientry);
      if(true){
	       histBeamPID->Fill(beam->aoq,beam->z);
	    }
    }
    cvs2->cd();
    histBeamPID->Draw("colz");


}
