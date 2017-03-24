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

const char *cutfname = NULL;
bool filecalled = false;


void ICplaCorr(){
    gROOT->ForceStyle();
    BeamRaw *raw = new BeamRaw();

    int first_run=2819;//2819 for 132Sn
    int last_run=3039;//3039 for 132Sn

    for (int ii=first_run; ii<=last_run; ii++){
        raw->fChain->AddFile(Form("data/run%i.ridf.root",ii),0,"raw");
    }

    raw->Init();

    Long64_t nentries = raw->fChain->GetEntriesFast();
    cout << " Entries = " << nentries << endl;

    TCanvas *cvs2 = new TCanvas("cvs2", "", 0, 0, 800, 580);
    TCanvas *cvs3 = new TCanvas("cvs3", "", 0, 0, 800, 580);



    TH2D *histF3ICCorr = new TH2D("histF3ICCorr","",600,400,1000,1000,0,6000);
    TH2D *histF7ICCorr = new TH2D("histF7ICCorr","",600,400,1000,1000,0,6000);





    Double_t F3q=0;
    Double_t F7q=0;
    Double_t ICq=0;


    for(int ientry=0; ientry < nentries; ientry++){
      if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
      raw->fChain->GetEvent(ientry);

      ICq=0;
      for(int ii=0;ii<12;ii++){
        ICq=ICq+raw->BigRIPSIC_fEnergy[2][ii];
      }

      F3q= raw->BigRIPSPlastic_fQRRaw[0]+raw->BigRIPSPlastic_fQLRaw[0];
      F7q= raw->BigRIPSPlastic_fQRRaw[2]+raw->BigRIPSPlastic_fQLRaw[2];

      if(true){

	histF3ICCorr->Fill(F3q,ICq);
	histF7ICCorr->Fill(F7q,ICq);


	}

    }

    histF3ICCorr->SetTitle("F3 / Ion Chamber correlation");
    histF3ICCorr->GetXaxis()->SetTitle("F3 charge");
    histF3ICCorr->GetXaxis()->CenterTitle();
    histF3ICCorr->GetYaxis()->SetTitleOffset(1.5);
    histF3ICCorr->GetYaxis()->SetTitle("IC charge");
    histF3ICCorr->GetYaxis()->CenterTitle();

    histF7ICCorr->SetTitle("F7 / Ion Chamber correlation");
    histF7ICCorr->GetXaxis()->SetTitle("F7 charge");
    histF7ICCorr->GetXaxis()->CenterTitle();
    histF7ICCorr->GetYaxis()->SetTitleOffset(1.5);
    histF7ICCorr->GetYaxis()->SetTitle("IC charge");
    histF7ICCorr->GetYaxis()->CenterTitle();

    cvs2->cd();
    histF3ICCorr->Draw("colz");
    cvs2->SaveAs(Form("./figures/Plastic/F3ICCorr.%i.%i.png",first_run,last_run));
    cvs2->SaveAs(Form("./figures/Plastic/F3ICCorr.%i.%i.root",first_run,last_run));

    cvs3->cd();
    histF7ICCorr->Draw("colz");
    cvs3->SaveAs(Form("./figures/Plastic/F7ICCorr.%i.%i.png",first_run,last_run));
    cvs3->SaveAs(Form("./figures/Plastic/F7ICCorr.%i.%i.root",first_run,last_run));






}
