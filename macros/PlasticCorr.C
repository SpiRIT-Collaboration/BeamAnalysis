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


void PlasticCorr(){
    gROOT->ForceStyle();
    BeamRaw *raw = new BeamRaw();

    int first_run=2542;//2819 for 132Sn
    int last_run=2542;//3039 for 132Sn

    for (int ii=first_run; ii<=last_run; ii++){
        raw->fChain->AddFile(Form("data/run%i.ridf.root",ii),0,"raw");
    }

    raw->Init();

    Long64_t nentries = raw->fChain->GetEntriesFast();
    cout << " Entries = " << nentries << endl;

    TCanvas *cvs2 = new TCanvas("cvsF3", "", 0, 0, 800, 580);
    TCanvas *cvs3 = new TCanvas("cvsF7", "", 0, 0, 800, 580);
    TCanvas *cvs4 = new TCanvas("cvsF13_1", "", 0, 0, 800, 580);
    TCanvas *cvs5 = new TCanvas("cvsF13_2", "", 0, 0, 800, 580);


    TH2D *histF3Corr = new TH2D("histF3Corr","",100,0,100,100,-0.5,0.5);
    TH2D *histF7Corr = new TH2D("histF7Corr","",100,0,150,100,-0.5,0.5);
    TH2D *histF13_1Corr = new TH2D("histF13_1Corr","",300,-200,100,200,-1,1);
    TH2D *histF13_2Corr = new TH2D("histF13_2Corr","",300,-360,-60,200,-1,1);



    Double_t F3dt=0;
    Double_t F3q=0;
    Double_t F5dt=0;
    Double_t F5q=0;
    Double_t F7dt=0;
    Double_t F7q=0;
    Double_t F13_1dt=0;
    Double_t F13_1q=0;
    Double_t F13_2dt=0;
    Double_t F13_2q=0;

    for(int ientry=0; ientry < nentries; ientry++){
      if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
      raw->fChain->GetEvent(ientry);


      F3dt= raw->BigRIPSPlastic_fTLRaw[0]-raw->BigRIPSPlastic_fTRRaw[0];
      F3q= log(raw->BigRIPSPlastic_fQRRaw[0])-log(raw->BigRIPSPlastic_fQLRaw[0]);
      F7dt= raw->BigRIPSPlastic_fTLRaw[2]-raw->BigRIPSPlastic_fTRRaw[2];
      F7q= log(raw->BigRIPSPlastic_fQRRaw[2])-log(raw->BigRIPSPlastic_fQLRaw[2]);
      F13_1dt= raw->BigRIPSPlastic_fTLRaw[3]-raw->BigRIPSPlastic_fTRRaw[3];
      F13_1q= log(raw->BigRIPSPlastic_fQRRaw[3])-log(raw->BigRIPSPlastic_fQLRaw[3]);
      F13_2dt= raw->BigRIPSPlastic_fTLRaw[4]-raw->BigRIPSPlastic_fTRRaw[4];
      F13_2q= log(raw->BigRIPSPlastic_fQRRaw[4])-log(raw->BigRIPSPlastic_fQLRaw[4]);

      if(true){

	histF3Corr->Fill(F3dt,F3q);
	histF7Corr->Fill(F7dt,F7q);
	histF13_1Corr->Fill(F13_1dt,F13_1q);
	histF13_2Corr->Fill(F13_2dt,F13_2q);

	}

    }

    histF3Corr->SetTitle("F3 correlation");
    histF3Corr->GetXaxis()->SetTitle("t_{L}-t_{R}");
    histF3Corr->GetXaxis()->CenterTitle();
    histF3Corr->GetYaxis()->SetTitleOffset(1.5);
    histF3Corr->GetYaxis()->SetTitle("ln(q_{R}/q_{L})");
    histF3Corr->GetYaxis()->CenterTitle();

    histF7Corr->SetTitle("F7 correlation");
    histF7Corr->GetXaxis()->SetTitle("t_{L}-t_{R}");
    histF7Corr->GetXaxis()->CenterTitle();
    histF7Corr->GetYaxis()->SetTitleOffset(1.5);
    histF7Corr->GetYaxis()->SetTitle("ln(q_{R}/q_{L})");
    histF7Corr->GetYaxis()->CenterTitle();

    histF13_1Corr->SetTitle("F13.1 correlation");
    histF13_1Corr->GetXaxis()->SetTitle("t_{L}-t_{R}");
    histF13_1Corr->GetXaxis()->CenterTitle();
    histF13_1Corr->GetYaxis()->SetTitleOffset(1.5);
    histF13_1Corr->GetYaxis()->SetTitle("ln(q_{R}/q_{L})");
    histF13_1Corr->GetYaxis()->CenterTitle();

    histF13_2Corr->SetTitle("F13.2 correlation");
    histF13_2Corr->GetXaxis()->SetTitle("t_{L}-t_{R}");
    histF13_2Corr->GetXaxis()->CenterTitle();
    histF13_2Corr->GetYaxis()->SetTitleOffset(1.5);
    histF13_2Corr->GetYaxis()->SetTitle("ln(q_{R}/q_{L})");
    histF13_2Corr->GetYaxis()->CenterTitle();

    cvs2->cd();
    histF3Corr->Draw("colz");
    cvs2->SaveAs(Form("./figures/Plastic/F3corr.%i.%i.png",first_run,last_run));
    cvs2->SaveAs(Form("./figures/Plastic/F3corr.%i.%i.root",first_run,last_run));

    cvs3->cd();
    histF7Corr->Draw("colz");
    cvs3->SaveAs(Form("./figures/Plastic/F7corr.%i.%i.png",first_run,last_run));
    cvs3->SaveAs(Form("./figures/Plastic/F7corr.%i.%i.root",first_run,last_run));



    cvs4->cd();
    histF13_1Corr->Draw("colz");
    cvs4->SaveAs(Form("./figures/Plastic/F13_1.%i.%i.png",first_run,last_run));
    cvs4->SaveAs(Form("./figures/Plastic/F13_1.%i.%i.root",first_run,last_run));

    cvs5->cd();
    histF13_2Corr->Draw("colz");
    cvs5->SaveAs(Form("./figures/Plastic/F13_2.%i.%i.png",first_run,last_run));
    cvs5->SaveAs(Form("./figures/Plastic/F13_2.%i.%i.root",first_run,last_run));



}
