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
#include "Riostream.h"

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "BeamRaw.h"
#include "BeamBeam.h"

using namespace std;

const char *cutfname = NULL;
bool filecalled = false;

TH1D *histTSUMX=new TH1D("histTSUMX","histTSUMX",1000,0,200);
TH1D *histTSUMY=new TH1D("histTSUMY","histTSUMY",1000,0,200);



void ResetPlots(){
  histTSUMX->Reset();
  histTSUMY->Reset();
}


void makeTSUM(){
    gROOT->ForceStyle();

    BeamBeam *beam = new BeamBeam();
    BeamRaw *raw = new BeamRaw();

    int first_run=2819;
    int last_run=3184;

    for (int ii=first_run; ii<=last_run; ii++){
        beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",ii),0,"beam");
        raw->fChain->AddFile(Form("data/run%i.ridf.root",ii),0,"raw");
    }
    beam->Init();
    raw->Init();



    Long64_t nentries = beam->fChainBeam->GetEntriesFast();
    cout << " Entries = " << nentries << endl;


    int NumPPAC=18;
    Double_t TSUMX;
    Double_t TSUMY;
    float PPAC_TSUMX[NumPPAC];
    float PPAC_TSUMY[NumPPAC];
    float PPAC_TSUMXsig[NumPPAC];
    float PPAC_TSUMYsig[NumPPAC];

    for (int ip = 0; ip < NumPPAC; ip ++){
      PPAC_TSUMX[ip]=0;
      PPAC_TSUMY[ip]=0;
      PPAC_TSUMXsig[ip]=0;
      PPAC_TSUMYsig[ip]=0;
    }


    for(int ippac = 12; ippac<NumPPAC; ippac++){
      ResetPlots();

      for(int ientry=0; ientry < nentries; ientry++){

        if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
        beam->fChainBeam->GetEvent(ientry);
        raw->fChain->GetEvent(ientry);
        TSUMX=raw->BigRIPSPPAC_fTSumX[ippac];
        TSUMY=raw->BigRIPSPPAC_fTSumY[ippac];
        histTSUMX->Fill(TSUMX);
        histTSUMY->Fill(TSUMY);

      }//loop thru events

      TF1 *fit1 = new TF1("fit1","gaus",10,200);
      TF1 *fit2 = new TF1("fit2","gaus",10,200);
      TCanvas *c1=new TCanvas("c1","c1",800,800);
      c1->cd();

      histTSUMX->Draw();
      histTSUMX->Fit("fit1");
      c1->SaveAs(Form("./figures/TSUM/%i.TSUMX.%i.%i.png",ippac,first_run,last_run));
      c1->Clear();
      histTSUMY->Draw();
      histTSUMY->Fit("fit2");
      c1->SaveAs(Form("./figures/TSUM/%i.TSUMY.%i.%i.png",ippac,first_run,last_run));
      c1->Clear();

      PPAC_TSUMX[ippac]=fit1->GetParameter(1);
      PPAC_TSUMY[ippac]=fit2->GetParameter(1);
      PPAC_TSUMXsig[ippac]=fit1->GetParameter(2);
      PPAC_TSUMYsig[ippac]=fit2->GetParameter(2);



    }//loop thru ppacs
    std::ofstream csvout (Form("./output/ppac_tsumoutput.%i.%i.csv",first_run,last_run),std::ofstream::out);
    csvout << "PPAC,TSUMX mean,TSUMX sigma,TSUMY mean,TSUMY sigma,txsum_min,txsum_max,tysum_min,tysum_max"<<endl;
//csvout << "PPAC,Xslope,Xoff,Yslope,Yoff,XvydiffSlope,Xvydiffoff"<<endl;
    for(int ippac = 0; ippac<NumPPAC; ippac++){
      csvout << ippac << "," << PPAC_TSUMX[ippac]<< "," << PPAC_TSUMXsig[ippac]<< "," << PPAC_TSUMY[ippac]<< "," << PPAC_TSUMYsig[ippac]
      << "," << PPAC_TSUMX[ippac]-3*PPAC_TSUMXsig[ippac]<< "," << PPAC_TSUMX[ippac]+3*PPAC_TSUMXsig[ippac]
      << "," << PPAC_TSUMY[ippac]-3*PPAC_TSUMYsig[ippac]<< "," << PPAC_TSUMY[ippac]+3*PPAC_TSUMYsig[ippac]
      << endl;
    }
    csvout.close();
}
