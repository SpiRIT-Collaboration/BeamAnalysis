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




void makeTSUM(){
    gROOT->ForceStyle();


    BeamRaw *raw = new BeamRaw();

    int first_run=2251;
    int last_run=2251;

    for (int ii=first_run; ii<=last_run; ii++){
        raw->fChain->AddFile(Form("data/run%i.ridf.root",ii),0,"raw");
    }
    raw->Init();



    Long64_t nentries = raw->fChain->GetEntriesFast();
    cout << " Entries = " << nentries << endl;


    int NumPPAC=18;
    Double_t TSUMX[NumPPAC];
    Double_t TSUMY[NumPPAC];
    Double_t F3X,F3Y,F5X,F5Y,F7X,F7Y;
    float PPAC_TSUMX[NumPPAC];
    float PPAC_TSUMY[NumPPAC];
    float PPAC_TSUMXsig[NumPPAC];
    float PPAC_TSUMYsig[NumPPAC];

    TH1D *histTSUMX[NumPPAC];
    TH1D *histTSUMY[NumPPAC];
    char nppac[150];
    for(int i=0;i<NumPPAC;i++){
        sprintf(nppac,"ppacx_%d",i);
        histTSUMX[i] = new TH1D(nppac,"ppacx",1000,0,200);
        sprintf(nppac,"ppacy_%d",i);
        histTSUMY[i] = new TH1D(nppac,"ppacy",1000,0,200);
    }


    for (int ip = 0; ip < NumPPAC; ip ++){
      PPAC_TSUMX[ip]=0;
      PPAC_TSUMY[ip]=0;
      PPAC_TSUMXsig[ip]=0;
      PPAC_TSUMYsig[ip]=0;
    }



      for(int ientry=0; ientry < nentries; ientry++){

        if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
        raw->fChain->GetEvent(ientry);

        F3X=raw->BigRIPSFocalPlane_X[3];
        F3Y=raw->BigRIPSFocalPlane_Y[3];
        F5X=raw->BigRIPSFocalPlane_X[5];
        F5Y=raw->BigRIPSFocalPlane_Y[5];
        F7X=raw->BigRIPSFocalPlane_X[7];
        F7Y=raw->BigRIPSFocalPlane_Y[7];


        for(int ippac = 0; ippac<NumPPAC; ippac++){
          //ResetPlots();
          TSUMX[ippac]=raw->BigRIPSPPAC_fTSumX[ippac];
          TSUMY[ippac]=raw->BigRIPSPPAC_fTSumY[ippac];
          //////////////////////////SET Focal Plane limits here////////////////////////////////
          //////an "OR" can be placed between true and focal plane conditions to allow all events into consideration///////
          ////// an "AND" should be blaced between the true and FP conditions to allow only events which fit focal plane considerations///////////
          if(true || (abs(F3X)<10 && abs(F3Y)<10 && abs(F5X)<10 && abs(F5Y)<10 && abs(F7X)<10 && abs(F7Y)<10 && F7X>0)){
              histTSUMX[ippac]->Fill(TSUMX[ippac]);
              histTSUMY[ippac]->Fill(TSUMY[ippac]);
          }
      }//loop thru ppacs
    }//

    TF1 *fit1 = new TF1("fit1","gaus",10,200);
    TF1 *fit2 = new TF1("fit2","gaus",10,200);
    TCanvas *c1=new TCanvas("c1","c1",800,800);

    for(int ippac = 0; ippac<NumPPAC; ippac++){

    c1->cd();

    histTSUMX[ippac]->Draw();
    histTSUMX[ippac]->Fit("fit1");
    c1->SaveAs(Form("./figures/TSUM/%i.TSUMX.%i.%i.png",ippac,first_run,last_run));
    c1->Clear();
    histTSUMY[ippac]->Draw();
    histTSUMY[ippac]->Fit("fit2");
    c1->SaveAs(Form("./figures/TSUM/%i.TSUMY.%i.%i.png",ippac,first_run,last_run));
    c1->Clear();

    PPAC_TSUMX[ippac]=fit1->GetParameter(1);
    PPAC_TSUMY[ippac]=fit2->GetParameter(1);
    PPAC_TSUMXsig[ippac]=fit1->GetParameter(2);
    PPAC_TSUMYsig[ippac]=fit2->GetParameter(2);
}

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
