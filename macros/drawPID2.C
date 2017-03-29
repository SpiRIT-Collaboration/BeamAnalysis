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


void drawPID2(){
    gROOT->ForceStyle();

    BeamBeam *beam = new BeamBeam();
    BeamRaw *raw = new BeamRaw();

    int first_run=2251;
    int last_run=first_run;//can set to be the first run, or extend to include a range of runs
    //add runs to the chain
    for (int ii=first_run; ii<=last_run; ii++){
        beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",ii),0,"beam");
        raw->fChain->AddFile(Form("data/run%i.ridf.root",ii),0,"raw");
    }

    beam->Init();
    raw->Init();

    Long64_t nentries = beam->fChainBeam->GetEntriesFast();
    cout << " Entries = " << nentries << endl;

    double goalAoQ=2.5;
    if(first_run>=2174 && last_run<=2509) goalAoQ=108./50.;
    if(first_run>=2520 && last_run<=2653) goalAoQ=112./50.;
    if(first_run>=3044 && last_run<=3184) goalAoQ=124./50.;
    if(first_run>=2819 && last_run<=3039) goalAoQ=132./50.;
    TCanvas *cvs = new TCanvas("cvspid", "", 0, 0, 800, 580);
    TCanvas *cvs2 = new TCanvas("cvspid2", "", 0, 0, 800, 580);
    auto *histBeamPID = new TH2D("histBeamPID", "", 500, goalAoQ-0.2, goalAoQ+0.2, 500, 40, 60);
    auto *histBeamPID2 = new TH2D("histBeamPID2", "", 1000, -100, -30, 1000, 0, 5000);
    auto *histIC=new TH1D("histIC","",1000,500,5000);
    Double_t F3PLA_time,F7PLA_time;//must use "Double_t" not double
    Double_t F7IC_energy;


    int i_z=50;
    int i_a=108;
    int goalA=(int)(goalAoQ*50);
    for(i_z = 49; i_z<=51;i_z++){
    for(i_a = goalA-4;i_a< goalA+3;i_a++){
      histIC->Reset();
      for(int ientry=0; ientry < nentries; ientry++){
      if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
      beam->fChainBeam->GetEvent(ientry);
      raw->fChain->GetEvent(ientry);

      F3PLA_time=raw->BigRIPSPlastic_fTime[0];
      F7PLA_time=raw->BigRIPSPlastic_fTime[2];
      F7IC_energy=0;
      /*
      for(int ii=0;ii<12;ii++){
        F7IC_energy=F7IC_energy+raw->BigRIPSIC_fEnergy[2][ii];
      }
      */
      F7IC_energy=raw->BigRIPSIC_fRawADCSqSum[2];
      double d_aoq=beam->aoq;
      double d_z=beam->z;

      if((d_z-(double)i_z)*(d_z-(double)i_z)/0.2+(d_aoq-(double)i_a/(double)i_z)*(d_aoq-(double)i_a/(double)i_z)/0.000025<1){//abs(d_aoq*d_z-(float)i_a)*abs(d_aoq*d_z-(float)i_a) + abs(d_z-(float)i_z)*abs(d_z-(float)i_z)<.25){
	    histBeamPID->Fill(beam->aoq,beam->z);
	    histBeamPID2->Fill(F7PLA_time-F3PLA_time,F7IC_energy);
	    histIC->Fill(F7IC_energy);
      }
      //cout << i_z << endl;
      //histIC->Fit("gaus");
      }
    cout << i_z << "," << i_a << endl;
    histIC->Fit("gaus");
    }
    }
    cvs->cd();
    histBeamPID->Draw("colz");
    cvs2->cd();
    //histBeamPID2->Draw("colz");
    //TF1 f1("f1","gaus",2000,3000);
    histIC->Draw();
    histIC->Fit("gaus");

}
