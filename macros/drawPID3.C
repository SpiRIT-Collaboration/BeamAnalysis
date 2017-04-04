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

using namespace std;


void drawPID3(){
    gROOT->ForceStyle();
 
    auto TBeam = new TChain("TBeam");
    
    int first_run=2251;
    int last_run=2251;//can set to be the first run, or extend to include a range of runs
    //add runs to the chain
    for (int ii=first_run; ii<=last_run; ii++){
      TBeam->AddFile(Form("./output/beam/beam_run%i.ridf.root",ii));
    }



    Long64_t nentries = TBeam->GetEntries();
    cout << " Entries = " << nentries << endl;

    double goalAoQ=2.5;
    if(first_run>=2174 && last_run<=2509) goalAoQ=108./50.;
    if(first_run>=2520 && last_run<=2653) goalAoQ=112./50.;
    if(first_run>=3044 && last_run<=3184) goalAoQ=124./50.;
    if(first_run>=2819 && last_run<=3039) goalAoQ=132./50.;
    TCanvas *cvs = new TCanvas("cvspid", "", 0, 0, 800, 580);

    auto *histBeamPID = new TH2D("histBeamPID", "", 500, goalAoQ-0.1, goalAoQ+0.1, 500, 45, 55);
    /*
      for(int ientry=0; ientry < nentries; ientry++){
      if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
      TBeam->GetEntry(ientry);
      //do stuff event by event
      }
*/
    cvs->cd();
    TBeam->Draw("z:aoq >> histBeamPID","isGood","colz");


}
