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


void SimpleOut(int runNo){
    gROOT->ForceStyle();
    cout << runNo << endl;
    //beambeam and beamraw objects to store information from first ridf root file
    BeamBeam *beam = new BeamBeam();
    BeamRaw *raw = new BeamRaw();
    //tree to store information about PPAC hits
    //TChain *tree = new TChain("tr");

    //Output file and Trees to write out
    TFile *fout = new TFile(Form("./output/SimpleData/beam_run.%i.root",runNo),"recreate");
    auto out = new TTree("out","out");
    auto BDC = new TTree("BDC","BDC");
    //auto beam_out = new TTree("beam_out","beam_out");


    beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",runNo),0,"beam");
    raw->fChain->AddFile(Form("data/run%i.ridf.root",runNo),0,"raw");

    beam->Init();
    raw->Init();

    //TCanvas *cvs2 = new TCanvas("cvspid", "", 0, 0, 800, 580);
    auto *histBeamPID = new TH2D("histBeamPID", "", 1000, 2.6, 2.7, 5000,45, 55);

    Long64_t nentries = beam->fChainBeam->GetEntriesFast();
    cout << " Entries = " << nentries << endl;

    //load cut files
    TFile* cFileoutF3cut = new TFile("./cut/plastic_cuts/F3cut.2819.root");
    TCutG *F3pl;
    cFileoutF3cut->GetObject("CUTG",F3pl);
    cFileoutF3cut->Close();

    TFile* cFileoutF7cut = new TFile("./cut/plastic_cuts/F7cut.2819.root");
    TCutG *F7pl;
    cFileoutF7cut->GetObject("CUTG",F7pl);
    cFileoutF7cut->Close();

    TFile* cFileoutF13cut1 = new TFile("./cut/plastic_cuts/F13_1cut.2819.root");
    TCutG *F13pl1;
    cFileoutF13cut1->GetObject("CUTG",F13pl1);
    cFileoutF13cut1->Close();

    TFile* cFileoutF13cut2 = new TFile("./cut/plastic_cuts/F13_2cut.2819.root");
    TCutG *F13pl2;
    cFileoutF13cut2->GetObject("CUTG",F13pl2);
    cFileoutF13cut2->Close();


    //define variables
    double F3dt,F3q,F7dt,F7q,F13_1dt,F13_1q,F13_2dt,F13_2q;
    bool F3in,F7in,F13_1in,F13_2in;
    //define variables and branches to be written out
    int neve = 0;
    out -> Branch("neve",&neve);
    Double_t AoQ = -9999;
    out -> Branch("aoq",&AoQ);
    Double_t charge = -9999;
    out -> Branch("z",&charge);
    out -> Branch("neve",&neve);

    for(int ientry=0; ientry < nentries; ientry++){
      //initialize output variables
      AoQ = -9999;
      charge = -9999;

      // initialize booleans about cuts
      F3in = false;
      F7in = false;
      F13_1in = false;
      F13_2in = false;


      if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
      beam->fChainBeam->GetEvent(ientry);
      raw->fChain->GetEvent(ientry);


      F3dt= raw->BigRIPSPlastic_fTLRaw[0]-raw->BigRIPSPlastic_fTRRaw[0];
      F3q= log(raw->BigRIPSPlastic_fQRRaw[0])-log(raw->BigRIPSPlastic_fQLRaw[0]);
      F7dt= raw->BigRIPSPlastic_fTLRaw[2]-raw->BigRIPSPlastic_fTRRaw[2];
      F7q= log(raw->BigRIPSPlastic_fQRRaw[2])-log(raw->BigRIPSPlastic_fQLRaw[2]);
      F13_1dt= raw->BigRIPSPlastic_fTLRaw[3]-raw->BigRIPSPlastic_fTRRaw[3];
      F13_1q= log(raw->BigRIPSPlastic_fQRRaw[3])-log(raw->BigRIPSPlastic_fQLRaw[3]);
      F13_2dt= raw->BigRIPSPlastic_fTLRaw[4]-raw->BigRIPSPlastic_fTRRaw[4];
      F13_2q= log(raw->BigRIPSPlastic_fQRRaw[4])-log(raw->BigRIPSPlastic_fQLRaw[4]);

      if(F3pl->IsInside(F3dt,F3q)) F3in = true;
      if(F7pl->IsInside(F7dt,F7q)) F7in = true;
      if(F13pl1->IsInside(F13_1dt,F13_1q)) F13_1in = true;
      if(F13pl2->IsInside(F13_2dt,F13_2q)) F13_2in = true;

      if( true){
      	histBeamPID->Fill(beam->aoq,beam->z);
      	AoQ = beam -> aoq;
      	charge = beam -> z;
      }
      //fill trees
      out -> Fill();
      BDC -> Fill();
      neve ++;

    }
//write out file
    fout->cd();
    out->Write();
    BDC->Write();
    fout->Write();
    fout->Close();
    std::cout << "neve: " << neve << ", nentries: " << nentries << std::endl;

}
