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

using namespace std;
Int_t *GetIsotope(Double_t myZ,Double_t myAoq){
  Double_t myMass=myAoq*myZ;
  Int_t static Arr[2]={0,0};//array of Z,mass
  int iZ1=(int)(myZ-1.);
  int iZ2=(int)(myZ+1.);
  int iM1=(int)(myMass-1.);
  int iM2=(int)(myMass+1.);
  Double_t ellipse=100;
  for(int iz=iZ1;iz<iZ2;iz++){
    for(int im=iM1;im<iM2;im++){
      ellipse = (myZ-(double)iz)*(myZ-(double)iz)/0.5/0.5+(myMass-(double)im)*(myMass-(double)im)/0.5/0.5;//the ellipse can be numerically changed here
      if(ellipse<1.){
        Arr[0]=iz;
        Arr[1]=im;
        //if(iz>50) cout << myZ << "," << iz << "," << myAoq << "," << im << endl;
        break;
      }
    }//mass loop
    if(ellipse<1.){
      break;
    }
  }//charge loop
  return Arr;
}

void SimpleOut2(int runNo=2843){
    gROOT->ForceStyle();
    cout << runNo << endl;
    //beambeam and beamraw objects to store information from first ridf root file
    BeamBeam *beam = new BeamBeam();
    BeamRaw *raw = new BeamRaw();
    //tree to store information about PPAC hits
    TChain *ppactree = new TChain("tr");
    TChain *bdcinfo = new TChain("bdcinfo");
    TChain *MAGframe = new TChain("MAGframe");



    //Output file and Trees to write out
    TFile *fout = new TFile(Form("./output/beam/beam_run.%i.root",runNo),"recreate");
    //a tree to hold charge, aoq, beta, brho, integer mass and integer charge, as well as "isGood" flag
    auto TBeam = new TTree("TBeam","TBeam");
    //a tree to hold BDC position and momentum information. Linear and magnetic projection
    auto TBDC = new TTree("TBDC","TBDC");
    //a tree to hold focal plane information
    auto TFocalPlane = new TTree("TFocalPlane","TFocalPlane");

    //auto beam_out = new TTree("beam_out","beam_out");

    //add root files to the trees
    beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",runNo),0,"beam");
    raw->fChain->AddFile(Form("data/run%i.ridf.root",runNo),0,"raw");
    ppactree->AddFile(Form("output/ppac/ppac_hit.%i.root",runNo));
    bdcinfo->AddFile(Form("output/BDC/BDCout.%i.ridf.root",runNo));
    MAGframe->AddFile(Form("output/BDC/BDCout.%i.ridf.root",runNo));

//need to implement addfriends

    beam->Init();
    raw->Init();


    //TCanvas *cvs2 = new TCanvas("cvspid", "", 0, 0, 800, 580);
    auto *histBeamPID = new TH2D("histBeamPID", "", 1000, 2., 2.7, 2000, 35, 55);

    Long64_t nentries = beam->fChainBeam->GetEntriesFast();
    cout << " Entries = " << nentries << endl;

    //load plastic cut files
    int plastic_cut_no=2819;
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


    //define variables to be used within this macro
    double F3dt,F3q,F7dt,F7q,F13_1dt,F13_1q,F13_2dt,F13_2q;
    bool F3in,F7in,F13_1in,F13_2in;

    //define variables and branches to be written out
    int neve = 0;
    TBeam -> Branch("neve",&neve);
    TBDC -> Branch("neve",&neve);

    //TBeam branches
    Double_t aoq=-9999; Double_t z=-9999; Double_t tof=-9999; Double_t beta78=-9999; Double_t brho78=-9999;
    Int_t intZ=-9999;
    Int_t intMass=-9999;
    Bool_t isGood=false;
    TBeam->Branch("aoq",&aoq,"aoq/D");
    TBeam->Branch("z",&z,"z/D");
    TBeam->Branch("tof",&tof,"tof/D");
    TBeam->Branch("brho78",&brho78,"brho78/D");
    TBeam->Branch("intZ",&intZ,"intZ/I");
    TBeam->Branch("intMass",&intMass,"intMass/I");
    TBeam->Branch("isGood",&isGood,"isGood/B");

    //TBDC branches
    Double_t brho,beta;
    TBDC->Branch("beta",&beta,"beta/D");//this beta is energy corrected
    TBDC->Branch("brho",&brho,"brho/D");//this brho is energy corrected
    TBeam->Branch("beta",&beta,"beta/D");//this beta is energy corrected
    TBeam->Branch("brho",&brho,"brho/D");//this brho is energy corrected
    Double_t bdc1x,bdc1y,bdc2x,bdc2y;
    TBDC->Branch("bdc1x",&bdc1x,"bdc1x/D");
    TBDC->Branch("bdc1y",&bdc1y,"bdc1y/D");
    TBDC->Branch("bdc2x",&bdc2x,"bdc2x/D");
    TBDC->Branch("bdc2y",&bdc2y,"bdc2y/D");

    Double_t projPX=-9999; Double_t projPY=-9999; Double_t projPZ=-9999; Double_t projX=-9999; Double_t projY=-9999; Double_t projZ=-9999; Double_t projA=-9999; Double_t projB;
    TBDC->Branch("projX",&projX,"projX/D");
    TBDC->Branch("projY",&projY,"projY/D");
    TBDC->Branch("projZ",&projY,"projZ/D");
    TBDC->Branch("projA",&projA,"projA/D");
    TBDC->Branch("projB",&projB,"projB/D");
    TBDC->Branch("projPX",&projPX,"projPX/D");
    TBDC->Branch("projPY",&projPY,"projPY/D");
    TBDC->Branch("projPZ",&projPZ,"projPZ/D");

    //The Focalplane:
    Double_t F3X=-9999; Double_t F3A=-9999; Double_t F3Y=-9999; Double_t F3B=-9999;
    Double_t F5X=-9999; Double_t F5A=-9999; Double_t F5Y=-9999; Double_t F5B=-9999;
    Double_t F7X=-9999; Double_t F7A=-9999; Double_t F7Y=-9999; Double_t F7B=-9999;

    TFocalPlane->Branch("F3X",&F3X,"F3X/D");
    TFocalPlane->Branch("F3A",&F3A,"F3A/D");
    TFocalPlane->Branch("F3Y",&F3Y,"F3Y/D");
    TFocalPlane->Branch("F3B",&F3B,"F3B/D");

    TFocalPlane->Branch("F5X",&F5X,"F5X/D");
    TFocalPlane->Branch("F5A",&F5A,"F5A/D");
    TFocalPlane->Branch("F5Y",&F5Y,"F5Y/D");
    TFocalPlane->Branch("F5B",&F5B,"F5B/D");

    TFocalPlane->Branch("F7X",&F7X,"F7X/D");
    TFocalPlane->Branch("F7A",&F7A,"F7A/D");
    TFocalPlane->Branch("F7Y",&F7Y,"F7Y/D");
    TFocalPlane->Branch("F7B",&F7B,"F7B/D");



    //loop through entries
    for(int ientry=0; ientry < nentries; ientry++){
      //output progress every 100,000 events
      if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
      //get current entry
      beam->fChainBeam->GetEvent(ientry);
      raw->fChain->GetEvent(ientry);
      ppactree->GetEntry(ientry);
      bdcinfo->GetEntry(ientry);
      MAGframe->GetEntry(ientry);

      //initialize output variables
      aoq=beam->aoq; z=beam->z; tof=beam->tof; beta78=beam->beta;//brho78=beam->brho;

      //projPX=MAGframe->MAGpx;projPY=MAGframe->MAGpy;projPZ=MAGframe->MAGpz;projX=MAGframe->MAGx;projY=MAGframe->MAGy;projZ=MAGframe->MAGz;
      //bdc1x=bdcinfo->bdc1trx;bdc1y=bdcinfo->bdc1try;bdc2x=bdcinfo->bdc2trx;bdc2y=bdcinfo->bdc2try;
      //beta=MAGframe->beta;
      //brho=brho78*beta/beta78;

      F3X=beam->F3X; F3A=beam->F3A; F3Y=beam->F3Y; F3B=beam->F3B;
      F5X=beam->F5X; F5A=beam->F5A; F5Y=beam->F5Y; F5B=beam->F5B;
      F7X=beam->F7X; F7A=beam->F7A; F7Y=beam->F7Y; F7B=beam->F7B;

      intZ=-9999;intMass=-9999;
      isGood=false;
      // initialize booleans about cuts
      F3in = false;
      F7in = false;
      F13_1in = false;
      F13_2in = false;


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
      //set conditions to require for "good" events
      if( F3in && F7in && F13_1in && F13_2in && aoq>0 && aoq<3 && z>0 && z<100){
        isGood=true;
      }
      //if good,
      if(isGood){
      	histBeamPID->Fill(aoq,z);
        intZ=GetIsotope(z,aoq)[0];
        intMass=GetIsotope(z,aoq)[1];
      }
      //fill trees
      TBeam -> Fill();
      TBDC -> Fill();
      TFocalPlane -> Fill();
      neve ++;

    }


//write out file
    fout->cd();
    TBeam->Write();
    TBDC->Write();
    TFocalPlane -> Write();
    fout->Write();
    fout->Close();
    std::cout << "neve: " << neve << ", nentries: " << nentries << std::endl;


}
