#include <fstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TSystem.h"
#include "TNtuple.h"
#include "TArtEventStore.hh"
#include "TArtStoreManager.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtSAMURAIParameters.hh"
#include "TArtCalibBDC1Hit.hh"
#include "TArtCalibBDC1Track.hh"
#include "TArtCalibBDC2Hit.hh"
#include "TArtCalibBDC2Track.hh"
#include "TArtDCTrack.hh"
#include "TArtCalibIC.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtFocalPlane.hh"
#include "TArtEventInfo.hh"
#include "TArtPlastic.hh"
#include "TArtPPAC.hh"
#include "TArtCalibPID.hh"
#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"
#include "TArtPPACPara.hh"
#include "TArtBeam.hh"
#include "TArtTOF.hh"
#include "TArtIC.hh"
#include "TArtRIPS.hh"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraph.h"
#include "TClonesArray.h"
#include "TCutG.h"
#include "TBeamEnergy.h"
#include "TBDCProjection.h"


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
        break;
      }
    }//mass loop
    if(ellipse<1.){
      break;
    }
  }//charge loop
  return Arr;
}

int main(int argc, char *argv[]) {
    float customTOF;
    if(argc < 2){
        cerr << "Missing RIDF file argument" << endl;
    }
    if(argc < 3){
        cerr << "Missing TOF argument, using 280 ns" << endl;
        customTOF=280.;
    }
    else{
        char* charTOF = argv[2];
        customTOF = atof(charTOF);
    }
    char* charRun = argv[1];
    int runNo;
    runNo = atoi(charRun);

    cerr << "Run Number =" << runNo << endl;


    Int_t maxEvents = 50000000;
    auto bigripsParameters = TArtBigRIPSParameters::Instance();
    bigripsParameters -> LoadParameter((Char_t *) "db/BigRIPSPPAC.xml");
    bigripsParameters -> LoadParameter((Char_t *) "db/BigRIPSPlastic.xml");
    bigripsParameters -> LoadParameter((Char_t *) "db/BigRIPSIC.xml");
    bigripsParameters -> LoadParameter((Char_t *) "db/FocalPlane.xml");

    auto samurai_prm = TArtSAMURAIParameters::Instance();
    samurai_prm -> LoadParameter((Char_t *) "db/SAMURAIBDC1.xml");
    samurai_prm -> LoadParameter((Char_t *) "db/SAMURAIBDC2.xml");

    auto eventStore = new TArtEventStore();
  TString ridfFile = Form("ridf/SMDAQ%d.ridf", runNo);

  eventStore -> Open(ridfFile.Data());
//  eventStore -> Open();// for online analysis

  auto calibbdc1hit = new TArtCalibBDC1Hit;
  auto calibbdc1tr = new TArtCalibBDC1Track;
  auto calibbdc2hit = new TArtCalibBDC2Hit;
  auto calibbdc2tr = new TArtCalibBDC2Track;
  calibbdc1tr -> SetTDCWindow(700,1500);
  calibbdc2tr -> SetTDCWindow(700,1500);

  auto storeMan = TArtStoreManager::Instance();
  TArtRawEventObject *rawevent = eventStore->GetRawEventObject();

  auto calibPID  = new TArtCalibPID();
  auto calibPPAC = calibPID -> GetCalibPPAC();
  auto calibPLA  = calibPID -> GetCalibPlastic();
  auto calibFP   = calibPID -> GetCalibFocalPlane();
  auto calibIC   = calibPID -> GetCalibIC();

  auto recoPID = new TArtRecoPID();
  // Construct BigRIPS segment
  // Construct transport matrix from the file
  // Set the BigRIPS parameters: up & downstream focal plane, center brho
  // Load the focal plane parameters
  auto rips35 = recoPID -> DefineNewRIPS(3, 5, (Char_t *) "matrix/mat1.mat", (Char_t *) "D3"); // F3 - F5
  auto rips57 = recoPID -> DefineNewRIPS(5, 7, (Char_t *) "matrix/mat2.mat", (Char_t *) "D5"); // F5 - F7

  // Construct ToF segment
  // Calculate the flight length from upstream to downstream focal plane.
  // Calculate the flight length of up to middle and middle to downstream focal plane (last argument is middle focal plane)
  // Load the focal plane parameters
  //////////////////////
  TArtTOF *tof37 = new TArtTOF();
  tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", customTOF, 5);
  //TArtTOF *tof57 = new TArtTOF();
  //tof57 = recoPID -> DefineNewTOF((Char_t *) "F5pl", (Char_t *) "F7pl", customTOF, 6);

  TArtBeam *beam_br37 = recoPID -> DefineNewBeam(rips35, rips57, tof37, (Char_t *) "F7IC");
  //TArtBeam *beam_br57 = recoPID -> DefineNewBeam(rips57, tof37, (Char_t *) "F7IC");
///////////////////////
  //auto tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", customTOF); // F3 - F7



  double goalAoQ;
  if(runNo>=2174 && runNo<=2509) goalAoQ=108./50.;
  if(runNo>=2520 && runNo<=2653) goalAoQ=112./50.;
  if(runNo>=3044 && runNo<=3184) goalAoQ=124./50.;
  if(runNo>=2819 && runNo<=3039) goalAoQ=132./50.;

  //histograms to be written with file
  auto *histBeamPID = new TH2D("histBeamPID", "", 2000, goalAoQ-0.2, goalAoQ+0.2, 2000, 35, 55);
  auto *histAllBeamPID = new TH2D("histAllBeamPID", "", 1000, goalAoQ-0.2, goalAoQ+0.2, 2000, 35, 55);
  auto *histBeamEnergy= new TH1D("histBeamEnergy", "", 1000, 100000, 200000 );
  auto *histBeamEnergy78= new TH1D("histBeamEnergy78", "", 1000, 100000, 200000 );
  auto *histBeta= new TH1D("histBeta", "", 1000, 0.6, 0.7);
  auto *histBeta37= new TH1D("histBeta37", "", 1000, 0.6, 0.7);
  auto *histProjection = new TH2D("histProjection", "", 1000, -50, 50, 2000, -50, 50);
  auto *histBrho= new TH1D("histBrho", "", 1000, 6, 8);
  auto *histBrho2= new TH1D("histBrho2", "", 1000, 6, 8);


  auto hF3corr = new TH2F("F3corr","F3 corr",2000,0,100,2000,0,0.5);
  auto hF7corr = new TH2F("F7corr","F7 corr",2000,0,100,2000,-0.3,0.2);
  auto hF131corr = new TH2F("F131corr","F13 1corr",2000,-200,100,2000,-1,1);
  auto hF132corr = new TH2F("F132corr","F13 2 corr",2000,-350,-100,2000,-1,1);
  //outbut file
  TFile *outfile = new TFile(Form("output/beam/beam_run%d.ridf.root", runNo), "recreate");
  //output trees


  auto cut_tree = new TTree("cut_tree", "cuts made");
  auto TBeam = new TTree("TBeam", "Beam information");
  auto TBDC = new TTree("TBDC", "BDC information");
  auto TFocalPlane = new TTree("TFocalPlane", "Focal Plane Information");

  //output branches
  bool fill_cuts=true;
  //optional cut branch
  Bool_t F3PLA_in;
  Bool_t F7PLA_in;
  Bool_t F13_1PLA_in;
  Bool_t F13_2PLA_in;
  Bool_t totalPPAChits_lte96;
  Bool_t F7A_PPAChits_lte8;
  Bool_t F7B_PPAChits_lte8;
  if(fill_cuts){
      cut_tree->Branch("F3PLA_in",&F3PLA_in,"F3PLA_in/B");cut_tree->Branch("F7PLA_in",&F7PLA_in,"F7PLA_in/B");
      cut_tree->Branch("F13_1PLA_in",&F13_1PLA_in,"F13_1PLA_in/B");cut_tree->Branch("F13_2PLA_in",&F13_2PLA_in,"F13_2PLA_in/B");
      cut_tree->Branch("totalPPAChits_lte96",&totalPPAChits_lte96,"totalPPAChits_lte96/B");
      cut_tree->Branch("F7A_PPAChits_lte8",&F7A_PPAChits_lte8,"F7A_PPAChits_lte8/B");
      cut_tree->Branch("F7B_PPAChits_lte8",&F7B_PPAChits_lte8,"F7B_PPAChits_lte8/B");
  }



  //optional raw branch
  Double_t F3PPAC1A_X,F3PPAC1A_Y,F3PPAC1B_X,F3PPAC1B_Y,F3PPAC2A_X, F3PPAC2A_Y,F3PPAC2B_X,F3PPAC2B_Y;
  Double_t F5PPAC1A_X,F5PPAC1A_Y,F5PPAC1B_X,F5PPAC1B_Y,F5PPAC2A_X,F5PPAC2A_Y, F5PPAC2B_X, F5PPAC2B_Y;
  Double_t F7PPAC1A_X,F7PPAC1A_Y,F7PPAC1B_X,F7PPAC1B_Y, F7PPAC2A_X,F7PPAC2A_Y, F7PPAC2B_X,F7PPAC2B_Y;
  Int_t ICQ[8];
  Int_t F3PLA_QL,F3PLA_QR;
  Double_t F3PLA_TL,F3PLA_TR,F3PLA_DT,F3PLA_Q_TEST;
  Int_t F7PLA_QL,F7PLA_QR;
  Double_t F7PLA_TL,F7PLA_TR,F7PLA_DT,F7PLA_Q_TEST;
  Int_t F13_1PLA_QL,F13_1PLA_QR;
  Double_t F13_1PLA_TL,F13_1PLA_TR,F13_1PLA_DT,F13_1PLA_Q_TEST;
  Int_t F13_2PLA_QL,F13_2PLA_QR;
  Double_t F13_2PLA_TL,F13_2PLA_TR,F13_2PLA_DT,F13_2PLA_Q_TEST;
/*
    auto raw = new TTree("raw", "rawdata");
    raw -> Branch("F3PPAC1A_X", &F3PPAC1A_X, "F3PPAC1A_X/D");raw -> Branch("F3PPAC1A_Y", &F3PPAC1A_Y, "F3PPAC1A_Y/D");
    raw -> Branch("F3PPAC1B_X", &F3PPAC1B_X, "F3PPAC1B_X/D");raw -> Branch("F3PPAC1B_Y", &F3PPAC1B_Y, "F3PPAC1B_Y/D");
    raw -> Branch("F3PPAC2A_X", &F3PPAC2A_X, "F3PPAC2A_X/D");raw -> Branch("F3PPAC2A_Y", &F3PPAC2A_Y, "F3PPAC2A_Y/D");
    raw -> Branch("F3PPAC2B_X", &F3PPAC2B_X, "F3PPAC2B_X/D");raw -> Branch("F3PPAC2B_Y", &F3PPAC2B_Y, "F3PPAC2B_Y/D");

    raw -> Branch("F5PPAC1A_X", &F5PPAC1A_X, "F5PPAC1A_X/D");raw -> Branch("F5PPAC1A_Y", &F5PPAC1A_Y, "F5PPAC1A_Y/D");
    raw -> Branch("F5PPAC1B_X", &F5PPAC1B_X, "F5PPAC1B_X/D");raw -> Branch("F5PPAC1B_Y", &F5PPAC1B_Y, "F5PPAC1B_Y/D");
    raw -> Branch("F5PPAC2A_X", &F5PPAC2A_X, "F5PPAC2A_X/D");raw -> Branch("F5PPAC2A_Y", &F5PPAC2A_Y, "F5PPAC2A_Y/D");
    raw -> Branch("F5PPAC2B_X", &F5PPAC2B_X, "F5PPAC2B_X/D");raw -> Branch("F5PPAC2B_Y", &F5PPAC2B_Y, "F5PPAC2B_Y/D");

    raw -> Branch("F7PPAC1A_X", &F7PPAC1A_X, "F7PPAC1A_X/D");raw -> Branch("F7PPAC1A_Y", &F7PPAC1A_Y, "F7PPAC1A_Y/D");
    raw -> Branch("F7PPAC1B_X", &F7PPAC1B_X, "F7PPAC1B_X/D");raw -> Branch("F7PPAC1B_Y", &F7PPAC1B_Y, "F7PPAC1B_Y/D");
    raw -> Branch("F7PPAC2A_X", &F7PPAC2A_X, "F7PPAC2A_X/D");raw -> Branch("F7PPAC2A_Y", &F7PPAC2A_Y, "F7PPAC2A_Y/D");
    raw -> Branch("F7PPAC2B_X", &F7PPAC2B_X, "F7PPAC2B_X/D");raw -> Branch("F7PPAC2B_Y", &F7PPAC2B_Y, "F7PPAC2B_Y/D");

    raw -> Branch("ICQ", &ICQ, "ICQ[8]/I");

    raw -> Branch("F3PLA_TL", &F3PLA_TL, "F3PLA_TL/D");raw -> Branch("F3PLA_TR", &F3PLA_TR, "F3PLA_TR/D");
    raw -> Branch("F3PLA_QL", &F3PLA_QL, "F3PLA_QL/I");raw -> Branch("F3PLA_QR", &F3PLA_QR, "F3PLA_QR/I");
    raw -> Branch("F3PLA_DT", &F3PLA_DT, "F3PLA_DT/D");raw -> Branch("F3PLA_Q_TEST", &F3PLA_Q_TEST, "F3PLA_Q_TEST/D");

    raw -> Branch("F7PLA_TL", &F7PLA_TL, "F7PLA_TL/D");raw -> Branch("F7PLA_TR", &F7PLA_TR, "F7PLA_TR/D");
    raw -> Branch("F7PLA_QL", &F7PLA_QL, "F7PLA_QL/I");raw -> Branch("F7PLA_QR", &F7PLA_QR, "F7PLA_QR/I");
    raw -> Branch("F7PLA_DT", &F7PLA_DT, "F7PLA_DT/D");raw -> Branch("F7PLA_Q_TEST", &F7PLA_Q_TEST, "F7PLA_Q_TEST/D");

    raw -> Branch("F13_1PLA_TL", &F13_1PLA_TL, "F13_1PLA_TL/D");raw -> Branch("F13_1PLA_TR", &F13_1PLA_TR, "F13_1PLA_TR/D");
    raw -> Branch("F13_1PLA_QL", &F13_1PLA_QL, "F13_1PLA_QL/I");raw -> Branch("F13_1PLA_QR", &F13_1PLA_QR, "F13_1PLA_QR/I");
    raw -> Branch("F13_1PLA_DT", &F13_1PLA_DT, "F13_1PLA_DT/D");raw -> Branch("F13_1PLA_Q_TEST", &F13_1PLA_Q_TEST, "F13_1PLA_Q_TEST/D");

    raw -> Branch("F13_2PLA_TL", &F13_2PLA_TL, "F13_2PLA_TL/D");raw -> Branch("F13_2PLA_TR", &F13_2PLA_TR, "F13_2PLA_TR/D");
    raw -> Branch("F13_2PLA_QL", &F13_2PLA_QL, "F13_2PLA_QL/I");raw -> Branch("F13_2PLA_QR", &F13_2PLA_QR, "F13_2PLA_QR/I");
    raw -> Branch("F13_2PLA_DT", &F13_2PLA_DT, "F13_2PLA_DT/D");raw -> Branch("F13_2PLA_Q_TEST", &F13_2PLA_Q_TEST, "F13_2PLA_Q_TEST/D");

  */
  ///multi-hit information from TDCs/////
  int nppachit = 0;
  int nppachit_f7_0 = 0;
  int nppachit_f7_1 = 0;
  int ppactdc_first[128];
  int ppactdc_second[128];
  int f7ppac1_nhit_first = 0;
  int f7ppac1_nhit_second = 0;
  Double_t f7ppac1_firsttdc_ave = 0;
  Double_t f7ppac1_secondtdc_ave = 0;



  //TBeam branches
  Int_t neve=0; TBeam-> Branch("neve", &neve, "neve/I");
  Double_t z=-9999; TBeam -> Branch("z", &z, "z/D");
  Double_t aoq=-9999; TBeam -> Branch("aoq", &aoq, "aoq/D");
  Double_t beta37=-9999; TBeam -> Branch("beta37", &beta37, "beta37/D");
  Double_t brho78=-9999; TBeam -> Branch("brho78", &brho78, "brho78/D");
  Double_t beta=-9999; TBeam -> Branch("beta", &beta, "beta/D");
  Double_t brho=-9999; TBeam -> Branch("brho", &brho, "brho/D");
  Bool_t isGood=false; TBeam -> Branch("isGood", &isGood, "isGood/B");
  Int_t intZ=-9999; TBeam-> Branch("intZ", &intZ, "intZ/I");
  Int_t intA=-9999; TBeam-> Branch("intA", &intA, "intA/I");

  //TBDC branches
  Double_t bdc1x=-9999; TBDC -> Branch("bdc1x",&bdc1x,"bdc1x/D");
  Double_t bdc1y=-9999; TBDC -> Branch("bdc1y",&bdc1y,"bdc1y/D");
  Double_t bdc2x=-9999; TBDC -> Branch("bdc2x",&bdc2x,"bdc2x/D");
  Double_t bdc2y=-9999; TBDC -> Branch("bdc2y",&bdc2y,"bdc2y/D");
  Double_t bdcax=-9999; TBDC -> Branch("bdcax",&bdcax,"bdcax/D");
  Double_t bdcby=-9999; TBDC -> Branch("bdcby",&bdcby,"bdcby/D");
  Double_t ProjX=-9999; TBDC -> Branch("ProjX",&ProjX,"ProjX/D");
  Double_t ProjY=-9999; TBDC -> Branch("ProjY",&ProjY,"ProjY/D");
  Double_t ProjZ=-9999; TBDC -> Branch("ProjZ",&ProjZ,"ProjZ/D");
  Double_t ProjA=-9999; TBDC -> Branch("ProjA",&ProjA,"ProjA/D");
  Double_t ProjB=-9999; TBDC -> Branch("ProjB",&ProjB,"ProjB/D");
  Double_t ProjE=-9999; TBDC -> Branch("ProjE",&ProjE,"ProjE/D");
  Double_t ProjBeta=-9999; TBDC -> Branch("ProjBeta",&ProjBeta,"ProjBeta/D");
  Double_t ProjP=-9999; TBDC -> Branch("ProjP",&ProjP,"ProjP/D");
  Double_t ProjPX=-9999; TBDC -> Branch("ProjPX",&ProjPX,"ProjPX/D");
  Double_t ProjPY=-9999; TBDC -> Branch("ProjPY",&ProjPY,"ProjPY/D");
  Double_t ProjPZ=-9999; TBDC -> Branch("ProjPZ",&ProjPZ,"ProjPZ/D");

  //Focal Plane Branch
  Double_t F7X; TFocalPlane -> Branch("F7X", &F7X, "F7X/D");
  Double_t F7A; TFocalPlane -> Branch("F7A", &F7A, "F7A/D");
  Double_t F7Y; TFocalPlane -> Branch("F7Y", &F7Y, "F7Y/D");
  Double_t F7B; TFocalPlane -> Branch("F7B", &F7B, "F7B/D");

  //////////Cut files//////////////////
  //load plastic cut files


  TFile* cF3cut = new TFile("cut/plastic_cuts/F3cut.root");
  TCutG *F3pl;
  cF3cut->GetObject("CUTG",F3pl);
  cF3cut->Close();

  TFile* cF7cut = new TFile("cut/plastic_cuts/F7cut.root");
  TCutG *F7pl;
  cF7cut->GetObject("CUTG",F7pl);
  cF7cut->Close();

  TFile* cF13cut1 = new TFile("cut/plastic_cuts/F13_1cut.root");
  TCutG *F13pl1;
  cF13cut1->GetObject("CUTG",F13pl1);
  cF13cut1->Close();

  TFile* cF13cut2 = new TFile("cut/plastic_cuts/F13_2cut.root");
  TCutG *F13pl2;
  cF13cut2->GetObject("CUTG",F13pl2);
  cF13cut2->Close();


  //Defining BDCs for SAMURAI
  char myname[128];
  TFile *bdcin = new TFile("./dctpf/dc_tpf.root", "READ");
  if (bdcin->IsOpen()){
    std::cout << "open dc_tpf.root" << std::endl;
    gROOT->cd();
    TH2* hist = NULL;

    for(int i=0;i<8;i++){
      hist = (TH2F*)bdcin->Get(Form("bdc1_tdc_l%02d",i));
      calibbdc1tr->SetTDCDistribution(hist->ProjectionX(),i);
      delete hist; hist = NULL;
    }

    for(int i=0;i<8;i++){
      sprintf(myname,"bdc2_tdc_l%02d",i);
      hist = (TH2F*)bdcin->Get(myname);
      calibbdc2tr->SetTDCDistribution(hist->ProjectionX(),i);
      delete hist; hist = NULL;
    }
  }
  delete bdcin;

  TBDCProjection *bdcProj = new TBDCProjection();
  bdcProj->setBeam(runNo);

  Int_t numEvents = 0;
  while (eventStore -> GetNextEvent() && numEvents < maxEvents) {
    if (numEvents%100 == 0) {
      std::cout << "\r event: " << numEvents <<" / " << maxEvents \
	              << " (" << 100.*numEvents/maxEvents << "%)" << std::flush;
    }

    calibPID -> ClearData();
     recoPID -> ClearData();
    calibPID -> ReconstructData();
     recoPID -> ReconstructData();

    calibbdc1hit->ClearData();
    calibbdc1tr->ClearData();
    calibbdc1hit->ReconstructData();
    calibbdc1tr->ReconstructData();
    //
    ///initialize variables for multi-hit information
    //////following routine provided by Aki///
    nppachit = 0;
    nppachit_f7_0 = 0;
    nppachit_f7_1 = 0;
    for(int k=0;k<128;k++) ppactdc_first[k] = -1;
    for(int k=0;k<128;k++) ppactdc_second[k] = -1;
    // start to scan raw event tree
    for(int i=0;i<rawevent->GetNumSeg();i++){
      TArtRawSegmentObject *seg = rawevent->GetSegment(i);
      //int device = seg->GetDevice();
      //int fp = seg->GetFP();
      int detector = seg->GetDetector();
      //int module = seg->GetModule();

      if(detector == 1) {//PPACT

        for(int j=0;j<seg->GetNumData();j++){
          TArtRawDataObject *d = seg->GetData(j);
          int ch = d->GetCh();
          int val = d->GetVal();
          if(val>95000) continue;
          if(ppactdc_first[ch]<0){
            ppactdc_first[ch]=val;
          }
          else if(ppactdc_second[ch]<0){
            ppactdc_second[ch]=val;
          }

          nppachit ++;
          //if(ch>=48&&ch<64) nppachit_f7 ++;
          if(ch>=48&&ch<56){
            nppachit_f7_0 ++;
            //hf7ppac1tdc->Fill(val);
          }
          if(ch>=56&&ch<64) nppachit_f7_1 ++;
        }
      } // end of PPACT

    }
    f7ppac1_nhit_first = 0;
    f7ppac1_nhit_second = 0;
    f7ppac1_firsttdc_ave = 0;
    f7ppac1_secondtdc_ave = 0;
    for(int c=48;c<56;c++){
      if(ppactdc_first[c]>0){
        f7ppac1_nhit_first ++;
        f7ppac1_firsttdc_ave += ppactdc_first[c];
      }
      if(ppactdc_second[c]>0){
        f7ppac1_nhit_second ++;
        f7ppac1_secondtdc_ave += ppactdc_second[c];
      }
    }
    if(f7ppac1_nhit_first>0)
    f7ppac1_firsttdc_ave /= (double)f7ppac1_nhit_first;
    if(f7ppac1_nhit_second>0)
    f7ppac1_secondtdc_ave /= (double)f7ppac1_nhit_second;
    ///end Aki's section//////////

    //fill some TBeam items
    intZ=-9999;
    intA=-9999;
    z=beam_br37->GetZet(); neve=numEvents; aoq=beam_br37 -> GetAoQ();
    beta37=tof37 -> GetBeta();
    brho78=rips57->GetBrho();//beam_br57->GetBrho();


    Float_t tx1 = -9999;
    Float_t ty1 = -9999;

    auto bdc1trks = (TClonesArray *)storeMan->FindDataContainer("SAMURAIBDC1Track");
    if (bdc1trks) {
      Int_t bdc1ntr = bdc1trks -> GetEntries();
      for (Int_t itr = 0; itr < bdc1ntr; ++itr) {
        auto trk = (TArtDCTrack *) bdc1trks -> At(itr);
        if (trk -> GetPosition(0) > -9999){
          tx1 = trk -> GetPosition(0);
	  bdc1x=tx1-0.72;
        } else if (trk -> GetPosition(1) > -9999){
          ty1 = trk -> GetPosition(1);
	  bdc1y=ty1;
        }
      }
    }

    calibbdc2hit->ClearData();
    calibbdc2tr->ClearData();
    calibbdc2hit->ReconstructData();
    calibbdc2tr->ReconstructData();

    Float_t tx2 = -9999;
    Float_t ty2 = -9999;
    auto bdc2trks = (TClonesArray *)storeMan->FindDataContainer("SAMURAIBDC2Track");
    if (bdc2trks) {
      Int_t bdc2ntr = bdc2trks -> GetEntries();
      for (Int_t itr = 0; itr < bdc2ntr; ++itr) {
        auto trk = (TArtDCTrack *) bdc2trks -> At(itr);
        if (trk -> GetPosition(0) > -9999){
          tx2 = trk -> GetPosition(0);
	  bdc2x=tx2-0.52;//include BDC offset
        }else if (trk -> GetPosition(1) > -9999){
          ty2 = trk -> GetPosition(1);
	  bdc2y=ty2;
        }
      }
    }
    TBeamEnergy *beamE = new TBeamEnergy(z,aoq,beta37);
    beamE->setBeam(runNo);
    beta=beamE->getBeta();
    brho=beamE->getBrho();
    bdcax=std::atan((bdc2x-bdc1x)/1000.)*1000.;
    bdcby=std::atan((bdc2y-bdc1y)/1000.)*1000.;

    ///insert BDC projection///
    ProjX=-9999;ProjY=-9999;ProjZ=-9999;ProjA=-9999;ProjB=-9999;
    ProjPX=-9999;ProjPY=-9999;ProjPZ=-9999;ProjP=-9999;ProjE=-9999;ProjBeta=-9999;
    double E1;
    E1=beamE->getCorrectedEnergy();
    if(z>0 && z<75 && aoq>1. && aoq<3){
    bdcProj->ProjectParticle(bdc2x, bdc2y, -2160., bdcax, bdcby, z, E1, -580.4,beamE->getMass());//-580.4
    ProjX=bdcProj->getX();
    ProjY=bdcProj->getY();
    ProjZ=bdcProj->getZ();
    ProjA=bdcProj->getA();
    ProjB=bdcProj->getB();
    ProjPX=bdcProj->getPX();
    ProjPY=bdcProj->getPY();
    ProjPZ=bdcProj->getPZ();
    ProjP=bdcProj->getP();
    ProjE=bdcProj->getE();
    ProjBeta=bdcProj->getBeta();
    }
    if(z>49.5 && z <50.5 && aoq>goalAoQ-0.01 && aoq < goalAoQ+0.01){
        ;
    }
    F3PPAC1A_X = -9999; F3PPAC1A_Y = -9999; F3PPAC1B_X = -9999; F3PPAC1B_Y = -9999;
    F3PPAC2A_X = -9999; F3PPAC2A_Y = -9999; F3PPAC2B_X = -9999; F3PPAC2B_Y = -9999;
    F5PPAC1A_X = -9999; F5PPAC1A_Y = -9999; F5PPAC1B_X = -9999; F5PPAC1B_Y = -9999;
    F5PPAC2A_X = -9999; F5PPAC2A_Y = -9999; F5PPAC2B_X = -9999; F5PPAC2B_Y = -9999;
    F7PPAC1A_X = -9999; F7PPAC1A_Y = -9999; F7PPAC1B_X = -9999; F7PPAC1B_Y = -9999;
    F7PPAC2A_X = -9999; F7PPAC2A_Y = -9999; F7PPAC2B_X = -9999; F7PPAC2B_Y = -9999;
    //F3PLA_TL = -9999; F3PLA_TR = -9999; F3PLA_QL = -9999; F3PLA_QR = -9999;
    F7PLA_TL = -9999; F7PLA_TR = -9999; F7PLA_QL = -9999; F7PLA_QR = -9999;
    F13_1PLA_TL = -9999; F13_1PLA_TR = -9999; F13_1PLA_QL = -9999; F13_1PLA_QR = -9999;
    F13_2PLA_TL = -9999; F13_2PLA_TR = -9999; F13_2PLA_QL = -9999; F13_2PLA_QR = -9999;
    F3PLA_DT=-9999;F3PLA_Q_TEST=-9999;F7PLA_DT=-9999;F7PLA_Q_TEST=-9999;
    F13_1PLA_DT=-9999;F13_1PLA_Q_TEST=-9999;F13_2PLA_DT=-9999;F13_2PLA_Q_TEST=-9999;
    memset(ICQ, -9999, sizeof(Int_t)*8);
    F7X = -9999; F7A = -9999; F7Y = -9999; F7B = -9999;

    TArtPPAC *ppac;
    ppac = calibPPAC -> FindPPAC((Char_t *) "F3PPAC-1A");
    if (ppac) { F3PPAC1A_X = ppac -> GetX(); F3PPAC1A_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F3PPAC-1B");
    if (ppac) { F3PPAC1B_X = ppac -> GetX(); F3PPAC1B_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F3PPAC-2A");
    if (ppac) { F3PPAC2A_X = ppac -> GetX(); F3PPAC2A_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F3PPAC-2B");
    if (ppac) { F3PPAC2B_X = ppac -> GetX(); F3PPAC2B_Y = ppac -> GetY(); }

    ppac = calibPPAC -> FindPPAC((Char_t *) "F5PPAC-1A");
    if (ppac) { F5PPAC1A_X = ppac -> GetX(); F5PPAC1A_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F5PPAC-1B");
    if (ppac) { F5PPAC1B_X = ppac -> GetX(); F5PPAC1B_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F5PPAC-2A");
    if (ppac) { F5PPAC2A_X = ppac -> GetX(); F5PPAC2A_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F5PPAC-2B");
    if (ppac) { F5PPAC2B_X = ppac -> GetX(); F5PPAC2B_Y = ppac -> GetY(); }

    ppac = calibPPAC -> FindPPAC((Char_t *) "F7PPAC-1A");
    if (ppac) { F7PPAC1A_X = ppac -> GetX(); F7PPAC1A_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F7PPAC-1B");
    if (ppac) { F7PPAC1B_X = ppac -> GetX(); F7PPAC1B_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F7PPAC-2A");
    if (ppac) { F7PPAC2A_X = ppac -> GetX(); F7PPAC2A_Y = ppac -> GetY(); }
    ppac = calibPPAC -> FindPPAC((Char_t *) "F7PPAC-2B");
    if (ppac) { F7PPAC2B_X = ppac -> GetX(); F7PPAC2B_Y = ppac -> GetY(); }

    //Where we fill plastic data
    TArtPlastic *pla;
    pla = calibPLA -> FindPlastic((Char_t *) "F3pl");
    if (pla) {
      F3PLA_TL = pla -> GetTLRaw(); F3PLA_TR = pla -> GetTRRaw();
      F3PLA_QL = pla -> GetQLRaw(); F3PLA_QR = pla -> GetQRRaw();
      F3PLA_DT=F3PLA_TL-F3PLA_TR; F3PLA_Q_TEST=std::log(F3PLA_QR)-std::log(F3PLA_QL);
    }
    pla = calibPLA -> FindPlastic((Char_t *) "F7pl");
    if (pla) {
      F7PLA_TL = pla -> GetTLRaw(); F7PLA_TR = pla -> GetTRRaw();
      F7PLA_QL = pla -> GetQLRaw(); F7PLA_QR = pla -> GetQRRaw();
      F7PLA_DT=F7PLA_TL-F7PLA_TR; F7PLA_Q_TEST=std::log(F7PLA_QR)-std::log(F7PLA_QL);
    }
    pla = calibPLA -> FindPlastic((Char_t *) "F13pl-1");
    if (pla) {
      F13_1PLA_TL = pla -> GetTLRaw(); F13_1PLA_TR = pla -> GetTRRaw();
      F13_1PLA_QL = pla -> GetQLRaw(); F13_1PLA_QR = pla -> GetQRRaw();
      F13_1PLA_DT=F13_1PLA_TL-F13_1PLA_TR; F13_1PLA_Q_TEST=std::log(F13_1PLA_QR)-std::log(F13_1PLA_QL);
    }
    pla = calibPLA -> FindPlastic((Char_t *) "F13pl-2");
    if (pla) {
      F13_2PLA_TL = pla -> GetTLRaw(); F13_2PLA_TR = pla -> GetTRRaw();
      F13_2PLA_QL = pla -> GetQLRaw(); F13_2PLA_QR = pla -> GetQRRaw();
      F13_2PLA_DT=F13_2PLA_TL-F13_2PLA_TR; F13_2PLA_Q_TEST=std::log(F13_2PLA_QR)-std::log(F13_2PLA_QL);
      }

      hF3corr->Fill(F3PLA_DT,F3PLA_Q_TEST);
      hF7corr->Fill(F7PLA_DT,F7PLA_Q_TEST);
      hF131corr->Fill(F13_1PLA_DT,F13_1PLA_Q_TEST);
      hF132corr->Fill(F13_2PLA_DT,F13_2PLA_Q_TEST);


    TArtIC *ic = calibIC -> FindIC((Char_t *) "F7IC");
    if (ic) for (Int_t i = 0; i < 6; i++) ICQ[i] = ic -> GetRawADC(i);

    TArtFocalPlane *fp = nullptr;
    TVectorD *vec;
    fp = calibFP -> FindFocalPlane(3);
    fp = calibFP -> FindFocalPlane(7);
    if (fp) {
      vec = fp -> GetOptVector();
      F7X = (*vec)(0); F7A = (*vec)(1);
      F7Y = (*vec)(2); F7B = (*vec)(3);
    }

    isGood=true;
    //apply cuts
    F3PLA_in=true;F7PLA_in=true;F13_1PLA_in=true;F13_2PLA_in=true;
    totalPPAChits_lte96=true;F7A_PPAChits_lte8=true;F7B_PPAChits_lte8=true;

    //if(not F3pl->IsInside(F3PLA_DT,F3PLA_Q_TEST)) {isGood=false;F3PLA_in=false;}
    //if(not F7pl->IsInside(F7PLA_DT,F7PLA_Q_TEST)) {isGood=false;F7PLA_in=false;}
    //if(not F13pl1->IsInside(F13_1PLA_DT,F13_1PLA_Q_TEST)) {isGood=false;F13_1PLA_in=false;}
    //if(not F13pl2->IsInside(F13_2PLA_DT,F13_2PLA_Q_TEST)) {isGood=false;F13_2PLA_in=false;}

    if(nppachit>96) {isGood=false;totalPPAChits_lte96=false;}
    if(nppachit_f7_0>8) {isGood=false;F7A_PPAChits_lte8=false;}
    if(nppachit_f7_0>8) {isGood=false;F7B_PPAChits_lte8=false;}


    histAllBeamPID->Fill(aoq,z);
    if(isGood){
      intZ=GetIsotope(z,aoq)[0];
      intA=GetIsotope(z,aoq)[1];
      histBeamPID->Fill(aoq,z);
      histProjection->Fill(ProjX,ProjY);

      histBeamEnergy78->Fill(E1);
      histBeamEnergy->Fill(ProjE);
      histBeta->Fill(ProjBeta);
      histBeta37->Fill(beta37);
      histBrho2->Fill(brho78);
      histBrho->Fill(brho);
    }

    if(fill_cuts) cut_tree->Fill();
    //raw -> Fill();
    TBeam -> Fill();
    TBDC -> Fill();
    TFocalPlane -> Fill();

    numEvents++;

  }

  outfile->cd();


  F3pl->Write();
  F7pl->Write();
  F13pl1->Write();
  F13pl2->Write();
  histBeamPID->Write();
  histAllBeamPID->Write();
  histBeamEnergy->Write();
  histBeamEnergy78->Write();
  histProjection->Write();
  histBeta->Write();
  histBeta37->Write();
  histBrho->Write();
  histBrho2->Write();
  if(fill_cuts) {cut_tree->Write();}
  //raw->Write();
  TBeam->Write();
  TBDC->Write();
  TFocalPlane -> Write();
  outfile->Write();
  outfile->Close();

  //gSystem -> Exec(Form("ln -sf ../data/run%d.ridf.root plaData/run%d.ridf.root", runNo, runNo));
}
