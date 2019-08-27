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

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "BeamRaw.h"
#include "BeamBeam.h"

using namespace std;

const char *cutfname = NULL;
bool filecalled = false;
/*
double getTOF(int RunNo){//return custom time of flight from ridf_events.csv
  int r,run;
  double Custom_TOF=280.;
  std::string ModLine;
  ifstream in;
  in.open("./ridf_events.csv");
  in.ignore(10000,'\n');
  cout << "start: "<< RunNo << endl;
  for( std::string line; getline( in, line ); ){
    r = line.find(',');
    run=stoi(line.substr(0,r));
    if(run == RunNo){
      ModLine= line.substr(0,line.rfind(','));
      //cout << ModLine << endl;
      if(ModLine.substr(ModLine.rfind(',')+1).length()>0){
	ModLine= ModLine.substr(0,ModLine.rfind(','));
	ModLine= ModLine.substr(0,ModLine.rfind(','));
	
	Custom_TOF=stod(ModLine.substr(ModLine.rfind(',')+1));
	break;
	
      }

    }
  }
  in.close();
  return Custom_TOF;
}
*/
void aoqByRun(){
  //gROOT->ForceStyle();
    //define the first and last run to evaluate
    int first_run=2819;//2542
    int last_run=3039;//2653

    // where we want the AoQ to be: for 132 Sn, this is 132/50 = 2.64
    double goalAoQ=2.48;//assume 132Sn by default
    if(first_run>=1735 && last_run <2174) goalAoQ=2.64;//Commissioning run
    if(first_run>=2174 && last_run <=2509) goalAoQ=108./50.;
    if(first_run>=2509 && last_run <=2805) goalAoQ=112./50.;
    if(first_run>=2805 && last_run <=3040) goalAoQ=132./50.;
    if(first_run>=3040 && last_run <=3184) goalAoQ=124./50.;


    //define the AoQ range to investigate
    double startAoQ,endAoQ;
    startAoQ=goalAoQ-0.005;
    endAoQ=goalAoQ+0.005;
    cout << startAoQ << "," << endAoQ << endl;
    
    Double_t nentries=0;    
    //the canvas is used to store a graph of the reconstructed AoQ, run by run

    for (int ii=first_run; ii<=last_run; ii++){
      ofstream myfile;
      if(ii==first_run){
	myfile.open (Form("output/AoqByRun.%i.%i.csv",first_run,last_run));
      }
      if(ii>first_run){
	myfile.open (Form("output/AoqByRun.%i.%i.csv",first_run,last_run),std::ios_base::app);
      }  

      ifstream checkfile(Form("./data/run%i.ridf.root",ii));
      if(!checkfile){
      myfile <<ii << "," <<  0.  <<endl;
      continue;
      }
      BeamBeam *beam = new BeamBeam();
      beam->fChainBeam->AddFile(Form("./data/run%i.ridf.root",ii),0,"beam");

      beam->Init();
      nentries=0;
      nentries = beam->fChainBeam->GetEntriesFast();
      //cout << " Entries = " << nentries << endl;

      TCanvas *cvs = new TCanvas("cvs", "", 0, 0, 1200, 600);
      cvs->Divide(2,1);
      auto *histAoQSn = new TH1D("histAoQSn", "", 1000, startAoQ, endAoQ);
      auto *histPID = new TH2D("histPID", "", 1000, startAoQ-0.05, endAoQ+0.05,1000,45,55);
      for(int ientry=0; ientry < nentries; ientry++){
        //if(ientry%100000 == 0) cout << "File " << ii << " read: " << ientry*100./nentries << "%" << endl;
        beam->fChainBeam->GetEvent(ientry);
        histPID->Fill(beam->aoq,beam->z);
        if(beam->z>=49.6 && beam->z<=50.4 && nentries>1) histAoQSn->Fill(beam->aoq);

    }
      
      cvs->cd(1);
      histAoQSn->Draw();
      TF1* fitFuncAoq = new TF1("fitFuncAoq","gaus",startAoQ,endAoQ);
      fitFuncAoq->SetParameter(1,goalAoQ);
      histAoQSn->Fit("fitFuncAoq","Q");
      cvs->cd(2);
      histPID->SetMinimum(1.);
      histPID->Draw("colz");
      cvs->Update();
      cvs->SaveAs(Form("./figures/aoqByRun/AOQ%i.png",ii));


      //myfile<<ii<<","<<nentries<<","<<fitFuncAoq->GetParameter(1)<<","<< (goalAoQ-fitFuncAoq->GetParameter(1))*50. <<","<<myTOF<<endl;

      myfile <<ii << "," <<  (goalAoQ-fitFuncAoq->GetParameter(1)) <<endl;
      cout <<ii << "," <<  (goalAoQ-fitFuncAoq->GetParameter(1))*50 <<endl;
      //if(fitFuncAoq->GetParameter(1)>0.1 && nentries>10000) graph->SetPoint(graph->GetN(),ii,fitFuncAoq->GetParameter(1));

      delete cvs;
      delete beam;
      delete histAoQSn;
      delete histPID;
      myfile.close();
    }
}
