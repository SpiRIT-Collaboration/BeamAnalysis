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
    if(run >= RunNo){
      ModLine= line.substr(0,line.rfind(','));
      if(ModLine.substr(ModLine.rfind(',')+1).length()>0) Custom_TOF=stod(ModLine.substr(ModLine.rfind(',')+1));
      break;
    }
  }
  in.close();
  return Custom_TOF;
}

void aoqByRun(){
    gROOT->ForceStyle();
    double myTOF;
    //define the first and last run to evaluate
    int first_run=2840;
    int last_run=2845;

    // where we want the AoQ to be: for 132 Sn, this is 132/50 = 2.64
    double goalAoQ=2.64;//assume 132Sn by default
    if(first_run>=1735 && last_run <2174) goalAoQ=2.64;//Commissioning run
    if(first_run>=2174 && last_run <2509) goalAoQ=108./50.;
    if(first_run>=2509 && last_run <2805) goalAoQ=112./50.;
    if(first_run>=2805 && last_run <3040) goalAoQ=132./50.;
    if(first_run>=3040 && last_run <3184) goalAoQ=124./50.;


    double guessTOF=0.;

    //define the AoQ range to investigate
    double startAoQ=goalAoQ+.2;
    double endAoQ=goalAoQ-.2;

    ofstream myfile;
    myfile.open (Form("output/AoqByRun.%i.%i.csv",first_run,last_run));
    //the canvas is used to store a graph of the reconstructed AoQ, run by run
    TCanvas *cvs1 = new TCanvas("cvs1","",0,0,800,1200);
    auto graph = new TGraphErrors();
    graph->SetMarkerStyle(5);
    graph->SetMarkerSize(0.8);
    graph->SetLineColor(2);
    graph->SetLineWidth(1);

    BeamBeam *beam = new BeamBeam();

    for (int ii=first_run; ii<=last_run; ii++){
      BeamBeam *beam = new BeamBeam();
      beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",ii),0,"beam");

      beam->Init();

      Long64_t nentries = beam->fChainBeam->GetEntriesFast();
      cout << " Entries = " << nentries << endl;

      TCanvas *cvs = new TCanvas("cvs", "", 0, 0, 1200, 600);
      cvs->Divide(2,1);
      auto *histAoQSn = new TH1D("histAoQSn", "", 1000, startAoQ, endAoQ);
      auto *histPID = new TH2D("histPID", "", 1000, startAoQ-0.05, endAoQ+0.05,1000,30,60);
      myTOF=getTOF(ii);
      for(int ientry=0; ientry < nentries; ientry++){
        if(ientry%100000 == 0) cout << "File " << ii << " read: " << ientry*100./nentries << "%" << endl;
        beam->fChainBeam->GetEvent(ientry);
        histPID->Fill(beam->aoq,beam->z);
        if(beam->z>=49.5 && beam->z<=50.5 && nentries>1) histAoQSn->Fill(beam->aoq);

    }

      cvs->cd(1);
      histAoQSn->Draw();
      TF1* fitFuncAoq = new TF1("fitFuncAoq","gaus",startAoQ,endAoQ);
      fitFuncAoq->SetParameter(1,goalAoQ);
      histAoQSn->Fit("fitFuncAoq");
      cvs->cd(2);
      histPID->Draw("colz");
      cvs->Update();
      cvs->SaveAs(Form("./figures/aoqByRun/AOQ%i.png",ii));


      if(goalAoQ-fitFuncAoq->GetParameter(1)<endAoQ-startAoQ){//using a rudimentary correction to suggest a TOF to bring the reconstructed AoQ closer to the desired AoQ
        guessTOF=(goalAoQ-fitFuncAoq->GetParameter(1))*50.+myTOF;

      }
      myfile<<ii<<","<<nentries<<","<<fitFuncAoq->GetParameter(1)<<","<< guessTOF<<","<<myTOF<<endl;

      if(fitFuncAoq->GetParameter(1)>0.1) graph->SetPoint(graph->GetN(),ii,fitFuncAoq->GetParameter(1));


    }
    cvs1->cd();
    graph->Draw();
    cvs1->SaveAs(Form("./figures/AoqByRun.%i.%i.png",first_run,last_run));
    myfile.close();
}
