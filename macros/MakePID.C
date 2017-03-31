//----------------------------------------------------------------------------------//
//-----------------! Macro for PID plot (w/ stats) & fit over A/Q !-----------------//
//---------------!  !!! Change the values depending on beam run !!! !---------------//
//------------------------! Written by C. Santamaria / NSCL !-----------------------//
//----------------------------------------------------------------------------------//


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TChain.h"
#include "TCutG.h"
#include "TEllipse.h"
#include <TFile.h>
#include "TCanvas.h"
#include "TH2D.h"
#include "TH2.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "BeamBeam.h"

using namespace std;

/*
/////////// 132Sn ///////////
double AoQmin = 2.6;
double AoQmax = 2.72;
double dAoQ= 0.005; //range for AoQ fit
double Zetmin = 46;
double Zetmax = 54;
int RunMin = 2840;
int RunMax = 3039;
int nNuclei = 8;
string NameNuclei[] = {"129In","130In","131In","131Sn","132Sn","133Sn","134Sb","135Sb"};
double ANuclei[]    = {129.,130.,131.,131.,132.,133.,134.,135.};
double ZNuclei[]    = {49,  49,  49,  50,  50,  50,  51,  51};
*/

/////////// 124Sn ///////////
double AoQmin = 2.33;
double AoQmax = 2.65;
double dAoQ= 0.005;
double Zetmin = 46;
double Zetmax = 54;
int RunMin = 3058;
int RunMax = 3184;
int nNuclei = 9;
string NameNuclei[] = {"122In","123In","124In","123Sn","124Sn","125Sn","124Sb","125Sb","126Sb"};
double ANuclei[]    = {122.,123.,124.,123.,124.,125.,124.,125.,126.};
double ZNuclei[]    = {49,  49,  49,  50,  50,  50,  51,  51,  51};


void MakePID(){
    gROOT->ForceStyle();
    
    BeamBeam *beam = new BeamBeam();
    for (int ii=RunMin; ii<=RunMax; ii++){
        //int ii = 3058;
        beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",ii),0,"beam");
    }
    beam->Init();
    //cout << "try" << endl;
    Long64_t nentries = beam->fChainBeam->GetEntriesFast();
    cout << " Entries = " << nentries << endl;
    
    
    TCanvas *cvsPID = new TCanvas("cvsPID", "", 0, 0, 800, 580);
    TCanvas *cvsAoQSn = new TCanvas("cvsAoQSn", "", 0, 0, 800, 580);
    TCanvas *cvsAoQIn = new TCanvas("cvsAoQIn", "", 0, 0, 800, 580);
    TCanvas *cvsAoQSb = new TCanvas("cvsAoQSb", "", 0, 0, 800, 580);
    TCanvas *cvsAoQdiff = new TCanvas("cvsAoQdiff", "", 0, 0, 800, 580);
    TCanvas *cvsAoQres = new TCanvas("cvsAoQres", "", 0, 0, 800, 580);
    
    auto *histBeamPID = new TH2D("histBeamPID", "", 1500, AoQmin, AoQmax, 1500, Zetmin, Zetmax);
    auto *histAoQSn = new TH1D("histAoQSn", "histAoQSn", 1500, AoQmin, AoQmax);
    auto *histAoQIn = new TH1D("histAoQIn", "histAoQIn", 1500, AoQmin, AoQmax);
    auto *histAoQSb = new TH1D("histAoQSb", "histAoQSb", 1500, AoQmin, AoQmax);
    
    //-------------------- Ellipse values --------------------//
    const Int_t np = 20; // number of points in the ellipse, this seems enough here, but we can decide when you do it for all nuclei...
    double kPI = TMath::Pi();
    double r1 = 0.003; // A/Q spread, shouldn't have to be changed, but should be checked once we have it for all nuclei
    double r2 = 0.5; // Z spread, same than above
    double theta = 0.; // tilt of the Ellipse in angle, don't need it here
    
    double AoqNuclei[nNuclei];
    int countsNuclei[nNuclei];
    for(int in=0; in<nNuclei; in++){
        AoqNuclei[in] = ANuclei[in]/ZNuclei[in];
        countsNuclei[in] = 0;
    }
    
    vector<TEllipse *> ellipseNuclei;
    vector<TCutG *> cutsNuclei;
    
    //----------------------- Loop over entries -----------------------//
    for(int ientry=0; ientry < nentries; ientry++){
        if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
        beam->fChainBeam->GetEvent(ientry);
        
        //--------------------- Create Histograms ---------------------//
        double AoQ = beam->aoq;
        double Zet = beam->z;
        histBeamPID->Fill(AoQ,Zet);
        if(Zet>=48.5 && Zet<=49.5) histAoQIn->Fill(AoQ);
        if(Zet>=49.5 && Zet<=50.5) histAoQSn->Fill(AoQ);
        if(Zet>=50.5 && Zet<=51.5) histAoQSb->Fill(AoQ);
        
        for(int NN=0;NN<nNuclei;NN++){ //loop over nuclei i.e. cuts
            //-------------------- Ellipse to TCutG --------------------//
            static Double_t x[np+3], y[np+3];
            double x1 = AoqNuclei[NN]; // A/Q mean value
            double y1 = ZNuclei[NN]; // Z value
            
            ellipseNuclei.push_back(new TEllipse(x1, y1, r1, r2));
            ellipseNuclei.back()->SetFillStyle(0);
            ellipseNuclei.back()->SetLineColor(2);
            ellipseNuclei.back()->SetLineWidth(2);
            
            cutsNuclei.push_back(new TCutG());
            cutsNuclei.back()->SetVarX("aoq");
            cutsNuclei.back()->SetVarY("z");
            
            Double_t circ = kPI*(r1+r2);
            Double_t angle,dx,dy;
            Double_t dphi = 2*kPI/np;
            Double_t ct   = TMath::Cos(kPI*theta/180);
            Double_t st   = TMath::Sin(kPI*theta/180);
            for (Int_t i=0;i<=np;i++) {
                angle = Double_t(i)*dphi;
                dx    = r1*TMath::Cos(angle);
                dy    = r2*TMath::Sin(angle);
                x[i]  = x1 + dx*ct - dy*st;
                y[i]  = y1 + dx*st + dy*ct;
                cutsNuclei.back()->SetPoint(i,x[i],y[i]);
            }
        
            if(cutsNuclei.back()->IsInside(beam->aoq,beam->z)) countsNuclei[NN]++;
        }
        
    } // end of loop over entries
    
    //----------------------- Plot Histograms -----------------------//
    cvsPID->cd();
    histBeamPID->Draw("colz");
    histBeamPID->GetXaxis()->SetTitle("A/Q");
    histBeamPID->GetXaxis()->SetTitleSize(0.05);
    histBeamPID->GetXaxis()->SetTitleOffset(0.75);
    histBeamPID->GetXaxis()->SetLabelSize(0.04);
    histBeamPID->GetXaxis()->SetTitleFont(22);
    histBeamPID->GetXaxis()->SetLabelFont(132);
    histBeamPID->GetXaxis()->CenterTitle(true);
    histBeamPID->GetYaxis()->SetTitle("Z");
    histBeamPID->GetYaxis()->SetTitleSize(0.05);
    histBeamPID->GetYaxis()->SetTitleOffset(0.7);
    histBeamPID->GetYaxis()->SetLabelSize(0.04);
    histBeamPID->GetYaxis()->SetTitleFont(22);
    histBeamPID->GetYaxis()->SetLabelFont(132);
    histBeamPID->GetYaxis()->CenterTitle(true);
    gStyle->SetPalette(55);
    cvsPID->Update();

    cvsAoQIn->cd();
    histAoQIn->Draw();
    histAoQIn->GetXaxis()->SetTitle("A/Q");
    histAoQIn->GetXaxis()->SetTitleSize(0.05);
    histAoQIn->GetXaxis()->SetTitleOffset(0.7);
    histAoQIn->GetXaxis()->SetLabelSize(0.04);
    histAoQIn->GetXaxis()->SetTitleFont(22);
    histAoQIn->GetXaxis()->SetLabelFont(132);
    histAoQIn->GetXaxis()->CenterTitle(true);
    histAoQIn->GetYaxis()->SetTitle("Number of counts");
    histAoQIn->GetYaxis()->SetTitleSize(0.05);
    histAoQIn->GetYaxis()->SetTitleOffset(0.7);
    histAoQIn->GetYaxis()->SetLabelSize(0.04);
    histAoQIn->GetYaxis()->SetTitleFont(22);
    histAoQIn->GetYaxis()->SetLabelFont(132);
    histAoQIn->GetYaxis()->CenterTitle(true);
    cvsAoQIn->Update();
    
    cvsAoQSn->cd();
    histAoQSn->Draw();
    histAoQSn->GetXaxis()->SetTitle("A/Q");
    histAoQSn->GetXaxis()->SetTitleSize(0.05);
    histAoQSn->GetXaxis()->SetTitleOffset(0.7);
    histAoQSn->GetXaxis()->SetLabelSize(0.04);
    histAoQSn->GetXaxis()->SetTitleFont(22);
    histAoQSn->GetXaxis()->SetLabelFont(132);
    histAoQSn->GetXaxis()->CenterTitle(true);
    histAoQSn->GetYaxis()->SetTitle("Number of counts");
    histAoQSn->GetYaxis()->SetTitleSize(0.05);
    histAoQSn->GetYaxis()->SetTitleOffset(0.7);
    histAoQSn->GetYaxis()->SetLabelSize(0.04);
    histAoQSn->GetYaxis()->SetTitleFont(22);
    histAoQSn->GetYaxis()->SetLabelFont(132);
    histAoQSn->GetYaxis()->CenterTitle(true);
    cvsAoQSn->Update();
    
    cvsAoQSb->cd();
    histAoQSb->Draw();
    histAoQSb->GetXaxis()->SetTitle("A/Q");
    histAoQSb->GetXaxis()->SetTitleSize(0.05);
    histAoQSb->GetXaxis()->SetTitleOffset(0.7);
    histAoQSb->GetXaxis()->SetLabelSize(0.04);
    histAoQSb->GetXaxis()->SetTitleFont(22);
    histAoQSb->GetXaxis()->SetLabelFont(132);
    histAoQSb->GetXaxis()->CenterTitle(true);
    histAoQSb->GetYaxis()->SetTitle("Number of counts");
    histAoQSb->GetYaxis()->SetTitleSize(0.05);
    histAoQSb->GetYaxis()->SetTitleOffset(0.7);
    histAoQSb->GetYaxis()->SetLabelSize(0.04);
    histAoQSb->GetYaxis()->SetTitleFont(22);
    histAoQSb->GetYaxis()->SetLabelFont(132);
    histAoQSb->GetYaxis()->CenterTitle(true);
    cvsAoQSb->Update();
    
    //--------------------- Show Nuclei Stats & A/Q mean vs Resolution ---------------------//
    double AoQdiff[nNuclei];
    double AoQres[nNuclei];
    double xx[nNuclei];
    for(int ii=0;ii<nNuclei;ii++) xx[ii] = ii;
    for (int in=0; in<nNuclei; in++) {
        cout << " ---------- " << NameNuclei[in] << " (A=" << ANuclei[in] << " & Z=" << ZNuclei[in] << ") ----------" << endl;
        cout << "       Counts = " << countsNuclei[in] <<  ";  " << countsNuclei[in]*100./nentries << " %" << endl;
        cvsPID->cd();
        ellipseNuclei[in]->Draw("same");
        cvsPID->Update();

        cout << "       Fit A/Q = " << endl;
        TF1 *GaussianFit = new TF1("GaussianFit","gaus", AoqNuclei[in]-dAoQ, AoqNuclei[in]+dAoQ);
        GaussianFit->SetParameter(1, AoqNuclei[in]);
        if(ZNuclei[in]==49){
            cvsAoQIn->cd();
            histAoQIn->Fit("GaussianFit","","",AoqNuclei[in]-dAoQ, AoqNuclei[in]+dAoQ);
        }
        else if(ZNuclei[in]==50){
            cvsAoQSn->cd();
            histAoQSn->Fit("GaussianFit","","",AoqNuclei[in]-dAoQ, AoqNuclei[in]+dAoQ);
        }
        else if(ZNuclei[in]==51){
            cvsAoQSb->cd();
            histAoQSb->Fit("GaussianFit","","",AoqNuclei[in]-dAoQ, AoqNuclei[in]+dAoQ);
        }
        //GaussianFit->Draw("same");
        AoQdiff[in] = abs(AoqNuclei[in] - GaussianFit->GetParameter(1))/AoqNuclei[in]*100.;
        AoQres[in] = GaussianFit->GetParameter(2)/AoqNuclei[in]*100.;
        cout << endl;
    }
    

    cvsAoQdiff->cd();
    TGraph *grAoQdiff = new TGraph(nNuclei,xx,AoQdiff);
    grAoQdiff->Draw("AP*");
    grAoQdiff->SetMarkerStyle(20);
    grAoQdiff->SetMarkerSize(2);
    grAoQdiff->GetYaxis()->SetTitle("A/Q difference (%)");
    cvsAoQdiff->Update();
    
    cvsAoQres->cd();
    TGraph *grAoQres = new TGraph(nNuclei,xx,AoQres);
    grAoQres->Draw("AP*");
    grAoQres->SetMarkerStyle(20);
    grAoQres->SetMarkerSize(2);
    grAoQres->GetYaxis()->SetTitle("A/Q resolution (%)");
    cvsAoQres->Update();

}