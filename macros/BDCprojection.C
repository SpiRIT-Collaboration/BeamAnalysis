//macro to asses BDC information. Starting point.

void BDCprojection(Int_t runNo = 3202, Int_t neve_max=30000000)
{
  //Parameters////////////

  Double_t BDC1_z=-3160.;//mm, center of BDC1 z in magnet frame
  Double_t BDC2_z=-2160.;//mm, center of BDC2 z in magnet frame
  Double_t TGT_z=-593.1;//mm, desired projection plane in magnet frame

  ////////////////////////

  TArtSAMURAIParameters *samurai_prm = new TArtSAMURAIParameters();
  samurai_prm->LoadParameter("db/SAMURAIBDC1.xml");
  samurai_prm->LoadParameter("db/SAMURAIBDC2.xml");

  TArtEventStore *estore = new TArtEventStore();
  //ridf file
  TString ridffilename = Form("ridf/SMDAQ%04d.ridf", runNo);
  estore->Open(ridffilename.Data());// for offline analysis
//  estore->Open(0);// for online analysis
  //drift chamber tpf file (soft linked)
  TFile *bdcin = new TFile("./dctpf/dc_tpf.root","READ");

  TArtCalibBDC1Hit *calibbdc1hit = new TArtCalibBDC1Hit;
  TArtCalibBDC1Track *calibbdc1tr = new TArtCalibBDC1Track;
  TArtCalibBDC2Hit *calibbdc2hit = new TArtCalibBDC2Hit;
  TArtCalibBDC2Track *calibbdc2tr = new TArtCalibBDC2Track;
  calibbdc1tr->SetTDCWindow(0, 2000);
  calibbdc2tr->SetTDCWindow(0, 2000);

  char myname[128];


  if (bdcin->IsOpen()){
    std::cout << "open dc_tpf.root" << std::endl;
    gROOT->cd();
    TH2* hist = NULL;

    for(int i=0;i<8;i++){
//      hist = (TH1D*)bdcin->Get(Form("hbdc1tdc%d",i));
      hist = (TH2F*)bdcin->Get(Form("bdc1_tdc_l%02d",i));
      calibbdc1tr->SetTDCDistribution(hist->ProjectionX(),i);
      delete hist; hist = NULL;
    }

    for(int i=0;i<8;i++){
//      sprintf(myname,"hbdc2tdc%d",i);
      sprintf(myname,"bdc2_tdc_l%02d",i);
      hist = (TH2F*)bdcin->Get(myname);
      calibbdc2tr->SetTDCDistribution(hist->ProjectionX(),i);
      delete hist; hist = NULL;
    }
  }

  gROOT -> cd();

  TH1* hbdc1idtdc = new TH2D("hbdc1idtdc","BDC1 TDC",128,0.5,128.5, 200,0,2000);
  TH1* hbdc2idtdc = new TH2D("hbdc2idtdc","BDC2 TDC",128,0.5,128.5, 200,0,2000);

  TH1* hbdc1dtdd[8];
  for (int i=0;i<8;++i){
    hbdc1dtdd[i] = new TH2D(Form("hbdc1dtdd%d",i),Form("BDC1 dt dd%d",i),
			    300,700,1000, 200,0,3);
  }
  TH1* hbdc2dtdd[8];
  for (int i=0;i<8;++i){
    hbdc2dtdd[i] = new TH2D(Form("hbdc2dtdd%d",i),Form("BDC2 dt dd%d",i),
			    300,700,1000, 200,0,3);
  }

  TH1* hbdc1xy = new TH2D("hbdc1xy", "BDC1 XY",200,-100,100, 200,-100,100); // mm
  TH1* hbdc1xa = new TH2D("hbdc1xa", "BDC1 XA",400,-100,100, 400,-100,100); // angle: mrad
  TH1* hbdc1yb = new TH2D("hbdc1yb", "BDC1 YB",400,-100,100, 400,-100,100); // angle: mrad

  TH1* hbdc2xy = new TH2D("hbdc2xy", "BDC2 XY",200,-100,100, 200,-100,100); // mm
  TH1* hbdc2xa = new TH2D("hbdc2xa", "BDC2 XA",400,-100,100, 400,-100,100); // angle: mrad
  TH1* hbdc2yb = new TH2D("hbdc2yb", "BDC2 YB",400,-100,100, 400,-100,100); // angle: mrad

  TH1* htgt2xy0T = new TH2D("htgt2xy0T", "TGT XY; x (mm); y (mm)",200,-100,100, 200,-100,100); // mm
  TH1* htgt2xa0T = new TH2D("htgt2xa0T", "TGT XA; x (mm); x' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* htgt2yb0T = new TH2D("htgt2yb0T", "TGT YB; y (mm); y' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* hbdc1abdc2a = new TH2D("hbdc1abdc2a", "BDC1A vs BDC2A",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* hbdc1bbdc2b = new TH2D("hbdc1bbdc2b", "BDC1B vs BDC2B",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* hbdc1xbdc2x = new TH2D("hbdc1xbdc2x", "BDC1Xvs BDC2X",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* hbdc1ybdc2y = new TH2D("hbdc1ybdc2y", "BDC1Y vs BDC2Y",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* hbdc1aexta = new TH2D("hbdc1aexta", "BDC1A vs extrapolated angle",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* hbdc1bextb = new TH2D("hbdc1bextb", "BDC1B vs extrapolated angle",200,-100,100, 200, -100, 100); // angle: mrad


  TArtStoreManager *sman = TArtStoreManager::Instance();

  auto cvs = new TCanvas("cvs", "", 1200, 500);
  cvs -> Divide(3, 1);

  auto cvs3 = new TCanvas("cvs3", "", 1200, 500);
  cvs3 -> Divide(3, 1);

  auto cvs4 = new TCanvas("cvs4", "", 1200, 500);
  cvs4 -> Divide(3, 1);

  auto cvs5 = new TCanvas("cvs5", "", 1000, 500);
  cvs5 -> Divide(2, 1);

  auto cvs6 = new TCanvas("cvs6", "", 1000, 500);
  cvs6 -> Divide(2, 1);

  auto cvs7 = new TCanvas("cvs7", "", 1000, 500);
  cvs7 -> Divide(2, 1);

  int neve = 0;
  while(estore->GetNextEvent() && neve<neve_max){
    if (neve%100==0){
      std::cout<<"\r event: "<<neve<<" / "
	       <<neve_max
	       <<" ("<<100.*neve/neve_max<<"%)"
	       <<std::flush;
    }
    //----------------------------------------------------------
    // BDC1
    calibbdc1hit->ClearData();
    calibbdc1tr->ClearData();
    calibbdc1hit->ReconstructData();
    calibbdc1tr->ReconstructData();

    TClonesArray *bdc1hits = (TClonesArray*)calibbdc1hit->GetDCHitArray();
    for (int i=0;i<bdc1hits->GetEntries();++i){
      TArtDCHit* hit = (TArtDCHit*) bdc1hits->At(i);
      hbdc1idtdc->Fill(hit->GetID(), hit->GetTDC());
    }

    Double_t bdc1trx, bdc1try, bdc1trax, bdc1tray;
    Double_t bdc2trx, bdc2try, bdc2trax, bdc2tray;

    TClonesArray *bdc1trks = (TClonesArray *)sman->FindDataContainer("SAMURAIBDC1Track");
    if(bdc1trks){
      //TArtCore::Info(__FILE__,"num bdc1tr %d", bdc1trks->GetEntries());
      Int_t bdc1ntr = bdc1trks->GetEntries();
      Double_t posx=-9999;
      Double_t posy=-9999;
      Double_t angx=-9999;
      Double_t angy=-9999;

      for (int itr=0;itr<bdc1ntr;++itr){
	TArtDCTrack *trk = (TArtDCTrack *)bdc1trks->At(itr);
	if (trk->GetPosition(0)>-9999){
	  posx = trk->GetPosition(0);
	  angx = trk->GetAngle(0);
	}else if (trk->GetPosition(1)>-9999){
	  posy = trk->GetPosition(1);
	  angy = trk->GetAngle(1);
	}
	Int_t nhit = trk->GetNumHitLayer();
	for (int ii=0;ii<nhit;++ii){
	  Int_t hitid = trk->GetHitID(ii);
	  TArtDCHit* hit = (TArtDCHit*)bdc1hits->At(hitid);
	  hbdc1dtdd[hit->GetLayer()]->Fill(hit->GetTDC(),trk->GetDriftLength(ii));
	}
      }

      bdc1trx = posx; bdc1try = posy;
      bdc1trax = angx*1000.; bdc1tray = angy*1000.;

      hbdc1xy->Fill(posx,posy);
      hbdc1xa->Fill(posx,bdc1trax);
      hbdc1yb->Fill(posy,bdc1tray);
    }

    //----------------------------------------------------------
    // BDC2
    calibbdc2hit->ClearData();
    calibbdc2tr->ClearData();
    calibbdc2hit->ReconstructData();
    calibbdc2tr->ReconstructData();
    TClonesArray *bdc2hits = (TClonesArray*)calibbdc2hit->GetDCHitArray();
    for (int i=0;i<bdc2hits->GetEntries();++i){
      TArtDCHit* hit = (TArtDCHit*) bdc2hits->At(i);
      hbdc2idtdc->Fill(hit->GetID(), hit->GetTDC());
    }
    TClonesArray *bdc2trks = (TClonesArray *)sman->FindDataContainer("SAMURAIBDC2Track");
    if(bdc2trks){
      //TArtCore::Info(__FILE__,"num bdc2tr %d", bdc2trks->GetEntries());
      Int_t bdc2ntr = bdc2trks->GetEntries();
      Double_t posx=-9999;
      Double_t posy=-9999;
      Double_t angx=-9999;
      Double_t angy=-9999;
      for (int itr=0;itr<bdc2ntr;++itr){
	TArtDCTrack *trk = (TArtDCTrack *)bdc2trks->At(itr);
	if (trk->GetPosition(0)>-9999){
	  posx = trk->GetPosition(0);
	  angx = trk->GetAngle(0);
	}else if (trk->GetPosition(1)>-9999){
	  posy = trk->GetPosition(1);
	  angy = trk->GetAngle(1);
	}
	Int_t nhit = trk->GetNumHitLayer();
	for (int ii=0;ii<nhit;++ii){
	  Int_t hitid = trk->GetHitID(ii);
	  TArtDCHit* hit = (TArtDCHit*)bdc2hits->At(hitid);
	  hbdc2dtdd[hit->GetLayer()]->Fill(hit->GetTDC(),trk->GetDriftLength(ii));
	}
      }

      bdc2trx = posx; bdc2try = posy;
      bdc2trax = angx*1000.; bdc2tray = angy*1000.;

      hbdc2xy->Fill(posx,posy);
      hbdc2xa->Fill(posx,bdc2trax);
      hbdc2yb->Fill(posy,bdc2tray);

    }
    //---------------------------
    //Investigating correlation between BDC's
    hbdc1abdc2a->Fill(bdc1trax,bdc2trax);
    hbdc1bbdc2b->Fill(bdc1tray,bdc2tray);

    hbdc1xbdc2x->Fill(bdc1trx,bdc2trx);
    hbdc1ybdc2y->Fill(bdc1try,bdc2try);

    hbdc1aexta->Fill(bdc1trax,(bdc2trx-bdc1trx));
    hbdc1bextb->Fill(bdc1tray,(bdc2try-bdc1try));
    //----------------------------------------------------------
    // Target



    Double_t dist_BDCs = BDC2_z-BDC1_z; //mm
    Double_t dist_BDC1_TGT = TGT_z-BDC1_z; //mm
    //produce linear projection
    Double_t TGT_x_0T=( bdc2trx-bdc1trx )/dist_BDCs*dist_BDC1_TGT + bdc1trx; //mm
    Double_t TGT_y_0T=( bdc2try-bdc1try )/dist_BDCs*dist_BDC1_TGT + bdc1try; //mm
    Double_t TGT_a_0T=atan(( bdc2trx-bdc1trx )/dist_BDCs)/1000.; //mrad
    Double_t TGT_b_0T=atan(( bdc2try-bdc1try )/dist_BDCs)/1000.; //mrad

    if( bdc1trx>-1000 && bdc1try>-1000 && bdc2trx>-1000 && bdc2try>-1000){
      htgt2xy0T -> Fill(TGT_x_0T,TGT_y_0T); // mm
      htgt2xa0T -> Fill(TGT_x_0T,TGT_a_0T); //mrad
      htgt2yb0T -> Fill(TGT_y_0T,TGT_b_0T); //mrad
    }

    estore->ClearData();
    ++neve;


  }//end of event loop

  cvs -> cd(1);
  htgt2xy -> Draw("colz");

  cvs -> cd(2);
  htgt2xa -> Draw("colz");

  cvs -> cd(3);
  htgt2yb -> Draw("colz");

  cvs3->cd(1);
  hbdc1xy->Draw("colz");

  cvs3->cd(2);
  hbdc1xa->Draw("colz");

  cvs3->cd(3);
  hbdc1yb->Draw("colz");

  cvs4->cd(1);
  hbdc2xy->Draw("colz");

  cvs4->cd(2);
  hbdc2xa->Draw("colz");

  cvs4->cd(3);
  hbdc2yb->Draw("colz");

  cvs5->cd(1);
  hbdc1abdc2a->Draw("colz");

  cvs5->cd(2);
  hbdc1bbdc2b->Draw("colz");

  cvs6->cd(1);
  hbdc1xbdc2x->Draw("colz");

  cvs6->cd(2);
  hbdc1ybdc2y->Draw("colz");

  cvs7->cd(1);
  hbdc1aexta->Draw("colz");

  cvs7->cd(2);
  hbdc1bextb->Draw("colz");

  TCanvas *cvs2 = new TCanvas("cvs2", "", 800, 400);
  cvs2 -> Divide(2, 1);

  TH1D *projX = ((TH2D *) htgt2xy0T) -> ProjectionX();
  TH1D *projY = ((TH2D *) htgt2xy0T) -> ProjectionY();

  cvs2 -> cd(1);
  projX -> Draw();
  projX -> Fit("gaus");
  Double_t mean = projX -> GetFunction("gaus") -> GetParameter(1);
  projX -> GetXaxis() -> SetRangeUser(mean - 30, mean + 30);
  cvs2 -> cd(2);
  projY -> Draw();
  projY -> Fit("gaus");
  mean = projY -> GetFunction("gaus") -> GetParameter(1);
  projY -> GetXaxis() -> SetRangeUser(mean - 30, mean + 30);

  cvs -> SaveAs(Form("figures/beam_%04d.png", runNo));
  cvs2 -> SaveAs(Form("figures/beam_proj_%04d.png", runNo));


  return;
}
