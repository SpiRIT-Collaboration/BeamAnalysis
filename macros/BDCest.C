//macro to asses BDC information. Starting point.

//parameters
Double_t BDC1_z=-3159.28;//mm, center of BDC1 z in magnet frame
Double_t BDC2_z=-2158.73;//mm, center of BDC2 z in magnet frame
Double_t BDC1_x=-0.563;//mm, center of BDC1 x in magnet frame
Double_t BDC2_x=0.436;//mm, center of BDC2 x in magnet frame
Double_t TGT_z=-593.1;//mm, desired projection plane in magnet frame
Double_t AC_z=-820.;//mm, desired projection plane in magnet frame
Double_t AC_left=-3.;//BL side of AC
Double_t AC_right=26.;//BR side of AC
Double_t AC_up=19.;//side of AC
Double_t AC_down=-19.;//side of AC
Double_t TGT_left=-15.;//side of TGT
Double_t TGT_right=15.;//side of TGT
Double_t TGT_up=20.;//side of TGT
Double_t TGT_down=-20.;//side of TGT
Double_t dist_BDCs = BDC2_z-BDC1_z; //mm
Double_t dist_BDC1_TGT = TGT_z-BDC1_z; //mm
Double_t dz=1.;


Double_t *MagStep(Double_t Mdz,Double_t MBrho,Double_t MB,Double_t Ma){
  //only for positive charge in +y magnetic field
  Double_t static Arr[2];//output:dx, da in mm, mrad
  Double_t Mya=Ma;//exit angle, mrad (linear approximation)
  Arr[0]=-Mdz*std::tan(Mya/1000.);//dx, mm - this is a linear approximation
  if(abs(MB)>0.){
    Double_t Mrho=MBrho/MB*1000.;//mm
    Mya=std::asin((Mdz+Mrho*std::sin(Ma/1000.))/Mrho)*1000.;//exit angle, mrad
    Arr[0]=Mrho*(std::cos(Mya/1000.)-std::cos(Ma/1000.));//std::sqrt(Mrho*Mrho-Mdzp*Mdzp)-Mrho*std::cos(Ma/1000.);//dx, mm - this is a detailed iteration
  }
  Arr[1]=Mya;
  return Arr;
}

void BDCest(Int_t runNo = 3202, Int_t neve_max=30000000)
{

  //Output file and Trees to write out
  ifstream Bfield;
  Bfield.open("../ReducedBMap.txt");
  Double_t xx[300],yy[300],zz[300],Bxx[300],Byy[300],Bzz[300];
  int ii=0;
  while(ii<=300){
	  Bfield >>xx[ii]>>yy[ii]>>zz[ii]>>Bxx[ii]>>Byy[ii]>>Bzz[ii];
    ii++;
    if (!Bfield.good()) break;
  }
  Bfield.close();
  TFile *fout = new TFile(Form("./output/BDC/BDCout.%i.root",runNo),"recreate");
  auto TGT_lin = new TTree("TGT_lin","TGT_lin");
  auto TGT_mag = new TTree("TGT_mag","TGT_mag");
  auto bdc_info = new TTree("bdc_info","bdc_info");

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
    hbdc1dtdd[i] = new TH2D(Form("hbdc1dtdd%d",i),Form("BDC1 dt dd%d",i),300,700,1000, 200,0,3);
  }
  TH1* hbdc2dtdd[8];
  for (int i=0;i<8;++i){
    hbdc2dtdd[i] = new TH2D(Form("hbdc2dtdd%d",i),Form("BDC2 dt dd%d",i),300,700,1000, 200,0,3);
  }

  TH1* hbdc1xy = new TH2D("hbdc1xy", "BDC1 XY",200,-100,100, 200,-100,100); // mm
  TH1* hbdc1xa = new TH2D("hbdc1xa", "BDC1 XA",400,-100,100, 400,-100,100); // angle: mrad
  TH1* hbdc1yb = new TH2D("hbdc1yb", "BDC1 YB",400,-100,100, 400,-100,100); // angle: mrad

  TH1* hbdc2xy = new TH2D("hbdc2xy", "BDC2 XY",200,-100,100, 200,-100,100); // mm
  TH1* hbdc2xa = new TH2D("hbdc2xa", "BDC2 XA",400,-100,100, 400,-100,100); // angle: mrad
  TH1* hbdc2yb = new TH2D("hbdc2yb", "BDC2 YB",400,-100,100, 400,-100,100); // angle: mrad

  TH1* htgt2xy0T = new TH2D("htgt2xy0T", "TGT XY; x (mm); y (mm)",12000,-60,60, 12000,-60,60); // mm
  TH1* htgt2xa0T = new TH2D("htgt2xa0T", "TGT XA; x (mm); x' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* htgt2yb0T = new TH2D("htgt2yb0T", "TGT YB; y (mm); y' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* htgt2xy0_5T = new TH2D("htgt2xy0_5T", "TGT XY; x (mm); y (mm)",200,-100,100, 200,-100,100); // mm
  TH1* htgt2xa0_5T = new TH2D("htgt2xa0T_5", "TGT XA; x (mm); x' (mrad)",200,-100,100, 200, 0, 200); // angle: mrad
  TH1* htgt2yb0_5T = new TH2D("htgt2yb0T_5", "TGT YB; y (mm); y' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* hACxy0_5T = new TH2D("hACxy0_5T", "TGT XY; x (mm); y (mm)",200,-100,100, 200,-100,100); // mm
  TH1* hACxa0_5T = new TH2D("hACxa0T_5", "TGT XA; x (mm); x' (mrad)",200,-100,100, 200, 0, 200); // angle: mrad
  TH1* hACyb0_5T = new TH2D("hACyb0T_5", "TGT YB; y (mm); y' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* hXvXb = new TH2D("hXvXB", "BDC X v TGT X; X tgt (mm); X bdc (mm)",200,-100,100, 200, -100, 100);
  TH1* hXBDCvXAC = new TH2D("hXBDCvXAC", "BDC X v AC X; X AC (mm); X bdc (mm)",200,-100,100, 200, -100, 100);

  TArtStoreManager *sman = TArtStoreManager::Instance();

  auto cvs = new TCanvas("cvs", "linear projection", 1200, 500);
  cvs -> Divide(3, 1);

  auto cvs2 = new TCanvas("cvs2", "Magnetic field projection", 1200, 500);
  cvs2 -> Divide(3, 1);

  auto cvs3 = new TCanvas("cvs3", "AC field projection", 1200, 500);
  cvs3 -> Divide(3, 1);

  auto cvs4 = new TCanvas("cvs4", "X correlation", 800, 800);

  auto cvs5 = new TCanvas("cvs5", "X correlation", 800, 800);


  int neve = 0;
  bdc_info -> Branch("neve",&neve);

  Double_t TGT_x_0T; //mm
  Double_t TGT_y_0T; //mm
  Double_t TGT_a_0T; //mrad
  Double_t TGT_b_0T; //mrad
  Double_t TGT_px_0T;//MeV/c
  Double_t TGT_py_0T;//MeV/c
  Double_t TGT_pz_0T;//MeV/c
  Double_t TGT_x_0_5T; //mm
  Double_t TGT_y_0_5T; //mm
  Double_t TGT_a_0_5T; //mrad
  Double_t TGT_b_0_5T; //mrad
  Double_t TGT_px_0_5T;//MeV/c
  Double_t TGT_py_0_5T;//MeV/c
  Double_t TGT_pz_0_5T;//MeV/c

  Double_t AC_x_0_5T; //mm
  Double_t AC_y_0_5T; //mm
  Double_t AC_a_0_5T; //mrad
  Double_t AC_b_0_5T; //mrad

  TGT_lin -> Branch("TGT_x_0T",&TGT_x_0T,"TGT_x_0T/D");
  TGT_lin -> Branch("TGT_y_0T",&TGT_y_0T,"TGT_y_0T/D");
  TGT_lin -> Branch("TGT_a_0T",&TGT_a_0T,"TGT_a_0T/D");
  TGT_lin -> Branch("TGT_b_0T",&TGT_b_0T,"TGT_b_0T/D");
  TGT_lin -> Branch("TGT_px_0T",&TGT_px_0T,"TGT_px_0T/D");
  TGT_lin -> Branch("TGT_py_0T",&TGT_py_0T,"TGT_py_0T/D");
  TGT_lin -> Branch("TGT_pz_0T",&TGT_pz_0T,"TGT_pz_0T/D");

  TGT_mag -> Branch("TGT_x_0_5T",&TGT_x_0_5T,"TGT_x_0_5T/D");
  TGT_mag -> Branch("TGT_y_0_5T",&TGT_y_0_5T,"TGT_y_0_5T/D");
  TGT_mag -> Branch("TGT_a_0_5T",&TGT_a_0_5T,"TGT_a_0_5T/D");
  TGT_mag -> Branch("TGT_b_0_5T",&TGT_b_0_5T,"TGT_b_0_5T/D");
  TGT_mag -> Branch("TGT_px_0_5T",&TGT_px_0_5T,"TGT_px_0_5T/D");
  TGT_mag -> Branch("TGT_py_0_5T",&TGT_py_0_5T,"TGT_py_0_5T/D");
  TGT_mag -> Branch("TGT_pz_0_5T",&TGT_pz_0_5T,"TGT_pz_0_5T/D");

  Double_t bdc1trax, bdc1tray;
  Double_t bdc1trx;
  Double_t bdc1try;
  Double_t bdc2trax, bdc2tray;
  Double_t bdc2trx;
  Double_t bdc2try;
  bdc_info -> Branch("bdc1trx",&bdc1trx,"bdc1trx/D");
  bdc_info -> Branch("bdc1try",&bdc1trx,"bdc1try/D");
  bdc_info -> Branch("bdc2trx",&bdc1trx,"bdc2trx/D");
  bdc_info -> Branch("bdc2try",&bdc1trx,"bdc2try/D");

  while(estore->GetNextEvent() && neve<neve_max){
    if (neve%100==0){
      std::cout<<"\r event: "<<neve<<" / "
	       <<neve_max
	       <<" ("<<100.*neve/neve_max<<"%)"
	       <<std::flush;
    }
    bdc1trax=-999;
    bdc1tray=-999;
    bdc1trx=-9999.;
    bdc1try=-9999.;
    bdc2trax=-999;
    bdc2tray=-999;
    bdc2trx=-9999.;
    bdc2try=-9999.;

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

      bdc1trx = posx+BDC1_x; bdc1try = posy;
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

      bdc2trx = posx+BDC2_x; bdc2try = posy;
      bdc2trax = angx*1000.; bdc2tray = angy*1000.;

      hbdc2xy->Fill(posx,posy);
      hbdc2xa->Fill(posx,bdc2trax);
      hbdc2yb->Fill(posy,bdc2tray);

    }

    //----------------------------------------------------------
    // Target
    TGT_x_0T=-9999; //mm
    TGT_y_0T=-9999; //mm
    TGT_a_0T=-999; //mrad
    TGT_b_0T=-999; //mrad
    TGT_px_0T=-9999;//MeV/c
    TGT_py_0T=-9999;//MeV/c
    TGT_pz_0T=-9999;//MeV/c
    TGT_x_0_5T=-9999; //mm
    TGT_y_0_5T=-9999; //mm
    TGT_a_0_5T=-999; //mrad
    TGT_b_0_5T=-999; //mrad
    TGT_px_0_5T=-9999;//MeV/c
    TGT_py_0_5T=-9999;//MeV/c
    TGT_pz_0_5T=-9999;//MeV/c

    AC_x_0_5T=-9999; //mm
    AC_y_0_5T=-9999; //mm
    AC_a_0_5T=-999; //mrad
    AC_b_0_5T=-999; //mrad

    //produce linear projection
    if(bdc1trks && bdc2trks){

      if( bdc1trx>-1000 && bdc1try>-1000 && bdc2trx>-1000 && bdc2try>-1000){
	TGT_x_0T=( bdc2trx-bdc1trx )/dist_BDCs*dist_BDC1_TGT + bdc1trx; //mm
  	TGT_y_0T=( bdc2try-bdc1try )/dist_BDCs*dist_BDC1_TGT + bdc1try; //mm
  	TGT_a_0T=atan(( bdc2trx-bdc1trx )/dist_BDCs)*1000.; //mrad
  	TGT_b_0T=atan(( bdc2try-bdc1try )/dist_BDCs)*1000.; //mrad
  	htgt2xy0T -> Fill(TGT_x_0T,TGT_y_0T); // mm
  	htgt2xa0T -> Fill(TGT_x_0T,TGT_a_0T); //mrad
  	htgt2yb0T -> Fill(TGT_y_0T,TGT_b_0T); //mrad
  	//magnetic field inclusion

	Double_t x,y,z,a,b;
	Double_t B;
	Double_t Brho=7.;//this is to be determined event by event in coming versions
	x=bdc2trx;
	y=bdc2try+(dist_BDC1_TGT-dist_BDCs)*std::tan(TGT_b_0T/1000.);
	z=BDC2_z;
	a=TGT_a_0T;
	//a=bdc2trax;
	b=TGT_b_0T;
	//b=bdc2tray;
	//TVector3 v1(x/10.,y/10.,z/10.);
	//TVector3 vec=mfield.GetField(v1);
	if(true){
	while(z<AC_z){
	  //v1.SetXYZ(x/10.,y/10.,z/10.);
	  //vec=mfield.GetField(v1);
	  //B=vec(2);
	  B=Byy[(int)(std::sqrt(z*z+x*x)/10.+0.5)];
	  x=x+MagStep(dz,Brho,B,a)[0];
	  a=MagStep(dz,Brho,B,a)[1];
	  z=z+dz;
	}
	AC_x_0_5T=x;
	AC_y_0_5T=y;
	AC_a_0_5T=a;
	AC_b_0_5T=b;
	while(z<TGT_z){
	  //v1.SetXYZ(x/10.,y/10.,z/10.);
	  //vec=mfield.GetField(v1);
    //B=vec(2);
    B=Byy[(int)(std::sqrt(z*z+x*x)/10.+0.5)];
	  x=x+MagStep(dz,Brho,B,a)[0];
	  a=MagStep(dz,Brho,B,a)[1];
	  z=z+dz;
	}

	TGT_x_0_5T=x;
	TGT_y_0_5T=y;
	TGT_a_0_5T=a;
	TGT_b_0_5T=b;

	if( abs(TGT_x_0_5T)>10000) TGT_x_0_5T=-9999;
	if(true){//provide filters for what to add to plot
	htgt2xy0_5T -> Fill(TGT_x_0_5T,TGT_y_0_5T); // mm
	htgt2xa0_5T -> Fill(TGT_x_0_5T,TGT_a_0_5T); //mrad
	htgt2yb0_5T -> Fill(TGT_y_0_5T,TGT_b_0_5T); //mrad
	hACxy0_5T -> Fill(AC_x_0_5T,AC_y_0_5T); // mm
	hACxa0_5T -> Fill(AC_x_0_5T,AC_a_0_5T); //mrad
	hACyb0_5T -> Fill(AC_y_0_5T,AC_b_0_5T); //mrad
	hXvXb->Fill(TGT_x_0_5T,bdc2trx);
	hXBDCvXAC->Fill(AC_x_0_5T,bdc2trx);
	}
	}
      }
    }

    TGT_lin -> Fill();
    TGT_mag -> Fill();
    bdc_info -> Fill();
    //move to next event
    estore->ClearData();
    ++neve;


  }//end of event loop

TLine *AC_up_line=new TLine(AC_left,AC_up,AC_right,AC_up);
TLine *AC_down_line=new TLine(AC_left,AC_down,AC_right,AC_down);
TLine *AC_left_line=new TLine(AC_left,AC_down,AC_left,AC_up);
TLine *AC_right_line=new TLine(AC_right,AC_down,AC_right,AC_up);
AC_up_line->SetLineColor(6);
AC_down_line->SetLineColor(6);
AC_left_line->SetLineColor(6);
AC_right_line->SetLineColor(6);


TLine *TGT_up_line=new TLine(TGT_left,TGT_up,TGT_right,TGT_up);
TLine *TGT_down_line=new TLine(TGT_left,TGT_down,TGT_right,TGT_down);
TLine *TGT_left_line=new TLine(TGT_left,TGT_down,TGT_left,TGT_up);
TLine *TGT_right_line=new TLine(TGT_right,TGT_down,TGT_right,TGT_up);
TGT_up_line->SetLineColor(6);
TGT_down_line->SetLineColor(6);
TGT_left_line->SetLineColor(6);
TGT_right_line->SetLineColor(6);

  cvs -> cd(1);
  htgt2xy0T -> Draw();

  cvs -> cd(2);
  htgt2xa0T -> Draw("colz");

  cvs -> cd(3);
  htgt2yb0T -> Draw("colz");

  cvs2 -> cd(1);
  htgt2xy0_5T -> Draw("colz");
  TGT_up_line->Draw("same");
  TGT_down_line->Draw("same");
  TGT_left_line->Draw("same");
  TGT_right_line->Draw("same");

  cvs2 -> cd(2);
  htgt2xa0_5T -> Draw("colz");

  cvs2 -> cd(3);
  htgt2yb0_5T -> Draw("colz");

  cvs3 -> cd(1);
  hACxy0_5T -> Draw("colz");
  AC_up_line->Draw("same");
  AC_down_line->Draw("same");
  AC_left_line->Draw("same");
  AC_right_line->Draw("same");

  cvs3 -> cd(2);
  hACxa0_5T -> Draw("colz");

  cvs3 -> cd(3);
  hACyb0_5T -> Draw("colz");

  cvs4->cd();
  hXvXb->Draw("colz");

  cvs5->cd();
  hXBDCvXAC->Draw("colz");

  fout->cd();
  TGT_lin->Write();
  TGT_mag->Write();
  bdc_info->Write();
  fout->Write();
  fout->Close();


  return;
}
