//macro to asses BDC information. Starting point.
//Parameters////////////

Double_t BDC1_z=-3160.;//mm, center of BDC1 z in magnet frame
Double_t BDC2_z=-2160.;//mm, center of BDC2 z in magnet frame
Double_t TGT_z=-593.1;//mm, desired projection plane in magnet frame
Double_t B_Field= 0.5; //Tesla
//temporary parameters, eventually to be included from ridf files. For now, assume every particle matches these parameters (not true)
Double_t brho=6.95; //T*m
Double_t beta=0.647;//dimensionless
Double_t Z = 50.;
Double_t AoQ = 2.64;
Double_t mp = 1.66053886e-27; //mass proton in kg
Double_t mass = Z*AoQ*mp;//mass in kg
Double_t c = 2.99e8;//m/s**2
Double_t dt = 1.e-5;//time step in seconds
////////////////////////
double * getBfield(double bx, double by, double bz){//return the magnetic field. Temporarily returns a static field
  static double Bfield[3]={0.,0.,0.};
  if (std::sqrt(bx**2 + bz**2) < 1000.){
    Bfield[0]=0.;
    Bfield[1]=0.5;
    Bfield[2]-0.;
  }
  return Bfield;
}

double * next_step(double mx, double my, double mz, double ma, double mb, double mE, double mmass, double mZ){
  static double pos_ang[6]={mx,my,mz,ma,mb,mE};
  double Efield[3]={0.,0.,0,};
  double mB=getBfield(mx,my,mz);
  double mp = std::sqrt(mE**2-(mmass*c**2)**2)/c;//calculate momentum from energy and mass
  double mpz = mp/std::sqrt(1+std::tan(ma)**2+std::tan(mb)**2);//z component of momentum
  double mom[3]={mpz*std::tan(ma),mpz*std::tan(mb),mpz};//momentum broken into three components
  double mv[3];
  mv[0]=mom[0]/mmass+mZ*(Efield[0]+(mom[1]*mB[2]-mom[2]*mB[1])/mmass)*dt/mmass;
  mv[1]=mom[1]/mmass+mZ*(Efield[1]+(mom[2]*mB[0]-mom[0]*mB[2])/mmass)*dt/mmass;
  mv[2]=mom[2]/mmass+mZ*(Efield[2]+(mom[0*mB[1]-mom[1]*mB[0])/mmass)*dt/mmass;


  return pos_ang;
}

void BDCprojection(Int_t runNo = 3202, Int_t neve_max=30000000)
{

  //Output file and Trees to write out
  TFile *fout = new TFile(Form("./output/BDC/BDCout.%i.root",runNo),"recreate");
  auto TGT_lin = new TTree("TGT_lin","TGT_lin");
  auto TGT_mag = new TTree("TGT_mag","TGT_mag");
  auto BDC = new TTree("BDC","BDC");

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


  TArtStoreManager *sman = TArtStoreManager::Instance();

  auto cvs = new TCanvas("cvs", "", 1200, 500);
  cvs -> Divide(3, 1);

  auto cvs2 = new TCanvas("cvs2", "", 1200, 500);
  cvs2 -> Divide(3, 1);



  int neve = 0;
  BDC -> Branch("neve",&neve);

  Double_t TGT_x_0T=-9999; //mm
  Double_t TGT_y_0T=-9999; //mm
  Double_t TGT_a_0T=-999; //mrad
  Double_t TGT_b_0T=-999; //mrad
  Double_t TGT_px_0=-9999;//MeV/c
  Double_t TGT_py_0T=-9999;//MeV/c
  Double_t TGT_pz_0T=-9999;//MeV/c
  Double_t TGT_x_0_5T=-9999; //mm
  Double_t TGT_y_0_5T=-9999; //mm
  Double_t TGT_a_0_5T=-999; //mrad
  Double_t TGT_b_0_5T=-999; //mrad
  Double_t TGT_px_0_5T=-9999;//MeV/c
  Double_t TGT_py_0_5T=-9999;//MeV/c
  Double_t TGT_pz_0_5T=-9999;//MeV/c

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
  Double_t bdc1trx=-9999.;
  Double_t bdc1try=-9999.;
  Double_t bdc2trax, bdc2tray;
  Double_t bdc2trx=-9999.;
  Double_t bdc2try=-9999.;
  BDC -> Branch("bdc1trx",&bdc1trx,"bdc1trx/D");
  BDC -> Branch("bdc1try",&bdc1trx,"bdc1try/D");
  BDC -> Branch("bdc2trx",&bdc1trx,"bdc2trx/D");
  BDC -> Branch("bdc2try",&bdc1trx,"bdc2try/D");

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

    //----------------------------------------------------------
    // Target
    TGT_x_0T=-9999; //mm
    TGT_y_0T=-9999; //mm
    TGT_a_0T=-999; //mrad
    TGT_b_0T=-999; //mrad
    TGT_px_0=-9999;//MeV/c
    TGT_py_0T=-9999;//MeV/c
    TGT_pz_0T=-9999;//MeV/c
    TGT_x_0_5T=-9999; //mm
    TGT_y_0_5T=-9999; //mm
    TGT_a_0_5T=-999; //mrad
    TGT_b_0_5T=-999; //mrad
    TGT_px_0_5T=-9999;//MeV/c
    TGT_py_0_5T=-9999;//MeV/c
    TGT_pz_0_5T=-9999;//MeV/c



    Double_t dist_BDCs = BDC2_z-BDC1_z; //mm
    Double_t dist_BDC1_TGT = TGT_z-BDC1_z; //mm
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

      }

    }

    TGT_lin -> Fill();
    TGT_mag -> Fill();
    BDC -> Fill();
    //move to next event
    estore->ClearData();
    ++neve;


  }//end of event loop

  cvs -> cd(1);
  htgt2xy0T -> Draw("colz");

  cvs -> cd(2);
  htgt2xa0T -> Draw("colz");

  cvs -> cd(3);
  htgt2yb0T -> Draw("colz");
  fout -> Write();
  
  return;
}