//macro to asses BDC information. Starting point.

//parameters
Double_t BDC1_z=-3160.;//mm, center of BDC1 z in magnet frame
Double_t BDC2_z=-2160.;//mm, center of BDC2 z in magnet frame
Double_t BDC1_x=-0.72;//-0.563;//mm, center of BDC1 x in magnet frame
Double_t BDC2_x=-0.52;//0.436;//mm, center of BDC2 x in magnet frame
Double_t TGT_z=-593.1;//mm, desired projection plane in magnet frame
Double_t AC_z=-820.;//mm, desired projection plane in magnet frame
Double_t dist_BDCs = BDC2_z-BDC1_z; //mm
Double_t dist_BDC1_TGT = TGT_z-BDC1_z; //mm
Double_t dz=1.;//step forward in mm
Double_t TPCzero[3]={2,204.3,-580.4};
TH1* hBeta = new TH1D("hBeta", "x(mm); beta",200,.5,1); // angle: mrad


Double_t Initial_momentum(Double_t myQ, Double_t myAoQ, Double_t myBeta){
  Double_t myE=0;
  Double_t myP;
  Double_t myMass;
  //hBeta->Fill(myBeta);
  myMass=myQ*myAoQ*931.494*0.9993774;
  //initial energy
  myE=myMass*(1/std::sqrt(1-myBeta*myBeta));
  //energy loss through F7PPAC1
  myE=myE+0.2228*myQ*myQ/myBeta/myBeta*(-7.67057+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.0045;
  //myE=myE+0.002289*myQ*myQ*(std::log(myBeta*myBeta/(1-myBeta*myBeta))/(myBeta*myBeta)-10.084);
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //Energy loss through Ion Chamber window
  myE=myE+0.223380*myQ*myQ/myBeta/myBeta*(-7.695977714+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.0125;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //Energy loss through Ion Chamber gas
  myE=myE+0.000182*myQ*myQ/myBeta/myBeta*(-7.806010114+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*58.6;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //Energy loss through Ion Chamber electrodes
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.644343291+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.01;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //Energy loss through Ion Chamber window
  myE=myE+0.223380*myQ*myQ/myBeta/myBeta*(-7.695977714+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.0125;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through F7PPAC2
  myE=myE+0.2228*myQ*myQ/myBeta/myBeta*(-7.67057+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.0045;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through F7Plas
  myE=myE+0.171309*myQ*myQ/myBeta/myBeta*(-90.60272126+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.002;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through SBT1mylar
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.300888145+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.0012;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through SBT1
  myE=myE+0.171309*myQ*myQ/myBeta/myBeta*(-15.04321305+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.012;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through SBT1+2mylar
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.300888145+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.0024;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through SBT2
  myE=myE+0.171309*myQ*myQ/myBeta/myBeta*(-15.04321305+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.012;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through SBT2mylar
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.300888145+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.0012;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through BDC1window
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.624580002+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.008;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through BDC1Gas
  myE=myE+0.000017*myQ*myQ/myBeta/myBeta*(-7.136558507+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*9;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through BDC1cathode
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.624580002+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.008;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through BDC2window
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.624580002+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.008;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through BDC2Gas
  myE=myE+0.000017*myQ*myQ/myBeta/myBeta*(-7.136558507+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*9;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  //energy loss through BDC2cathode
  myE=myE+0.222831*myQ*myQ/myBeta/myBeta*(-7.624580002+std::log(myBeta*myBeta/(1-myBeta*myBeta))-1)*0.008;
  myBeta=std::sqrt(1-(myMass/(myE))*(myMass/(myE)));
  myP=myMass*myBeta/std::sqrt(std::abs(1-myBeta*myBeta));
  //End up with momentum after the BDCs
  hBeta->Fill(myBeta);
  return myP;
}


Double_t *MagStep(Double_t Mdz,Double_t MBrho,Double_t MB,Double_t Ma){
  //only for positive charge in +y magnetic field
  Double_t static Arr[2];//output:dx, a in mm, mrad
  Double_t Mya=Ma;
  Arr[0]=-Mdz*std::tan(Mya/1000.);//dx, mm - this is a linear approximation
  if(abs(MB)>0.){
    Double_t Mrho=MBrho/MB*1000.;//mm
    Mya=(Ma-(std::asin(Mdz/Mrho)*1000.));//da, mrad
    Arr[0]=Mrho*(std::cos(Mya/1000.)-std::cos(Ma/1000.));//This is the dx for the step in the given magnetic field
  }
  Arr[1]=Mya;
  return Arr;
}

void BDCprojection(Int_t runNo = 3202, Int_t neve_max=3000000)
{
  //parameters
  Double_t Enc_x_Offset=2.;//offset of enclosure, x axis, in mm
  Double_t FC_x_Offset = 3.5;//offset of field cage, x axis, in mm
  Double_t AC_left=-5.;//BL side of AC
  Double_t AC_right=21;//BR side of AC
  Double_t AC_up=15.;//side of AC
  Double_t AC_down=-17.5;//side of AC
  Double_t TGT_left=-15.;//side of TGT
  Double_t TGT_right=15.;//side of TGT
  Double_t TGT_up=20.;//side of TGT
  Double_t TGT_down=-20.;//side of TGT

  //Apply offsets to parameters
  AC_left=AC_left+Enc_x_Offset;
  AC_right=AC_right+Enc_x_Offset;
  TGT_left=TGT_left+FC_x_Offset;
  TGT_right=TGT_right+FC_x_Offset;

  //create magnetic field map
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

  ///Input ridf.root file to extract beta. If this macro is incorporated into the ridftoroot macro, we should extract beta directly there.
  BeamBeam *beam = new BeamBeam();
  beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",runNo),0,"beam");
  beam->Init();

  //Output file and Trees to write out
  TFile *fout = new TFile(Form("./output/BDC/BDCout.%i.ridf.root",runNo),"recreate");
  auto TPCframe = new TTree("TPCframe","TPCframe");
  auto MAGframe = new TTree("MAGframe","MAGframe");
  auto bdcinfo = new TTree("bdcinfo","bdcinfo");

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

  TH1* htgt2xy0T = new TH2D("htgt2xy0T", "TGT XY; x (mm); y (mm)",200,-100,100, 200,-100,100); // mm
  TH1* htgt2xa0T = new TH2D("htgt2xa0T", "TGT XA; x (mm); x' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* htgt2yb0T = new TH2D("htgt2yb0T", "TGT YB; y (mm); y' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* htgt2xy0_5T = new TH2D("htgt2xy0_5T", "TGT XY; x (mm); y (mm)",200,-100,100, 200,-100,100); // mm
  TH1* htgt2xa0_5T = new TH2D("htgt2xa0T_5", "TGT XA; x (mm); x' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* htgt2yb0_5T = new TH2D("htgt2yb0T_5", "TGT YB; y (mm); y' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad

  TH1* hACxy0_5T = new TH2D("hACxy0_5T", "AC XY; x (mm); y (mm)",200,-100,100, 200,-100,100); // mm
  TH1* hACxa0_5T = new TH2D("hACxa0T_5", "AC XA; x (mm); x' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad
  TH1* hACyb0_5T = new TH2D("hACyb0T_5", "AC YB; y (mm); y' (mrad)",200,-100,100, 200, -100, 100); // angle: mrad

  TArtStoreManager *sman = TArtStoreManager::Instance();


  int inAC=0;
  int inTGT=0;
  //Define variables to write out
  int neve = 0;
  bdcinfo -> Branch("neve",&neve);
  TPCframe -> Branch("neve",&neve);
  MAGframe -> Branch("neve",&neve);

  Double_t AC_x_0_5T,AC_y_0_5T,AC_a_0_5T,AC_b_0_5T;

  Double_t TPCx,TPCy,TPCz,TPCpx,TPCpy,TPCpz,TPCa,TPCb;
  Double_t MAGx,MAGy,MAGz,MAGpx,MAGpy,MAGpz,MAGa,MAGb;
  Double_t beta,beta78;
  TVector3 TPCpos,TPCmom;
  TVector3 MAGpos,MAGmom;

  TPCframe -> Branch("TPCx",&TPCx,"TPCx/D");
  TPCframe -> Branch("TPCy",&TPCy,"TPCy/D");
  TPCframe -> Branch("TPCz",&TPCz,"TPCz/D");
  TPCframe -> Branch("TPCpx",&TPCpx,"TPCpx/D");
  TPCframe -> Branch("TPCpy",&TPCpy,"TPCpy/D");
  TPCframe -> Branch("TPCpz",&TPCpz,"TPCpz/D");
  TPCframe -> Branch("TPCa",&TPCa,"TPCa/D");
  TPCframe -> Branch("TPCb",&TPCb,"TPCb/D");
  TPCframe -> Branch("TPCpos",&TPCpos,"TPCpos/D");
  TPCframe -> Branch("TPCmom",&TPCmom,"TPCmom/D");

  MAGframe -> Branch("MAGx",&MAGx,"MAGx/D");
  MAGframe -> Branch("MAGy",&MAGy,"MAGy/D");
  MAGframe -> Branch("MAGz",&MAGz,"MAGz/D");
  MAGframe -> Branch("MAGpx",&MAGpx,"MAGpx/D");
  MAGframe -> Branch("MAGpy",&MAGpy,"MAGpy/D");
  MAGframe -> Branch("MAGpz",&MAGpz,"MAGpz/D");
  MAGframe -> Branch("MAGa",&MAGa,"MAGa/D");
  MAGframe -> Branch("MAGb",&MAGb,"MAGb/D");
  MAGframe -> Branch("MAGpos",&MAGpos,"MAGpos/D");
  MAGframe -> Branch("MAGmom",&MAGmom,"MAGmom/D");




  TPCframe-> Branch("beta",&beta,"beta/D");
  TPCframe -> Branch("beta",&beta,"beta/D");
  MAGframe -> Branch("beta78",&beta78,"beta78/D");//written in two branches for simplicity. should be same across both branches!
  MAGframe -> Branch("beta78",&beta78,"beta78/D");//written in two branches for simplicity. should be same across both branches!

  Double_t bdc1trax, bdc1tray;
  Double_t bdc1trx,bdc1try;
  Double_t bdc2trax, bdc2tray;
  Double_t bdc2trx,bdc2try;
  bdcinfo -> Branch("bdc1trx",&bdc1trx,"bdc1trx/D");
  bdcinfo -> Branch("bdc1try",&bdc1try,"bdc1try/D");
  bdcinfo -> Branch("bdc2trx",&bdc2trx,"bdc2trx/D");
  bdcinfo -> Branch("bdc2try",&bdc2try,"bdc2try/D");

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

    //initialize variables
    TPCx=-9999; //mm
    TPCy=-9999; //mm
    TPCa=-999; //mrad
    TPCb=-999; //mrad
    TPCpx=-9999;//MeV/c
    TPCpy=-9999;//MeV/c
    TPCpz=-9999;//MeV/c
    MAGx=-9999; //mm
    MAGy=-9999; //mm
    MAGa=-999; //mrad
    MAGb=-999; //mrad
    MAGpx=-9999;//MeV/c
    MAGpy=-9999;//MeV/c
    MAGpz=-9999;//MeV/c

    TPCpos.SetXYZ(-9999,-9999,-9999);
    TPCmom.SetXYZ(-99999,-99999,-99999);

    AC_x_0_5T=-9999; //mm
    AC_y_0_5T=-9999; //mm
    AC_a_0_5T=-999; //mrad
    AC_b_0_5T=-999; //mrad

    Double_t x,y,z,a,b;
    Double_t p=-9999;
    beam->fChainBeam->GetEvent(neve);
    beta78=beam->beta;
    //if(beta78 > 0. && beta78<1.) beta=beta78*0.973;//manually set normalization

    if(beta78 > 0. && beta78<1.)
    {
    p=Initial_momentum(beam->z,beam->aoq, beta78);//in Mev/c
    }
    //produce linear projection
    if(bdc1trks && bdc2trks && beam->z >20 && beam->z < 75 && beam->aoq > 0 && beam->aoq<3){
      if( bdc1trx>-1000 && bdc1try>-1000 && bdc2trx>-1000 && bdc2try>-1000){
        x=( bdc2trx-bdc1trx )/dist_BDCs*dist_BDC1_TGT + bdc1trx; //mm
	      y=( bdc2try-bdc1try )/dist_BDCs*dist_BDC1_TGT + bdc1try; //mm
	      a=atan(( bdc2trx-bdc1trx )/dist_BDCs)*1000.; //mrad
        b=atan(( bdc2try-bdc1try )/dist_BDCs)*1000.; //mrad
      	htgt2xy0T -> Fill(x,y); // mm
      	htgt2xa0T -> Fill(x,a); //mrad
      	htgt2yb0T -> Fill(y,b); //mrad
	//magnetic field inclusion



	Double_t B;
	Double_t Brho;//this is to be determined event by event in coming versions
	Brho=3.3356*p/(std::abs(beam->z))/1000.;//in Tm
	x=bdc2trx;
	y=bdc2try+(AC_z-BDC2_z)*std::tan(b/1000.);
	z=BDC2_z;
	//



	while(z<AC_z){
	  B=Byy[(int)(std::sqrt(z*z+x*x)/10.+0.5)];//pull magnetic field from the previously created map
	  x=x+MagStep(dz,Brho,B,a)[0];//add dx over this step of dz
	  a=MagStep(dz,Brho,B,a)[1];//recalculate angle a after this step of dz
	  z=z+dz;
	}
	AC_x_0_5T=x;
	AC_y_0_5T=y;
	AC_a_0_5T=a;
	AC_b_0_5T=b;

	while(z<TGT_z){
	  B=Byy[(int)(std::sqrt(z*z+x*x)/10.+0.5)];//pull magnetic field from the previously created map
	  x=x+MagStep(dz,Brho,B,a)[0];
	  a=MagStep(dz,Brho,B,a)[1];
	  z=z+dz;
	}
	y=bdc2try+(dist_BDC1_TGT-dist_BDCs)*std::tan(b/1000.);
	TPCpy=p*std::sin(b/1000.);
	TPCpx=std::sqrt(p*p-TPCpy*TPCpy)*std::sin(a/1000.);
	TPCpz=std::sqrt(p*p-TPCpx*TPCpx-TPCpy*TPCpy);

  MAGpos.SetXYZ(x,y,z);
  MAGmom.SetXYZ(TPCpx,TPCpy,TPCpz);

  MAGx=x;
  MAGy=y;
  MAGz=z;

  TPCpos.SetXYZ(x+TPCzero[0],y+TPCzero[1],z+TPCzero[2]);
  TPCmom.SetXYZ(TPCpx,TPCpy,TPCpz);

	TPCx=x+TPCzero[0];
	TPCy=y+TPCzero[1];
  TPCz=z+TPCzero[2];
	TPCa=a;
	TPCb=b;

	if( abs(x)>10000) x=-9999;
	htgt2xy0_5T -> Fill(x,y); // mm
	htgt2xa0_5T -> Fill(x,a); //mrad
	htgt2yb0_5T -> Fill(y,b); //mrad
	hACxy0_5T -> Fill(AC_x_0_5T,AC_y_0_5T); // mm
	hACxa0_5T -> Fill(AC_x_0_5T,AC_a_0_5T); //mrad
	hACyb0_5T -> Fill(AC_y_0_5T,AC_b_0_5T); //mrad

  if(AC_x_0_5T > AC_left && AC_x_0_5T < AC_right && AC_y_0_5T > AC_down && AC_y_0_5T < AC_up ) inAC++;
  if(x > TGT_left && x < TGT_right && y > TGT_down && y < TGT_up ) inTGT++;

      }

    }

    TPCframe -> Fill();
    MAGframe -> Fill();
    bdcinfo -> Fill();
    //move to next event
    estore->ClearData();
    ++neve;


  }//end of event loop
/////////////////create lines to draw target and AC profile/////////////////////////
  auto cvs = new TCanvas("cvs", "linear projection", 1200, 500);
  cvs -> Divide(3, 1);

  auto cvs2 = new TCanvas("cvs2", "Magnetic field projection", 1200, 500);
  cvs2 -> Divide(3, 1);

  auto cvs3 = new TCanvas("cvs3", "AC field projection", 1200, 500);
  cvs3 -> Divide(3, 1);

  auto cvs4 = new TCanvas("cvs4", "Beta distribution after BDCs", 1200, 500);


  cout << Form("Inside AC: %i",inAC)<< endl;
  cout << Form("Inside TGT: %i",inTGT) << endl;


  TLine *AC_up_line=new TLine(AC_left,AC_up,AC_right,AC_up);
  TLine *AC_down_line=new TLine(AC_left,AC_down,AC_right,AC_down);
  TLine *AC_left_line=new TLine(AC_left,AC_down,AC_left,AC_up);
  TLine *AC_right_line=new TLine(AC_right,AC_down,AC_right,AC_up);

  AC_up_line->SetLineColor(6);
  AC_down_line->SetLineColor(6);
  AC_left_line->SetLineColor(6);
  AC_right_line->SetLineColor(6);
  AC_up_line->SetLineWidth(2);
  AC_down_line->SetLineWidth(2);
  AC_left_line->SetLineWidth(2);
  AC_right_line->SetLineWidth(2);

  TLine *TGT_up_line=new TLine(TGT_left,TGT_up,TGT_right,TGT_up);
  TLine *TGT_down_line=new TLine(TGT_left,TGT_down,TGT_right,TGT_down);
  TLine *TGT_left_line=new TLine(TGT_left,TGT_down,TGT_left,TGT_up);
  TLine *TGT_right_line=new TLine(TGT_right,TGT_down,TGT_right,TGT_up);
  TGT_up_line->SetLineColor(6);
  TGT_down_line->SetLineColor(6);
  TGT_left_line->SetLineColor(6);
  TGT_right_line->SetLineColor(6);
  TGT_up_line->SetLineWidth(2);
  TGT_down_line->SetLineWidth(2);
  TGT_left_line->SetLineWidth(2);
  TGT_right_line->SetLineWidth(2);


  cvs -> cd(1);
  htgt2xy0T -> Draw("colz");

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
  hBeta->Draw();

  fout->cd();
  TPCframe->Write();
  MAGframe->Write();
  bdcinfo->Write();
  fout->Write();
  fout->Close();


  return;
}
