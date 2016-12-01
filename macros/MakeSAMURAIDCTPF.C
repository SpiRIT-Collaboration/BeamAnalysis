void MakeSAMURAIDCTPF(Int_t runNo, Int_t nanaeve=500000){
  TString ridffile = Form("ridf/sdaq02/SMDAQ%04d.ridf", runNo);

  TArtEventStore *estore = new TArtEventStore;
  estore->Open(ridffile.Data());
//  estore->Open(0);

  TArtSAMURAIParameters *samuraiparameters = TArtSAMURAIParameters::Instance();
  samuraiparameters->LoadParameter("db/SAMURAIBDC1.xml");
  samuraiparameters->LoadParameter("db/SAMURAIBDC2.xml");

  TArtCalibBDC1Hit *fCalibBDC1Hit = new TArtCalibBDC1Hit;
  TArtCalibBDC2Hit *fCalibBDC2Hit = new TArtCalibBDC2Hit;

  TArtStoreManager * sman = TArtStoreManager::Instance();
  TClonesArray *bdc1hits = (TClonesArray *)sman->FindDataContainer("SAMURAIBDC1Hit");
  TClonesArray *bdc2hits = (TClonesArray *)sman->FindDataContainer("SAMURAIBDC2Hit");

  TFile *fout = new TFile(Form("./dctpf/dc_tpf_%04d.root", runNo),"recreate");

  char dummy[32];
  char title[64];

  // define bpc histograms
  int numlayer = 2;
  int numwire = 64;

  // define bdc1 histograms
  numlayer = 8;
  numwire = 16;

  TH2 *hbdc1_ch; // channel distribution
  TH2 *hbdc1_nhit; // number of hits for 8 layers
  TH2 *hbdc1_nwire; // number of wire which hits for 8 layers
  TH2 *hbdc1_nclus; // number of clusters for 8 layers
  TH2 *hbdc1_fch_corr[4]; // corr. of first hit ch between neighboring layer
  TH1 *hbdc1_fch_diff[4]; // diff. of first hit ch between neighboring layer
  TH2 *hbdc1_ftdc_corr[4]; // corr. of first tdc between neighboring layer
  TH1 *hbdc1_ftdc_sum[4]; // sum. of first tdc between neighboring layer
  TH2 *hbdc1_tdc[8]; // tdc distribution for 8 layers

  hbdc1_ch  = new TH2F("bdc1_ch","bdc1 channel distribution",numwire,-0.5,numwire-0.5,numlayer,-0.5,numlayer-0.5);
  hbdc1_nhit  = new TH2F("bdc1_nhit","bdc1 number of hits for 14 layers",10,-0.5,9.5,numlayer,-0.5,numlayer-0.5);
  hbdc1_nwire = new TH2F("bdc1_nwire","bdc1 number of wire which hits for 8 layers",10,-0.5,9.5,numlayer,-0.5,numlayer-0.5);
  hbdc1_nclus = new TH2F("bdc1_nclus","bdc1 number of clusters for 8 layers",10,-0.5,9.5,numlayer,-0.5,numlayer-0.5);

  for(int i=0;i<numlayer/2;i++){
    sprintf(dummy,"bdc1_fch_corr_%d",i);
    sprintf(title,"tdc distribution v.s. channel for layer-%d",i);
    hbdc1_fch_corr[i]  = new TH2F(dummy,title,numwire,-0.5,numwire-0.5,numwire,-0.5,numwire-0.5);
    sprintf(dummy,"bdc1_fch_diff_%d",i);
    sprintf(title,"difference of first hit channel for layer-%d",i);
    hbdc1_fch_diff[i]  = new TH1F(dummy,title,21,-10.5,10.5);
    sprintf(dummy,"bdc1_ftdc_corr_%d",i);
    sprintf(title,"tdc correlation of neighboring channel for layer-%d-%d",i*2,i*2+1);
    hbdc1_ftdc_corr[i]  = new TH2F(dummy,title,2000,0,2000,2000,0,2000);
    sprintf(dummy,"bdc1_ftdc_sum_%d",i);
    sprintf(title,"corr. tdc-sum for layer-%d",i);
    hbdc1_ftdc_sum[i]  = new TH1F(dummy,title,2000,0,4000);
  }
  for(int i=0;i<numlayer;i++){
    sprintf(dummy,"bdc1_tdc_l%02d",i);
    sprintf(title,"bdc1 tdc distribution for layer-%02d",i);
    hbdc1_tdc[i]  = new TH2F(dummy,title,2000,0,2000,numwire,-0.5,numwire-0.5);
  }

  // making bdc2 histograms

  TH2 *hbdc2_ch; // channel distribution
  TH2 *hbdc2_nhit; // number of hits for 8 layers
  TH2 *hbdc2_nwire; // number of wire which hits for 8 layers
  TH2 *hbdc2_nclus; // number of clusters for 8 layers
  TH2 *hbdc2_fch_corr[4]; // corr. of first hit ch between neighboring layer
  TH1 *hbdc2_fch_diff[4]; // diff. of first hit ch between neighboring layer
  TH2 *hbdc2_ftdc_corr[4]; // corr. of first tdc between neighboring layer
  TH1 *hbdc2_ftdc_sum[4]; // sum. of first tdc between neighboring layer
  TH2 *hbdc2_tdc[8]; // tdc distribution for 8 layers

  hbdc2_ch  = new TH2F("bdc2_ch","bdc2 channel distribution",numwire,-0.5,numwire-0.5,numlayer,-0.5,numlayer-0.5);
  hbdc2_nhit  = new TH2F("bdc2_nhit","bdc2 number of hits for 14 layers",10,-0.5,9.5,numlayer,-0.5,numlayer-0.5);
  hbdc2_nwire = new TH2F("bdc2_nwire","bdc2 number of wire which hits for 8 layers",10,-0.5,9.5,numlayer,-0.5,numlayer-0.5);
  hbdc2_nclus = new TH2F("bdc2_nclus","bdc2 number of clusters for 8 layers",10,-0.5,9.5,numlayer,-0.5,numlayer-0.5);

  for(int i=0;i<numlayer/2;i++){
    sprintf(dummy,"bdc2_fch_corr_%d",i);
    sprintf(title,"tdc distribution v.s. channel for layer-%d",i);
    hbdc2_fch_corr[i]  = new TH2F(dummy,title,numwire,-0.5,numwire-0.5,numwire,-0.5,numwire-0.5);
    sprintf(dummy,"bdc2_fch_diff_%d",i);
    sprintf(title,"difference of first hit channel for layer-%d",i);
    hbdc2_fch_diff[i]  = new TH1F(dummy,title,21,-10.5,10.5);
    sprintf(dummy,"bdc2_ftdc_corr_%d",i);
    sprintf(title,"tdc distribution v.s. channel for layer-%d",i);
    hbdc2_ftdc_corr[i]  = new TH2F(dummy,title,2000,0,2000,2000,0,2000);
    sprintf(dummy,"bdc2_ftdc_sum_%d",i);
    sprintf(title,"corr. tdc-sum for layer-%d",i);
    hbdc2_ftdc_sum[i]  = new TH1F(dummy,title,2000,0,4000);
  }
  for(int i=0;i<numlayer;i++){
    sprintf(dummy,"bdc2_tdc_l%02d",i);
    sprintf(title,"bdc2 tdc distribution for layer-%02d",i);
    hbdc2_tdc[i]  = new TH2F(dummy,title,2000,0,2000,numwire,-0.5,numwire-0.5);
  }

  // start to analyze

  Int_t neve = 0;

  while(estore->GetNextEvent() && neve<nanaeve){
    if (neve%100==0){
      std::cout<<"\r event: "<<neve<<" / "
	       <<nanaeve
	       <<" ("<<100.*neve/nanaeve<<"%)"
	       <<std::flush;
    }

    fCalibBDC1Hit->ClearData();
    fCalibBDC2Hit->ClearData();

    fCalibBDC1Hit->ReconstructData();
    fCalibBDC2Hit->ReconstructData();

    //
    // bdc1
    //

    if(bdc1hits){
      int bdc1_nhit[8];
      int bdc1_nwire[8];
      int bdc1_nclus[8];
      bool isbdc1Hit[8][16+2];

      int bdc1_fch[8];
      int bdc1_ftdc[8];
      for(int i=0;i<8;i++){
	bdc1_nhit[i] = 0; bdc1_nwire[i] = 0; bdc1_nclus[i] = 0;
	bdc1_fch[i] = -1; bdc1_ftdc[i] = 99999;
	for(int j=0;j<16+2;j++) isbdc1Hit[i][j] = false;
      }

      for(int i=0;i<bdc1hits->GetEntries();i++){
	TArtDCHit *hit = (TArtDCHit *)bdc1hits->At(i);
	int layer = hit->GetLayer();
	int wireid = hit->GetWireID();
	int val = hit->GetTDC();
	Double_t pos = hit->GetPosition();
	if(val<bdc1_ftdc[layer]){
	  bdc1_ftdc[layer] = val;
	  bdc1_fch[layer] = wireid;
	}

	hbdc1_tdc[layer]->Fill(val,wireid);
	hbdc1_ch->Fill(wireid,layer);
	bdc1_nhit[layer] ++;
	isbdc1Hit[layer][wireid] = true;

      }
      for(int i=0;i<4;i++){
	hbdc1_fch_corr[i]->Fill(bdc1_fch[i*2],bdc1_fch[i*2+1]);
	if(TMath::Abs(bdc1_fch[i*2] - 4) <= 4)
	  hbdc1_fch_diff[i]->Fill(bdc1_fch[i*2]-bdc1_fch[i*2+1]);
	if(bdc1_fch[i*2] == bdc1_fch[i*2+1] && TMath::Abs(bdc1_fch[i*2] - 4) <= 4){
	  hbdc1_ftdc_corr[i]->Fill(bdc1_ftdc[i*2],bdc1_ftdc[i*2+1]);
	  hbdc1_ftdc_sum[i]->Fill(bdc1_ftdc[i*2]+bdc1_ftdc[i*2+1]);
	}
      }

      for(int i=0;i<8;i++)
	for(int j=0;j<16;j++)
	  if(isbdc1Hit[i][j])
	    bdc1_nwire[i] ++;

      for(int i=0;i<8;i++)
	for(int j=0;j<16;j++)
	  if(isbdc1Hit[i][j]){
	    bdc1_nclus[i] ++;
	    if(isbdc1Hit[i][j+1]){
	      if(isbdc1Hit[i][j+2]){
		j += 2;
	      }
	      else{
		j += 1;
	      }
	    }
	  }

      for(int i=0;i<8;i++){
	hbdc1_nhit->Fill(bdc1_nhit[i],i);
	hbdc1_nwire->Fill(bdc1_nwire[i],i);
	hbdc1_nclus->Fill(bdc1_nclus[i],i);
      }
    } // end of if(bdc1hits){

    //
    // bdc2
    //

    if(bdc2hits){
      int bdc2_nhit[8];
      int bdc2_nwire[8];
      int bdc2_nclus[8];
      bool isbdc2Hit[8][16+2];

      int bdc2_fch[8];
      int bdc2_ftdc[8];
      for(int i=0;i<8;i++){
	bdc2_nhit[i] = 0; bdc2_nwire[i] = 0; bdc2_nclus[i] = 0;
	bdc2_fch[i] = -1; bdc2_ftdc[i] = 99999;
	for(int j=0;j<16+2;j++) isbdc2Hit[i][j] = false;
      }

      for(int i=0;i<bdc2hits->GetEntries();i++){
	TArtDCHit *hit = (TArtDCHit *)bdc2hits->At(i);
	int layer = hit->GetLayer();
	int wireid = hit->GetWireID();
	int val = hit->GetTDC();
	Double_t pos = hit->GetPosition();
	if(val<bdc2_ftdc[layer]){
	  bdc2_ftdc[layer] = val;
	  bdc2_fch[layer] = wireid;
	}

	hbdc2_tdc[layer]->Fill(val,wireid);
	hbdc2_ch->Fill(wireid,layer);
	bdc2_nhit[layer] ++;
	isbdc2Hit[layer][wireid] = true;

      }
      for(int i=0;i<4;i++){
	hbdc2_fch_corr[i]->Fill(bdc2_fch[i*2],bdc2_fch[i*2+1]);
	if(TMath::Abs(bdc2_fch[i*2] - 4) <= 4)
	  hbdc2_fch_diff[i]->Fill(bdc2_fch[i*2]-bdc2_fch[i*2+1]);
	if(bdc2_fch[i*2] == bdc2_fch[i*2+1] && TMath::Abs(bdc2_fch[i*2] - 4) <= 4){
	  hbdc2_ftdc_corr[i]->Fill(bdc2_ftdc[i*2],bdc2_ftdc[i*2+1]);
	  hbdc2_ftdc_sum[i]->Fill(bdc2_ftdc[i*2]+bdc2_ftdc[i*2+1]);
	}
      }

      for(int i=0;i<8;i++)
	for(int j=0;j<16;j++)
	  if(isbdc2Hit[i][j])
	    bdc2_nwire[i] ++;

      for(int i=0;i<8;i++)
	for(int j=0;j<16;j++)
	  if(isbdc2Hit[i][j]){
	    bdc2_nclus[i] ++;
	    if(isbdc2Hit[i][j+1]){
	      if(isbdc2Hit[i][j+2]){
		j += 2;
	      }
	      else{
		j += 1;
	      }
	    }
	  }

      for(int i=0;i<8;i++){
	hbdc2_nhit->Fill(bdc2_nhit[i],i);
	hbdc2_nwire->Fill(bdc2_nwire[i],i);
	hbdc2_nclus->Fill(bdc2_nclus[i],i);
      }
    } // end of if(bdc2hits){


  neve++;
  } // end of while(estore->GetNextEvent() && neve<nanaeve){

  fout->Write();
  fout->Close();

}
