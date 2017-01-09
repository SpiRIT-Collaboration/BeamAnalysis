// PPAC: detector-id = 1(PPACT), module-id = 24(V1190)
// BigRIPS Plastic: fpl-id = 13, detector-id = 3(PLAT), module-id = 34(MTDC32)
// BigRIPS IC: detector-id = 4(ICE), module-id = 32(MADC32)

void PPAChits(Int_t runNo){

TString infile =Form("ridf/SMDAQ%i.ridf",runNo);

  gSystem->Load("libanaroot.so");

  TArtEventStore *estore = new TArtEventStore();
  estore->Open(infile);
  TArtRawEventObject *rawevent = estore->GetRawEventObject();

  TH1F *hrpv130bit = new TH1F("rpv130","rpv130",24,-0.5,23.5);
  TH1F *hf7ppac1tdc = new TH1F("f7ppac1tdc","f7ppac1tdc",100000,0,100000);

  int neve = 0;
  int nppachit = 0;
  int nppachit_f7_0 = 0;
  int nppachit_f7_1 = 0;
  int ppactdc_first[128];
  int ppactdc_second[128];
  int f7ppac1_nhit_first = 0;
  int f7ppac1_nhit_second = 0;
  Double_t f7ppac1_firsttdc_ave = 0;
  Double_t f7ppac1_secondtdc_ave = 0;
  Bool_t isGGEarlyClose;
  Bool_t rpv130_bit[24]; // rpv130 bit information

  auto tr = new TTree("tr","tr");
  tr -> Branch("neve", &neve);
  tr -> Branch("nppachit", &nppachit);
  tr -> Branch("nf7ppac1hit", &nppachit_f7_0);
  tr -> Branch("nf7ppac2hit", &nppachit_f7_1);
  //  tr -> Branch("ppactdc_first", ppactdc_first, "ppactdc_first[128]/I");
  //  tr -> Branch("ppactdc_second", ppactdc_second, "ppactdc_second[128]/I");
  tr -> Branch("f7ppac1_nhit_first", &f7ppac1_nhit_first);
  tr -> Branch("f7ppac1_nhit_second", &f7ppac1_nhit_second);
  tr -> Branch("f7ppac1_firsttdc_ave", &f7ppac1_firsttdc_ave);
  tr -> Branch("f7ppac1_secondtdc_ave", &f7ppac1_secondtdc_ave);
  tr -> Branch("ggec", &isGGEarlyClose);

  while(estore->GetNextEvent()){
    //    std::cout << "Event: " << neve << std::endl;
    nppachit = 0;
    nppachit_f7_0 = 0;
    nppachit_f7_1 = 0;
    isGGEarlyClose = kFALSE;
    for(int k=0;k<128;k++) ppactdc_first[k] = -1;
    for(int k=0;k<128;k++) ppactdc_second[k] = -1;
    for(int k=0;k<24;k++) rpv130_bit[k] = kFALSE;

    // start to scan raw event tree
    for(int i=0;i<rawevent->GetNumSeg();i++){
      TArtRawSegmentObject *seg = rawevent->GetSegment(i);
      int device = seg->GetDevice();
      int fp = seg->GetFP();
      int detector = seg->GetDetector();
      int module = seg->GetModule();

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
	    hf7ppac1tdc->Fill(val);
	  }
	  if(ch>=56&&ch<64) nppachit_f7_1 ++;

	}

      } // end of PPACT

      if(device != USERGR) continue;

      if(detector == 10) { //rpv130

	TArtRawDataObject *d = seg->GetData(0); //0x8000
	unsigned char mybit = d->GetVal() & 0xff;
	for(int k=0;k<8;k++){
	  if(mybit%2==1) rpv130_bit[k] = kTRUE;
	  mybit /= 2;
	}

	d = seg->GetData(1); //0x3300
	mybit = d->GetVal() & 0xff;
	for(int k=0;k<8;k++){
	  if(mybit%2==1) rpv130_bit[k+8] = kTRUE;
	  mybit /= 2;
	}

	d = seg->GetData(2); //0xf000
	mybit = d->GetVal() & 0xff;
	for(int k=0;k<8;k++){
	  if(mybit%2==1) rpv130_bit[k+16] = kTRUE;
	  mybit /= 2;
	}

	for(int k=0;k<24;k++)if(rpv130_bit[k])hrpv130bit->Fill(k);

      }

      if(rpv130_bit[15]) isGGEarlyClose = kTRUE;

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

    estore->ClearData();
    tr->Fill();
    neve ++;

  }
  TFile *fout = new TFile(Form("output/ppac/ppac_hit.%i.root",runNo),"recreate");
  hrpv130bit->Write();
  hf7ppac1tdc->Write();
  tr->Write();
  fout->Write();
  fout->Close();

  std::cout << "neve: " << neve << std::endl;

}
