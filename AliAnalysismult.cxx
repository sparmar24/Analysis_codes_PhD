#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include <TRandom.h>

#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TAxis.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVZERO.h"
#include "AliESDFMD.h"
# include "AliMultiplicity.h"
# include "AliESDVertex.h"

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>

#include "AliESDInputHandler.h"
#include "AliESDPmdTrack.h"
#include "AliCentrality.h"
#include "AliESDtrackCuts.h"
//#include "AliPhysicsSelection.h"
#include "AliAnalysismult.h"
#include <iostream>
ofstream of;
ofstream of1;
ofstream of2,of3,of4;
Int_t nc=0;
Int_t secm1[22],secm5[22],secm5f[22],secm4[22];
Int_t secp1[22],secp5[22],secp5f[22],secp4[22];
Float_t p[22],c[22];
Float_t outm[10][20];


//void amax(float &max,float &min,int &ma,int &mi,int eta,int bin);
void amax(int eta);
// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysismult) 
//________________________________________________________________________
//AliAnalysismult::AliAnalysismult(const char *name)  : AliAnalysisTaskSE(name), mcEvent(0), fOutputList(0),Vertex_dist(),pmdxy(),fHistEventCount(0),fHistPmdFmdmult(0),fHistFmdmult(0),fHistFmdmul2i(0),fHistFmdPhi2i(0),fHistFmdEta2i(0),fCentralityEstimator("V0M"),
//fHistCentrality(0),pmdxy1(0),fHistPmdPhi(0)
AliAnalysismult::AliAnalysismult(const char *name)  : AliAnalysisTaskSE(name),fESD(0), fTreeOut(0),fOutputList(0),fCentralityEstimator("V0M")
{
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    fHistPhiBin[iBin-1] = NULL;
    fHistPhiFull[iBin-1] = NULL;
    fHistPhiFMD[iBin-1] = NULL;
    fHistPhiFMD253[iBin-1] = NULL;
    fHistFMDmult[iBin-1] = NULL;
  }
  // Constructor
  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysismult::UserCreateOutputObjects()
{
  fOutputList = new TList;
  // Create output TTree
  fTreeOut = new TTree("fTreeOut", "Reconst ntuple");
  fTreeOut->Branch("outh", &outh,"outh[23]/I");  
  fTreeOut->Branch("out", &out,"out[8][10]/I");  
  fTreeOut->Branch("oute",&oute,"oute[20]/I");
  
  of.open("test", ios::out);
  of1.open("test1", ios::out);
  of2.open("test2", ios::out);
  of3.open("test3", ios::out);
  of4.open("test4", ios::out);
  // Create histograms
  TString histName;
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    histName = "fHistPhiBin"; histName += "Centrality"; histName += iBin;
    fHistPhiBin[iBin-1] = new TH1F(histName.Data(),"Phi Distribution",20,0.,360.);
    fOutputList->Add(fHistPhiBin[iBin-1]);
  } 
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    histName = "fHistPhiFull"; histName += "Centrality"; histName += iBin;
    fHistPhiFull[iBin-1] = new TH1F(histName.Data(),"Phi Distribution",20,0.,360.);
    fOutputList->Add(fHistPhiFull[iBin-1]);
  } 
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    histName = "fHistPhiFMD"; histName += "Centrality"; histName += iBin;
    fHistPhiFMD[iBin-1] = new TH1F(histName.Data(),"Phi Distribution",20,0.,20);
    fOutputList->Add(fHistPhiFMD[iBin-1]);
  } 
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    histName = "fHistPhi253"; histName += "Centrality"; histName += iBin;
    fHistPhiFMD253[iBin-1] = new TH1F(histName.Data(),"Phi Distribution",20,0.,20.);
    fOutputList->Add(fHistPhiFMD253[iBin-1]);
  } 
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    histName = "fHistFMDmult"; histName += "Centrality"; histName += iBin;
    fHistFMDmult[iBin-1] = new TH1F(histName.Data(),"mult Distribution",50,0.,2.);
    fOutputList->Add(fHistFMDmult[iBin-1]);
  } 
}

//________________________________________________________________________
void AliAnalysismult::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  //  fTreeOut->Branch("outm", &outm,"outm[10][20]/F");  
  // Post output data.
  
  Float_t pf[20],mf[20];
  /* fpmd(0)=0.877096;fpmd(1)=1.16159;fpmd(2)=1.03782;fpmd(3)=1.48217;
     fpmd(4)=1.044609;fpmd(5)=0.87845;fpmd(6)=0.937379;fpmd(7)=0.927648;
     fpmd(8)=1.03132;fpmd(9)=0.869793;fpmd(10)=0.876925;fpmd(11)=1.24582;
     fpmd(12)=1.02576;fpmd(13)=0.922184;
     ffmd(0)=0.997453;ffmd(1)=0.966253;ffmd(2)=0.961018;ffmd(3)=0.945919;
     ffmd(4)=0.95285;ffmd(5)=0.95593;ffmd(6)=0.979913;ffmd(7)=0.993723;
     ffmd(8)=1.02923;ffmd(9)=1.04131;ffmd(10)=1.06827;ffmd(11)=1.07069;
     ffmd(12)=1.04789;ffmd(13)=1.0141; */
  /*
    pf[0]=0.877096;pf[1]=1.16159;pf[2]=1.03782;pf[3]=1.48217;
    pf[4]=1.04609;pf[5]=0.87845;pf[6]=0.937379;pf[7]=0.927648;
    pf[8]=1.03132;pf[9]=0.869793;pf[10]=0.876925;pf[11]=1.24582;
    pf[12]=1.022576;pf[13]=0.922184;
    mf[0]=0.997453;mf[1]=0.966253;mf[2]=0.961018;mf[3]=0.945919;
    mf[4]=0.95285;mf[5]=0.95593;mf[6]=0.979913;mf[7]=0.993723;
    mf[8]=1.02923;mf[9]=1.04131;mf[10]=1.06827;mf[11]=1.07069;
    mf[12]=1.04789;mf[13]=1.0141;
  */
  
  //Weight factors for MC sample of root file 138534 & 138225 ,,,,By normalizing it to the highest bin content
  
  mf[0]=1.04088;mf[1]=1.01775;mf[2]=1.01688;mf[3]=1.00541;mf[4]=1.00373;mf[5]=1;mf[6]=1.01624;mf[7]=1.01541;
  mf[8]=1.03513;mf[9]=1.04301;mf[10]=1.06344;mf[11]=1.0722;mf[12]=1.07441;mf[13]=1.05367;
  
  pf[0]=1.04586;pf[1]=1.04017;pf[2]=1.00709;pf[3]=1.12888;pf[4]=1.08158;pf[5]=1.02348;pf[6]=1.02248;
  pf[7]=1;pf[8]=1.15356;pf[9]=1.08182;pf[10]=1.04051;pf[11]=1.05372;pf[12]=1.17521;pf[13]=1.06774;
  
  Int_t  mc,pmd1,pmdg2,multpmd144,multfmd144;
  
  Int_t temp[9][3];
  float pi=3.14156;
  float mPVz,voc;
  int num,nvertexCon;
  Int_t RunNo,eventno,timestamp,npmdcl,npmdcl1,track1,cl1c,tpcc,voc10;
  Int_t OrbitNo,BunchCros,PeriodNo;
  
  AliVEvent* event = InputEvent();
  const char * filename1 ="/Users/madanaggarwal/Downloads/";
  const char * filename2;
  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(event);
  
  if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return;
  const AliESDVertex *PrimaryVertex = fESD->GetPrimaryVertex();
  if(!PrimaryVertex) return;
  if(PrimaryVertex->GetNContributors() < 1)return;
  nvertexCon = PrimaryVertex->GetNContributors();
  if(PrimaryVertex->GetZRes() == 0)return; 
  if(PrimaryVertex->GetX()<10e-5 && PrimaryVertex->GetY()<10e-5 &&  PrimaryVertex->GetZ()<10e-5) return;
  //                if(PrimaryVertex->GetX()<3.0) return;
  //                if(PrimaryVertex->GetY()<3.0) return;
  mPVz = PrimaryVertex->GetZ(); 
  if(fabs(mPVz)>10.0)return;

  //  SPD
  if (fESD->IsPileupFromSPD(3,0.8)) return;
  const AliESDVertex* vtx = fESD->GetPrimaryVertexSPD();
  if (!vtx || !vtx->GetStatus()) return;
  if (vtx->IsFromVertexerZ() &&
      (vtx->GetDispersion() > 0.2 || vtx->GetZRes() > 1.25 * 0.2))
    return;
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  if (!mult) return;
  Double_t vz = vtx->GetZ();
  Int_t nTracklets = mult->GetNumberOfTracklets();
  //cout<<nTracklets<<endl;
  //for (Int_t i = 0; i < nTracklets; i++)
  //fMult->Fill(mult->GetEta(i), vz);
  //cout<<mult->GetEta(i);
  //  fHistPt->Fill(mult->GetEta(i)); 
  // end spd
  
  
  AliCentrality *centrality = fESD->GetCentrality();
  voc = centrality->GetCentralityPercentile("V0M");
  Int_t nCentrality = 0,nCentrality1 = 0,nCentrality2 = 0,nCentrality3 = 0,nCentrality5 = 0,nCentrality4 = 0;
  nCentrality = centrality->GetCentralityClass10(fCentralityEstimator.Data());
  nCentrality1 = centrality->GetCentralityClass5(fCentralityEstimator.Data());
  nCentrality2 = centrality->GetCentralityClass10("CL1");
  nCentrality3 = centrality->GetCentralityClass5("CL1");
  nCentrality4 = centrality->GetCentralityClass10("TRK");
  nCentrality5 = centrality->GetCentralityClass5("TRK");
  //cout<<"*****************    "<<nCentrality<<'\t'<<nCentrality<<endl;
  //  if(nCentrality > 5)return;
  for(Int_t i=0;i<10;i++)
    for(Int_t k=0;k<21;k++){
      outm[i][k]=0.0;}
  if( fESD->GetNumberOfTracks()> 0){printf("There are %d tracks in this event\n", fESD->GetNumberOfTracks());
  }
  //     if(fESD->GetNumberOfTracks() <200)return;
  npmdcl=fESD->GetNumberOfPmdTracks();
  track1=fESD->GetNumberOfTracks();
  nc++;
  //      if(nc >1 )return;
  //      if(track1 != 234 || track1 != 121)return;
  npmdcl1=npmdcl;
  RunNo =fESD->GetRunNumber();
  eventno =fESD->GetEventNumberInFile();
  //     if(eventno !=21)return;
  OrbitNo =fESD->GetOrbitNumber();
  BunchCros =fESD-> GetBunchCrossNumber();
  PeriodNo  =fESD->GetPeriodNumber();
  timestamp = fESD->GetTimeStamp();
  //cout<<RunNo<<'\t'<<eventno<<'\t'<<OrbitNo<<'\t'<<BunchCros<<'\t'<<PeriodNo<<'\t'<<timestamp<<endl;
  
  // Process FMD event summary data
  AliESDFMD *fmd = fESD->GetFMDData();
  Int_t fmdmul=0,fmdmul2i=0,fmdmult=0;
  for(Int_t i =0;i<22;i++){
    secm1[i]=0;
    secm5[i]=0;
    secm4[i]=0;
    secm5f[i]=0;
    secp1[i]=0;
    secp5[i]=0;
    secp4[i]=0;
    secp5f[i]=0;
  } 
  
  //  if (!fESDEvent) return kFALSE;
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* rng = rings; *rng != '\0'; rng++) {
      UShort_t nsec = (*rng == 'I' ?  20 :  40);
      UShort_t nstr = (*rng == 'I' ? 512 : 256);
      for (UShort_t sec = 0; sec < nsec; sec++) {
        for (UShort_t str = 0; str < nstr; str++) {
          Float_t mult = fmd->Multiplicity(det,*rng,sec,str);
	  if (mult == AliESDFMD::kInvalidMult || mult >20) continue;
          Float_t eta  = fmd->Eta(det,*rng,sec,str);
          if (!fmd->IsAngleCorrected())
            mult *= TMath::Abs(TMath::Cos(2.*TMath::ATan(TMath::Exp(-eta))));
	  if(mult >0.3)fmdmult++;
          if(det == 2 && *rng == 'I' ){
	    if(nCentrality<=4)fHistFMDmult[nCentrality]->Fill(mult);
            if(mult>0.3){
	      fmdmul2i++;
	      if(nCentrality<=4){
		if(eta >=2.8 && eta <=3.6)fHistPhiFMD[nCentrality]->Fill(sec+0.5);
		if(eta >=2.28 && eta <2.8)fHistPhiFMD253[nCentrality]->Fill(sec+0.5);
	      }          
	      if(eta >=2.8 && eta <= 3.6)secm1[sec+1]++;
	      if(eta >=2.8 && eta <= 3.6 && mult > 1.5)secm1[sec+1]++;
	      if(eta >=2.45 && eta <= 2.95)secm5[sec+1]++;
	      if(eta >=2.45 && eta <= 2.85)secm4[sec+1]++;
	      if(eta >=2.95 && eta <= 3.45)secm5f[sec+1]++; 
	      
	      if(eta >=2.28 && eta <= 3.68)//of4<<sec<<endl;}
	        if(eta >=2.8 && eta <= 3.6){
		  of1<<nCentrality<<'\t'<<fmdmul<<'\t'<<det<<'\t'<<sec<<'\t'<<str<<'\t'<<eta<<'\t'<<mult<<endl;
		  
		}
	    }
	  }
	  //          if (!ProcessESD(det, *rng, sec, str, eta, mult)) continue;
        }
      }
    }
  }
  Int_t iphi;
  Float_t phi;
  mc=0,pmd1=0,pmdg2=0;
  while (npmdcl--) {
    AliESDPmdTrack *pmdtr = fESD->GetPmdTrack(npmdcl);
    
    Int_t   det   = pmdtr->GetDetector();
    Float_t clsX  = pmdtr->GetClusterX();
    Float_t clsY  = pmdtr->GetClusterY();
    Float_t clsZ  = pmdtr->GetClusterZ();
    Float_t ncell = pmdtr->GetClusterCells();
    Float_t adc   = pmdtr->GetClusterADC();
    Float_t pid   = pmdtr->GetClusterPID();
    Float_t sigx = pmdtr->GetClusterSigmaX();
    Float_t sigy = pmdtr->GetClusterSigmaY();
    Float_t  phi1=TMath::ATan(clsY/clsX);
    // pmdxy->Fill(clsX,clsY);
    // if(phi <0)phi=-phi;
    if(clsX >=0 && clsY >=0)phi=phi1;
    if(clsX <0 && clsY<0)phi=pi+phi1;
    if(clsX <0 && clsY>=0)phi=pi+phi1;
    if(clsY<0 && clsX>=0)phi=2.0*pi+phi1;
    if(det == 1)continue;
    of3<<det<<'\t'<<clsX<<'\t'<<clsY<<'\t'<<phi<<'\t'<<phi1<<'\t'<<ncell<<'\t'<<adc<<'\t'<<pid<<'\t'<<sigx<<'\t'<<sigy<<endl;
    if(ncell == 1 )pmd1 +=1;
    if(adc <432 && ncell >1)pmdg2 +=1;
    if(adc >432 && ncell >1){
      //if(adc >0 && ncell >0){
      //if(adc >216 ){
      //if(adc <216 && ncell >1)pmdg2 +=1;
      //if(adc >216 && ncell >1){
      Float_t xy2 =TMath::Sqrt(clsX*clsX+clsY*clsY);
      Float_t theta = TMath::ATan(xy2/clsZ);
      Float_t rap = -TMath::Log(TMath::Tan(theta/2.0));
      if(rap >=2.28 && rap <= 3.68)mc++;
      phi=phi*180.0/3.14;
      if(rap >=2.8 && rap <= 3.6){
	of4<<clsX<<'\t'<<clsY<<'\t'<<rap<<endl;
	if(adc <432)of2<<mc<<'\t'<<clsZ<<'\t'<<clsX<<'\t'<<clsY<<'\t'<<theta<<'\t'<<rap<<'\t'<<phi<<'\t'<<adc<<'\t'<<ncell<<endl;
	if(adc >432 && ncell>1)of<<mc<<'\t'<<clsZ<<'\t'<<clsX<<'\t'<<clsY<<'\t'<<theta<<'\t'<<rap<<'\t'<<phi<<'\t'<<ncell<<endl;  }
      iphi = phi/18.0+1;
      //     of<<mc<<'\t'<<nCentrality<<'\t'<<clsZ<<'\t'<<clsX<<'\t'<<clsY<<'\t'<<theta<<'\t'<<rap<<'\t'<<iphi<<endl;  
      if(iphi <=20){
	if(nCentrality<=4){
	  fHistPhiFull[nCentrality]->Fill(phi);
	  if(rap>=2.8 && rap <=3.6)fHistPhiBin[nCentrality]->Fill(phi);
	}
	if(rap >=2.8 && rap <= 3.6)secp1[iphi]++;
	if(rap >=2.45 && rap <= 2.95)secp5[iphi]++;
	if(rap >=2.45 && rap <= 2.85)secp4[iphi]++;
	if(rap >=2.95 && rap <= 3.45)secp5f[iphi]++;
      }
      //     of<<mc<<'\t'<<clsZ<<'\t'<<clsX<<'\t'<<clsY<<'\t'<<theta<<'\t'<<rap<<endl;  
    }
  }
  int ma,mi,eta,bin;
  eta=1;
  //  cout<<max1<<'\t'<<ma<<'\t'<<secm1[ma]<<'\t'<<secp1[ma]<<'\t'<<min1<<'\t'<<mi<<'\t'<<secm1[mi]<<'\t'<<secp1[mi]<<endl;
  //     cout<<RunNo<<'\t'<<OrbitNo<<'\t'<<eventno<<'\t'<<timestamp<<'\t'<<npmdcl1<<'\t'<<mc<<'\t'<<track1<<endl;
  Int_t csec=0,k,i;
  for( k =1;k<13;k++)
    {
      secp1[k] *=pf[k-1];
      secm1[k] *=mf[k-1];
      p[k+2]=secp1[k];
      c[k+2]=secm1[k];
    }
  secp1[19] *=pf[12];
  secp1[20] *=pf[13];
  secm1[19] *=mf[12];
  secm1[20] *=mf[13];  
  p[2]=secp1[20];
  c[2]=secm1[20];
  p[1]=secp1[19];
  c[1]=secm1[19];
  multpmd144=0;
  multfmd144=0;
  for( k =1;k<15;k++)
    {
      multpmd144 +=p[k];
      multfmd144 +=c[k];
    }  
  /*
    Float_t ratio;
    for(Int_t i=0;i<20;i++){
    if(c[i]>0){
    ratio=(Float_t)p[i]/(Float_t)c[i];}
    of3<<i<<'\t'<<p[i]<<'\t'<<c[i]<<'\t'<<ratio<<endl;
  */
  //  p[i] =p[i]*1000;
  //   of4<<p[i]<<endl;
  //   p[i]=p[i]+c[i];
  //   of4<<p[i]<<endl;
  //   oute[i]=p[i];
  // }
  //   of4<<oute[i]<<endl;
  /*   oute[0]=p[0];
       oute[1]=p[1];
       oute[2]=p[2];
       oute[3]=p[3];
       oute[4]=p[4];
       oute[5]=p[5];
       oute[6]=p[6];
       oute[7]=p[7];
       oute[8]=p[8];
       oute[9]=p[9];
       oute[10]=p10];
       /*   oute[11]=p[11];
       oute[12]=p[12];
       oute[13]=p[13];
       oute[14]=p[14];
       oute[15]=p[15];
       oute[16]=p[16];
       oute[17]=p[17];
       oute[18]=p[18];
       oute[19]=p[19];
   //   oute[20]=p[20];
   */
  //   of4<<p[1]<<'\t'<<p[2]<<endl;

  
  Int_t rand = gRandom->Uniform(1,14);
  temp[1][1]=p[rand];
  temp[1][2]=c[rand];
  rand = gRandom->Uniform(1,13);
  temp[2][1]=0;temp[2][2]=0;
  for(i=rand;i<rand+2;i++){
    temp[2][1] +=p[i];
  temp[2][2] +=c[i];}
  rand = gRandom->Uniform(1,12);
  temp[3][1]=0;temp[3][2]=0;
  for(i=rand;i<rand+3;i++){
    temp[3][1] +=p[i];
    temp[3][2] +=c[i];}
  rand = gRandom->Uniform(1,11);
  temp[4][1]=0;temp[4][2]=0;
  for(i=rand;i<rand+4;i++){
    temp[4][1] +=p[i];
    temp[4][2] +=c[i];}
  rand = gRandom->Uniform(1,10);
  temp[5][1]=0;temp[5][2]=0;
  for(i=rand;i<rand+5;i++){
    temp[5][1] +=p[i];
    temp[5][2] +=c[i];}
  rand = gRandom->Uniform(1,9);
  temp[6][1]=0;temp[6][2]=0;
  for(i=rand;i<rand+6;i++){
    temp[6][1] +=p[i];
    temp[6][2] +=c[i];}
  rand = gRandom->Uniform(1,8);
  temp[7][1]=0;temp[7][2]=0;
  for(i=rand;i<rand+7;i++){
    temp[7][1] +=p[i];
    temp[7][2] +=c[i];}
  rand = gRandom->Uniform(1,7);
  temp[8][1]=0;temp[8][2]=0;
  for(i=rand;i<rand+8;i++){
    temp[8][1] +=p[i];
    temp[8][2] +=c[i];}
  amax(1);
  for(Int_t w =1;w<9;w++)
    {
      out[w-1][0] = outm[w][1]*1000+outm[w][3];
      out[w-1][1] = outm[w][4];
      out[w-1][2] = outm[w][5]*1000+outm[w][7];
      out[w-1][3] = outm[w][8];
      out[w-1][4] = outm[w][11]*1000+outm[w][13];
      out[w-1][5] = outm[w][14];
      out[w-1][6] = outm[w][15]*1000+outm[w][17];
      out[w-1][7] = outm[w][18];
      out[w-1][8] = temp[w][1];
      out[w-1][9] = temp[w][2];
    } 
  outh[0] = RunNo;
  outh[1] = eventno;
  outh[2] = nvertexCon;
  outh[3] = multpmd144;
  outh[4] = multfmd144;
  outh[5] = mPVz*1000.0;
  outh[6] = track1;
  outh[7] = npmdcl1;
  outh[8] = nCentrality;
  outh[9] = nCentrality1;
  outh[10] = nCentrality2;
  outh[11] = nCentrality3;
  outh[12] = nCentrality4;
  outh[13] = nCentrality5;
  outh[14] = mc;
  outh[15] = pmdg2;
  outh[16] = pmd1;
  outh[17] = fmdmult;
  outh[18]= fmdmul2i; 
  outh[19]= nTracklets;
  //       outh[20]= dub;
  //       outh[21]= dub1;
  //       outh[22]= dub2;
  
  //of<<voc<<'\t'<<nCentrality<<endl;
  for(Int_t w =1;w<6;w++){
    //  cout<<w<<'\t'<<outm[w][1]<<'\t'<<outm[w][2]<<'\t'<<outm[w][3]<<'\t'<<outm[w][4]<<'\t'<<outm[w][5]<<'\t'<<outm[w][6]<<'\t'<<outm[w][7]<<'\t'<<outm[w][8]<<endl;
    //cout<<w<<'\t'<<outm[w][11]<<'\t'<<outm[w][12]<<'\t'<<outm[w][13]<<'\t'<<outm[w][14]<<'\t'<<outm[w][15]<<'\t'<<outm[w][16]<<'\t'<<outm[w][17]<<'\t'<<outm[w][18]<<endl;
  }
  //  of<<w<<'\t'<<max1<<'\t'<<ma<<'\t'<<p[ma]<<'\t'<<c[ma]<<'\t'<<min1<<'\t'<<mi<<'\t'<<p[mi]<<'\t'<<c[mi]<<endl;
  for(Int_t i=1;i<21;i++){
    csec++;
    //   of<<csec<<'\t'<<secp1[i]<<'\t'<<secm1[i]<<endl;
  }
  //  fEvent=new Event(Float_t outm[10][20]);
  fTreeOut->Fill();
  PostData(1, fTreeOut);
  PostData(2, fOutputList);
}      
//void amax(float &max,float &min,int &ma,int &mi,int eta,int bin)
void amax(int eta)
{
  Float_t f[22];
  Int_t tbin,bin;
  Float_t p1[22],c1[22];
//  for(Int_t i =0;i<10;i++)
//  for(Int_t k =0;k<21;k++)
//  {
//  outm[i][k]=0.0;}
  
  for(Int_t l=1;l<3;l++){
    
    for(Int_t k =1;k<13;k++)
      {
	p1[k+2]=secp1[k];
	c1[k+2]=secm1[k]; 
      }
    p1[2]=secp1[20];
    c1[2]=secm1[20];
    p1[1]=secp1[19];
    c1[1]=secm1[19];
    if(l ==2){
      for(Int_t k =1;k<11;k++)
	{
	  c1[k+2]=secm1[11-k];
	}
      c1[1]=secm1[12];
      c1[2]=secm1[11];
      c1[13]=secm1[20];
      c1[14]=secm1[19];
    }  
    for(Int_t w =1;w<=8;w++){
      if(w == 1)
	{
	  tbin=14;
	  for(bin =1;bin<15;bin++)
	    {
	      p[bin]=p1[bin];
	      c[bin]=c1[bin];
	    }
	}
      if(w>1) 
	{
	  tbin=1;
	  for(bin =1;bin<15;bin++)
	    {
	      p[bin]=0;
	      c[bin]=0;
	    }
	  for (bin =1;bin<=w;bin++)
	    {
	      p[1] +=p1[bin];
	      c[1] +=c1[bin];
	    }
	  for(Int_t j=2;j<=14-w+1;j++)
	    {
	      tbin++;
	      p[j]=p[j-1]-p1[j-1]+p1[j+w-1];
	      c[j]=c[j-1]-c1[j-1]+c1[j+w-1];
	    }
	}
      
      
      Float_t max=0.0,min=999.0,pmax,cmax,pmin,cmin;
      Float_t ma=0,mi=0;
      
      for(Int_t i =1; i<=tbin;i++)
	{
	  if(p[i] == 0 && c[i]== 0) continue;
	  //       f[i] = secp1[i]/2.0/(secp1[i]/2.0+secm1[i]/1.5);
	  f[i] = p[i]/2.0/(p[i]/2.0+c[i]);
	  if(f[i] > max) {
	    max=f[i];
	    ma=i;
	    pmax=p[i];
	    cmax=c[i];}
	  if(f[i] <min){
	    min=f[i];
	    mi=i;
	    pmin=p[i];
	    cmin=c[i];}
	  
	}
      if(l == 1){
	outm[w][1]=ma;
	outm[w][2]=max;
	outm[w][3]=pmax;
	outm[w][4]=cmax;
	outm[w][5]=mi;
	outm[w][6]=min;
	outm[w][7]=pmin;
	outm[w][8]=cmin;
      }
      if(l == 2){
	outm[w][11]=ma;
	outm[w][12]=max;
	outm[w][13]=pmax;
	outm[w][14]=cmax;
	outm[w][15]=mi;
	outm[w][16]=min;
	outm[w][17]=pmin;
	outm[w][18]=cmin;
      }
    }
  }
  return;
} 
//________________________________________________________________________
void AliAnalysismult::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(2));
  fTreeOut = dynamic_cast<TTree*> (GetOutputData(1));
  if (!fTreeOut) {
    printf("ERROR: Output list not available\n");
    return;
  }
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
  
} 
