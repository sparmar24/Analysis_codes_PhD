#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDPmdTrack.h"
#include "AliESDVertex.h"
#include "AliESDVZERO.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskPMDEff.h"
#include <TParticle.h>
#include <AliLog.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <string>
#include <sstream>
#include <iostream>
ClassImp(AliAnalysisTaskPMDEff)

//________________________________________________________________________
AliAnalysisTaskPMDEff::AliAnalysisTaskPMDEff(const char *name) 
: AliAnalysisTaskSE(name), 
  fESD(0), 
  fOutputList(0),  
  fHistTotEvent(0), 
  fHistXYPre(0),
  fHistXYCpv(0),
  fHistCent(0),
  fHistGammaLikeMipCut1Ncell2(0),
  fHistGammaTrueMipCut1Ncell2(0),
  fHistGammaLikeMipCut2Ncell2(0),
  fHistGammaTrueMipCut2Ncell2(0),
  fHistGammaLikeMipCut1Ncell0(0),
  fHistGammaTrueMipCut1Ncell0(0),
  fHistGammaLikeMipCut2Ncell0(0),
  fHistGammaTrueMipCut2Ncell0(0),
  fHistMultPhotonInc(0),
  fCentralityEstimator("V0M"){
  for(Int_t i=0; i<8; i++){
    fHistGammaLikeEtaBinMipCut1Ncell2[i] = 0;
    fHistGammaTrueEtaBinMipCut1Ncell2[i] = 0;
    fHistGammaLikeEtaBinMipCut2Ncell2[i] = 0;
    fHistGammaTrueEtaBinMipCut2Ncell2[i] = 0;
    fHistGammaLikeEtaBinMipCut1Ncell0[i] = 0;
    fHistGammaTrueEtaBinMipCut1Ncell0[i] = 0;
    fHistGammaLikeEtaBinMipCut2Ncell0[i] = 0;
    fHistGammaTrueEtaBinMipCut2Ncell0[i] = 0;
    fHistMultPhotonIncEtaBin[i] = 0;
    fHistXYPreEtaBin[i]=0; 
    fHistXYCpvEtaBin[i]=0;
  }
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPMDEff::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
 
  fOutputList = new TList();
  Int_t kNBinsEvent = 10; Float_t XminEvent = 0; Float_t XmaxEvent = 10;
  Int_t kNbinsXY = 200; Float_t XminXY = -100.0; Float_t XmaxXY  = 100.0;
  fHistTotEvent = new TH1F("TotEvent","TotEvent",kNBinsEvent,XminEvent,XmaxEvent);  
  fHistXYPre = new TH2F("XYPre","XYPre",100,-100,100,100,-100,100);
  fHistXYCpv = new TH2F("XYCpv","XYCpv",100,-100,100,100,-100,100);
  fHistGammaLikeMipCut1Ncell2 = new TH1F("GammaLikeMipCut1Ncell2","GammaLikeMipCut1Ncell2",5000,0,5000);
  fHistGammaTrueMipCut1Ncell2 = new TH1F("GammaTrueMipCut1Ncell2","GammaTrueMipCut1Ncell2",5000,0,5000);
  fHistGammaLikeMipCut2Ncell2 = new TH1F("GammaLikeMipCut2Ncell2","GammaLikeMipCut2Ncell2",5000,0,5000);
  fHistGammaTrueMipCut2Ncell2 = new TH1F("GammaTrueMipCut2Ncell2","GammaTrueMipCut2Ncell2",5000,0,5000);
  fHistGammaLikeMipCut1Ncell0 = new TH1F("GammaLikeMipCut1Ncell0","GammaLikeMipCut1Ncell0",5000,0,5000);
  fHistGammaTrueMipCut1Ncell0 = new TH1F("GammaTrueMipCut1Ncell0","GammaTrueMipCut1Ncell0",5000,0,5000);
  fHistGammaLikeMipCut2Ncell0 = new TH1F("GammaLikeMipCut2Ncell0","GammaLikeMipCut2Ncell0",5000,0,5000);
  fHistGammaTrueMipCut2Ncell0 = new TH1F("GammaTrueMipCut2Ncell0","GammaTrueMipCut2Ncell0",5000,0,5000);
  fHistMultPhotonInc = new TH1F("MultPhotonInc","MultPhotonInc",10000,0,10000);
  fOutputList->Add(fHistTotEvent);
  fOutputList->Add(fHistXYPre);
  fOutputList->Add(fHistXYCpv);
  fOutputList->Add(fHistGammaLikeMipCut1Ncell2);
  fOutputList->Add(fHistGammaTrueMipCut1Ncell2);
  fOutputList->Add(fHistGammaLikeMipCut2Ncell2);
  fOutputList->Add(fHistGammaTrueMipCut2Ncell2);
  fOutputList->Add(fHistGammaLikeMipCut1Ncell0);
  fOutputList->Add(fHistGammaTrueMipCut1Ncell0);
  fOutputList->Add(fHistGammaLikeMipCut2Ncell0);
  fOutputList->Add(fHistGammaTrueMipCut2Ncell0);
  fOutputList->Add(fHistMultPhotonInc);


  fHistCent = new TH1F("Cent","Cent",101, -1, 100);
  fOutputList->Add(fHistCent);
  
  Char_t name1[256], name2[256],name3[256], name4[256], name11[256], name22[256],name33[256], name44[256], name5[256], name6[256], name7[256];
  for(Int_t i=0; i<8; i++){
    sprintf(name1,"GammaLikeMipCut1Ncell2_EtaBin%d",i+1);
    sprintf(name2,"GammaTrueMipCut1Ncell2_EtaBin%d",i+1);
    sprintf(name3,"GammaLikeMipCut2Ncell2_EtaBin%d",i+1);
    sprintf(name4,"GammaTrueMipCut2Ncell2_EtaBin%d",i+1);
    
    sprintf(name11,"GammaLikeMipCut1Ncell0_EtaBin%d",i+1);
    sprintf(name22,"GammaTrueMipCut1Ncell0_EtaBin%d",i+1);
    sprintf(name33,"GammaLikeMipCut2Ncell0_EtaBin%d",i+1);
    sprintf(name44,"GammaTrueMipCut2Ncell0_EtaBin%d",i+1);
    
    sprintf(name5,"PhotonInc_EtaBin%d",i+1);
    sprintf(name6,"XYPre_EtaBin%d",i+1);
    sprintf(name7,"XYCpv_EtaBin%d",i+1);
    
    fHistGammaLikeEtaBinMipCut1Ncell2[i] = new TH1F(name1,name1,1500,0,1500);
    fHistGammaTrueEtaBinMipCut1Ncell2[i] = new TH1F(name2,name2,1500,0,1500);
    fHistGammaLikeEtaBinMipCut2Ncell2[i] = new TH1F(name3,name3,1500,0,1500);
    fHistGammaTrueEtaBinMipCut2Ncell2[i] = new TH1F(name4,name4,1500,0,1500);

    fHistGammaLikeEtaBinMipCut1Ncell0[i] = new TH1F(name11,name11,1500,0,1500);
    fHistGammaTrueEtaBinMipCut1Ncell0[i] = new TH1F(name22,name22,1500,0,1500);
    fHistGammaLikeEtaBinMipCut2Ncell0[i] = new TH1F(name33,name33,1500,0,1500);
    fHistGammaTrueEtaBinMipCut2Ncell0[i] = new TH1F(name44,name44,1500,0,1500);
    
    fHistMultPhotonIncEtaBin[i] = new TH1F(name5,name5,1500,0,1500);
    fHistXYPreEtaBin[i] = new TH2F(name6,name6,100,-100,100,100,-100,100);
    fHistXYCpvEtaBin[i] = new TH2F(name7,name7,100,-100,100,100,-100,100);
    
    fOutputList->Add(fHistGammaLikeEtaBinMipCut1Ncell2[i]);  
    fOutputList->Add(fHistGammaTrueEtaBinMipCut1Ncell2[i]);
    fOutputList->Add(fHistGammaLikeEtaBinMipCut2Ncell2[i]);  
    fOutputList->Add(fHistGammaTrueEtaBinMipCut2Ncell2[i]);

    fOutputList->Add(fHistGammaLikeEtaBinMipCut1Ncell0[i]);  
    fOutputList->Add(fHistGammaTrueEtaBinMipCut1Ncell0[i]);
    fOutputList->Add(fHistGammaLikeEtaBinMipCut2Ncell0[i]);  
    fOutputList->Add(fHistGammaTrueEtaBinMipCut2Ncell0[i]);
    
    fOutputList->Add(fHistMultPhotonIncEtaBin[i]);
    fOutputList->Add(fHistXYPreEtaBin[i]);
    fOutputList->Add(fHistXYCpvEtaBin[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPMDEff::UserExec(Option_t *) 
{
  // Main loop
  Float_t MipCut1 = 432; // 6 MIP
  Float_t MipCut2 = 648; // 9 MIP
  Float_t etacls, theta, rdist;
  Int_t GammaLikeClsMipCut1Ncell2 = 0,MultPhotonInc = 0, GammaTrueMipCut1Ncell2 =0;
  Int_t GammaLikeClsMipCut1Ncell0 = 0, GammaTrueMipCut1Ncell0 =0;
  Int_t GammaLikeClsMipCut2Ncell2 = 0, GammaTrueMipCut2Ncell2 =0;
  Int_t GammaLikeClsMipCut2Ncell0 = 0, GammaTrueMipCut2Ncell0 =0;
  Int_t Cls = 0,ngammal=0, ngammat=0;
  Int_t nch = 0, PhotonInc=0;

  Int_t GammaLikeClsEtaBinMipCut1Ncell2[8] = {0};
  Int_t GammaTrueEtaBinMipCut1Ncell2[8]={0};
  Int_t GammaLikeClsEtaBinMipCut2Ncell2[8] = {0};
  Int_t GammaTrueEtaBinMipCut2Ncell2[8]={0};

  Int_t GammaLikeClsEtaBinMipCut1Ncell0[8] = {0};
  Int_t GammaTrueEtaBinMipCut1Ncell0[8]={0};
  Int_t GammaLikeClsEtaBinMipCut2Ncell0[8] = {0};
  Int_t GammaTrueEtaBinMipCut2Ncell0[8]={0};

  Int_t PhotonIncEtaBin[8] = {0};
  
  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  fHistTotEvent->Fill(5);
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if (! isSelected) return;
  printf("There are %d tracks in this event\n", fESD->GetNumberOfTracks());
  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    nch++;
  } //track loop 
  

    AliCentrality *centrality = fESD->GetCentrality();
    Float_t nCentrality = 0, nCentrality1 =0;
    nCentrality = centrality->GetCentralityPercentile("V0M");
    //    Printf("Centrality ============================= %f",nCentrality);       
    //    if(nCentrality > 0  && nCentrality <= 5){
    //    Printf("Centrality ============================= %f",nCentrality);
    nCentrality1 = centrality->GetCentralityClass10(fCentralityEstimator.Data());
    fHistCent->Fill(nCentrality);


  const AliESDVertex *vertex = fESD->GetPrimaryVertex();
  Float_t Vz = vertex->GetZRes();    
  Bool_t zVerStatus = vertex->GetStatus();
  if(zVerStatus){
    // fHistVtxZ->Fill(Vz);
    if(TMath::Abs(Vz)<10){
      //Pmd Track Loop
      Int_t ptracks = fESD->GetNumberOfPmdTracks();
      for(Int_t kk=0;kk<ptracks;kk++){
	AliESDPmdTrack *pmdtr = fESD->GetPmdTrack(kk);
	Int_t   det   = pmdtr->GetDetector();
	Float_t clsX  = pmdtr->GetClusterX();
	Float_t clsY  = pmdtr->GetClusterY();
	Float_t clsZ  = pmdtr->GetClusterZ();
	clsZ -= Vz;
	Float_t ncell = pmdtr->GetClusterCells();
	Float_t adc   = pmdtr->GetClusterADC();
	Int_t pid = pmdtr->GetClusterTrackPid();  	 
	Int_t smn = pmdtr->GetSmn();
	if(det==1) smn += 24;
	rdist = TMath::Sqrt(clsX*clsX + clsY*clsY);
	if(clsZ!=0) theta = TMath::ATan2(rdist,clsZ);
	etacls  = -TMath::Log(TMath::Tan(0.5*theta));
	if(det==0 && adc>6*72)fHistXYPre->Fill(clsX,clsY);
	if(det==1 && adc>6*72)fHistXYCpv->Fill(clsX,clsY);
	/*	
		if(det==0){  //Gamma Like Loop Starts
		if(adc>MipCut1Ncell2 && ncell>2){ 
		if(etacls>2.3 && etacls<3.9){
		//	fHistEta->Fill(etacls);
		ngammal++;
		}
		}
		}//Gamma Like Loop Ends
	*/
	if(det==0){  //Gamma Like Loop Starts                     
	  if(adc>MipCut1 && ncell>2){
	    if(etacls>2.3 && etacls<=2.5) {GammaLikeClsEtaBinMipCut1Ncell2[0]++;
	      fHistXYPreEtaBin[0]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[0]->Fill(clsX,clsY);
	    }
	    if(etacls>2.5 && etacls<=2.7) {GammaLikeClsEtaBinMipCut1Ncell2[1]++;
	      fHistXYPreEtaBin[1]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[1]->Fill(clsX,clsY);
	    }
	    if(etacls>2.7 && etacls<=2.9) {GammaLikeClsEtaBinMipCut1Ncell2[2]++;
	      fHistXYPreEtaBin[2]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[2]->Fill(clsX,clsY);
	    }
	    if(etacls>2.9 && etacls<=3.1) {GammaLikeClsEtaBinMipCut1Ncell2[3]++;
	      fHistXYPreEtaBin[3]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[3]->Fill(clsX,clsY);
	    }
	    if(etacls>3.1 && etacls<=3.3) {GammaLikeClsEtaBinMipCut1Ncell2[4]++;
	      fHistXYPreEtaBin[4]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[4]->Fill(clsX,clsY);
	    }
	    if(etacls>3.3 && etacls<=3.5) {GammaLikeClsEtaBinMipCut1Ncell2[5]++;
	      fHistXYPreEtaBin[5]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[5]->Fill(clsX,clsY);
	    }
	    if(etacls>3.5 && etacls<=3.7) {GammaLikeClsEtaBinMipCut1Ncell2[6]++;
	      fHistXYPreEtaBin[6]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[6]->Fill(clsX,clsY);
	    }
	    if(etacls>3.7 && etacls<=3.9) {GammaLikeClsEtaBinMipCut1Ncell2[7]++;
	      fHistXYPreEtaBin[7]->Fill(clsX,clsY);
	      fHistXYCpvEtaBin[7]->Fill(clsX,clsY);
	    }
	    if(etacls>2.3 && etacls<3.9) GammaLikeClsMipCut1Ncell2++;
	    
	  }
	  
	   if(adc>MipCut1 && ncell>0){
	    if(etacls>2.3 && etacls<=2.5) GammaLikeClsEtaBinMipCut1Ncell0[0]++;
	    if(etacls>2.5 && etacls<=2.7) GammaLikeClsEtaBinMipCut1Ncell0[1]++;
	    if(etacls>2.7 && etacls<=2.9) GammaLikeClsEtaBinMipCut1Ncell0[2]++;
	    if(etacls>2.9 && etacls<=3.1) GammaLikeClsEtaBinMipCut1Ncell0[3]++;
	    if(etacls>3.1 && etacls<=3.3) GammaLikeClsEtaBinMipCut1Ncell0[4]++;
	    if(etacls>3.3 && etacls<=3.5) GammaLikeClsEtaBinMipCut1Ncell0[5]++;
	    if(etacls>3.5 && etacls<=3.7) GammaLikeClsEtaBinMipCut1Ncell0[6]++;
	    if(etacls>3.7 && etacls<=3.9) GammaLikeClsEtaBinMipCut1Ncell0[7]++;
	    if(etacls>2.3 && etacls<3.9) GammaLikeClsMipCut1Ncell0++;
	   }

	    if(adc>MipCut2 && ncell>2){
	    if(etacls>2.3 && etacls<=2.5) GammaLikeClsEtaBinMipCut2Ncell2[0]++;
	    if(etacls>2.5 && etacls<=2.7) GammaLikeClsEtaBinMipCut2Ncell2[1]++;
	    if(etacls>2.7 && etacls<=2.9) GammaLikeClsEtaBinMipCut2Ncell2[2]++;
	    if(etacls>2.9 && etacls<=3.1) GammaLikeClsEtaBinMipCut2Ncell2[3]++;
	    if(etacls>3.1 && etacls<=3.3) GammaLikeClsEtaBinMipCut2Ncell2[4]++;
	    if(etacls>3.3 && etacls<=3.5) GammaLikeClsEtaBinMipCut2Ncell2[5]++;
	    if(etacls>3.5 && etacls<=3.7) GammaLikeClsEtaBinMipCut2Ncell2[6]++;
	    if(etacls>3.7 && etacls<=3.9) GammaLikeClsEtaBinMipCut2Ncell2[7]++;
	    if(etacls>2.3 && etacls<3.9) GammaLikeClsMipCut2Ncell2++;
	   }

	    if(adc>MipCut2 && ncell>0){
	    if(etacls>2.3 && etacls<=2.5) GammaLikeClsEtaBinMipCut2Ncell0[0]++;
	    if(etacls>2.5 && etacls<=2.7) GammaLikeClsEtaBinMipCut2Ncell0[1]++;
	    if(etacls>2.7 && etacls<=2.9) GammaLikeClsEtaBinMipCut2Ncell0[2]++;
	    if(etacls>2.9 && etacls<=3.1) GammaLikeClsEtaBinMipCut2Ncell0[3]++;
	    if(etacls>3.1 && etacls<=3.3) GammaLikeClsEtaBinMipCut2Ncell0[4]++;
	    if(etacls>3.3 && etacls<=3.5) GammaLikeClsEtaBinMipCut2Ncell0[5]++;
	    if(etacls>3.5 && etacls<=3.7) GammaLikeClsEtaBinMipCut2Ncell0[6]++;
	    if(etacls>3.7 && etacls<=3.9) GammaLikeClsEtaBinMipCut2Ncell0[7]++;
	    if(etacls>2.3 && etacls<3.9) GammaLikeClsMipCut2Ncell0++;
	   }
		  /*	if(pid==22){ //Gamma True loop starts
		if(adc>MipCut1Ncell2 && ncell>2){
		if(etacls>2.3 && etacls<3.9){
		ngammat++;
		}
		}
		}// Gamma True loop ends
		  */ 
	  if(pid==22){ //Gamma True loop starts
	    if(adc>MipCut1 && ncell>2){
	      if(etacls>2.3 && etacls<=2.5) GammaTrueEtaBinMipCut1Ncell2[0]++;
	      if(etacls>2.5 && etacls<=2.7) GammaTrueEtaBinMipCut1Ncell2[1]++;
	      if(etacls>2.7 && etacls<=2.9) GammaTrueEtaBinMipCut1Ncell2[2]++;
	      if(etacls>2.9 && etacls<=3.1) GammaTrueEtaBinMipCut1Ncell2[3]++;
	      if(etacls>3.1 && etacls<=3.3) GammaTrueEtaBinMipCut1Ncell2[4]++;
	      if(etacls>3.3 && etacls<=3.5) GammaTrueEtaBinMipCut1Ncell2[5]++;
	      if(etacls>3.5 && etacls<=3.7) GammaTrueEtaBinMipCut1Ncell2[6]++;
	      if(etacls>3.7 && etacls<=3.9) GammaTrueEtaBinMipCut1Ncell2[7]++;
	      if(etacls>2.3 && etacls<3.9)  GammaTrueMipCut1Ncell2++;
	    }
	     if(adc>MipCut1 && ncell>0){
	      if(etacls>2.3 && etacls<=2.5) GammaTrueEtaBinMipCut1Ncell0[0]++;
	      if(etacls>2.5 && etacls<=2.7) GammaTrueEtaBinMipCut1Ncell0[1]++;
	      if(etacls>2.7 && etacls<=2.9) GammaTrueEtaBinMipCut1Ncell0[2]++;
	      if(etacls>2.9 && etacls<=3.1) GammaTrueEtaBinMipCut1Ncell0[3]++;
	      if(etacls>3.1 && etacls<=3.3) GammaTrueEtaBinMipCut1Ncell0[4]++;
	      if(etacls>3.3 && etacls<=3.5) GammaTrueEtaBinMipCut1Ncell0[5]++;
	      if(etacls>3.5 && etacls<=3.7) GammaTrueEtaBinMipCut1Ncell0[6]++;
	      if(etacls>3.7 && etacls<=3.9) GammaTrueEtaBinMipCut1Ncell0[7]++;
	      if(etacls>2.3 && etacls<3.9)  GammaTrueMipCut1Ncell0++;
	    }
	    if(adc>MipCut2 && ncell>2){
	      if(etacls>2.3 && etacls<=2.5) GammaTrueEtaBinMipCut2Ncell2[0]++;
	      if(etacls>2.5 && etacls<=2.7) GammaTrueEtaBinMipCut2Ncell2[1]++;
	      if(etacls>2.7 && etacls<=2.9) GammaTrueEtaBinMipCut2Ncell2[2]++;
	      if(etacls>2.9 && etacls<=3.1) GammaTrueEtaBinMipCut2Ncell2[3]++;
	      if(etacls>3.1 && etacls<=3.3) GammaTrueEtaBinMipCut2Ncell2[4]++;
	      if(etacls>3.3 && etacls<=3.5) GammaTrueEtaBinMipCut2Ncell2[5]++;
	      if(etacls>3.5 && etacls<=3.7) GammaTrueEtaBinMipCut2Ncell2[6]++;
	      if(etacls>3.7 && etacls<=3.9) GammaTrueEtaBinMipCut2Ncell2[7]++;
	      if(etacls>2.3 && etacls<3.9)  GammaTrueMipCut2Ncell2++;
	    }
	     if(adc>MipCut2 && ncell>0){
	      if(etacls>2.3 && etacls<=2.5) GammaTrueEtaBinMipCut2Ncell0[0]++;
	      if(etacls>2.5 && etacls<=2.7) GammaTrueEtaBinMipCut2Ncell0[1]++;
	      if(etacls>2.7 && etacls<=2.9) GammaTrueEtaBinMipCut2Ncell0[2]++;
	      if(etacls>2.9 && etacls<=3.1) GammaTrueEtaBinMipCut2Ncell0[3]++;
	      if(etacls>3.1 && etacls<=3.3) GammaTrueEtaBinMipCut2Ncell0[4]++;
	      if(etacls>3.3 && etacls<=3.5) GammaTrueEtaBinMipCut2Ncell0[5]++;
	      if(etacls>3.5 && etacls<=3.7) GammaTrueEtaBinMipCut2Ncell0[6]++;
	      if(etacls>3.7 && etacls<=3.9) GammaTrueEtaBinMipCut2Ncell0[7]++;
	      if(etacls>2.3 && etacls<3.9)  GammaTrueMipCut2Ncell0++;
	    }
	  }
	}
      }
    }
  }
	for(Int_t i=0; i<8; i++){
	  fHistGammaLikeEtaBinMipCut1Ncell2[i]->Fill(GammaLikeClsEtaBinMipCut1Ncell2[i]);
	  fHistGammaTrueEtaBinMipCut1Ncell2[i]->Fill(GammaTrueEtaBinMipCut1Ncell2[i]);
	  fHistGammaLikeEtaBinMipCut2Ncell2[i]->Fill(GammaLikeClsEtaBinMipCut2Ncell2[i]);
	  fHistGammaTrueEtaBinMipCut2Ncell2[i]->Fill(GammaTrueEtaBinMipCut2Ncell2[i]);
	  
	  fHistGammaLikeEtaBinMipCut1Ncell0[i]->Fill(GammaLikeClsEtaBinMipCut1Ncell0[i]);
	  fHistGammaTrueEtaBinMipCut1Ncell0[i]->Fill(GammaTrueEtaBinMipCut1Ncell0[i]);
	  fHistGammaLikeEtaBinMipCut2Ncell0[i]->Fill(GammaLikeClsEtaBinMipCut2Ncell0[i]);
	  fHistGammaTrueEtaBinMipCut2Ncell0[i]->Fill(GammaTrueEtaBinMipCut2Ncell0[i]);
	}

	/*	cout<< GammaLikeClsEtaBinMipCut1Ncell2[0]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[0]<<"   "<< PhotonIncEtaBin[0] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[1]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[1]<<"   "<< PhotonIncEtaBin[1] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[2]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[2]<<"   "<< PhotonIncEtaBin[2] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[3]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[3]<<"   "<< PhotonIncEtaBin[3] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[4]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[4]<<"   "<< PhotonIncEtaBin[4] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[5]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[5]<<"   "<< PhotonIncEtaBin[5] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[6]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[6]<<"   "<< PhotonIncEtaBin[6] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[7]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[7]<<"   "<< PhotonIncEtaBin[7] << endl;*/

	fHistGammaLikeMipCut1Ncell2->Fill(GammaLikeClsMipCut1Ncell2);
	fHistGammaTrueMipCut1Ncell2->Fill(GammaTrueMipCut1Ncell2);
	fHistGammaLikeMipCut2Ncell2->Fill(GammaLikeClsMipCut2Ncell2);
	fHistGammaTrueMipCut2Ncell2->Fill(GammaTrueMipCut2Ncell2);

	fHistGammaLikeMipCut1Ncell0->Fill(GammaLikeClsMipCut1Ncell0);
	fHistGammaTrueMipCut1Ncell0->Fill(GammaTrueMipCut1Ncell0);
	fHistGammaLikeMipCut2Ncell0->Fill(GammaLikeClsMipCut2Ncell0);
	fHistGammaTrueMipCut2Ncell0->Fill(GammaTrueMipCut2Ncell0);

    
      // }
      //  }

      
      //reading MC info_____________________.///// Generator Level 
      AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>
	(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
	Printf("ERROR: Could not retrieve MC event handler");
    return;
      }
      AliMCEvent* mcEvent = eventHandler->MCEvent();
      if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
      }
      AliStack* stack = mcEvent->Stack();
      if (!stack)
	{
	  AliDebug(AliLog::kError, "Stack not available");
	  return;
	}
      for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
	AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
	if (!track) {
	  Printf("ERROR: Could not receive track %d", iTrack);
	  continue;
	}
	if(!stack->IsPhysicalPrimary(iTrack)) continue;
	Int_t mpart  = track->PdgCode();
	Float_t eta   = track->Eta();
	if(mpart==22){
	  if(eta>2.3 && eta<=2.5) PhotonIncEtaBin[0]++;
	  if(eta>2.5 && eta<=2.7) PhotonIncEtaBin[1]++;
	  if(eta>2.7 && eta<=2.9) PhotonIncEtaBin[2]++;
	  if(eta>2.9 && eta<=3.1) PhotonIncEtaBin[3]++;
	  if(eta>3.1 && eta<=3.3) PhotonIncEtaBin[4]++;
	  if(eta>3.3 && eta<=3.5) PhotonIncEtaBin[5]++;
	  if(eta>3.5 && eta<=3.7) PhotonIncEtaBin[6]++;
	  if(eta>3.7 && eta<=3.9) PhotonIncEtaBin[7]++;
	  if(eta>2.3 && eta<3.9)  PhotonInc++;
	}
      }//mpart
	for(Int_t i=0; i<8; i++){
	  //      if(PhotonIncEtaBin[i] !=0){
	  fHistMultPhotonIncEtaBin[i]->Fill(PhotonIncEtaBin[i]); 
	}//i loop                                                              
	fHistMultPhotonInc->Fill(PhotonInc);
	/*
	cout<< GammaLikeClsEtaBinMipCut1Ncell2[0]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[0]<<"   "<< PhotonIncEtaBin[0] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[1]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[1]<<"   "<< PhotonIncEtaBin[1] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[2]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[2]<<"   "<< PhotonIncEtaBin[2] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[3]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[3]<<"   "<< PhotonIncEtaBin[3] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[4]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[4]<<"   "<< PhotonIncEtaBin[4] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[5]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[5]<<"   "<< PhotonIncEtaBin[5] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[6]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[6]<<"   "<< PhotonIncEtaBin[6] << endl;
        cout<< GammaLikeClsEtaBinMipCut1Ncell2[7]<<"   "<<GammaTrueEtaBinMipCut1Ncell2[7]<<"   "<< PhotonIncEtaBin[7] << endl;
	cout<< "=========================================================================================================" << endl;
	cout<< GammaLikeClsMipCut1Ncell2 <<"    "<< GammaTrueMipCut2Ncell2<< "   "<< PhotonInc << endl;
	*/
	//   }//centrality
  PostData(1, fOutputList);     
}
  //________________________________________________________________________
  void AliAnalysisTaskPMDEff::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}
  
