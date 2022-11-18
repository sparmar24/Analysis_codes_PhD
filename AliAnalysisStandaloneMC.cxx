
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TRandom.h>
#include <TMath.h>
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDPmdTrack.h"
#include "AliESDVertex.h"
#include "AliESDVZERO.h"
#include "AliESDFMD.h"
#include "AliCentrality.h"
#include "AliAnalysisStandaloneMC.h"
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

Int_t nc=0;
Int_t secm1[22], secp1[22];
Float_t p[22], c[22];
Float_t outm[10][20];

//void amax(int eta);

ClassImp(AliAnalysisStandaloneMC)

//________________________________________________________________________
AliAnalysisStandaloneMC::AliAnalysisStandaloneMC(const char *name) 
: AliAnalysisTaskSE(name), 
  fESD(0), 
  fTreeOut(0),
  fOutputList(0),  
  fHistTotEvent(0), 
  fHistEta(0),
  fCentrality(0),
  fCentrality1(0),
  fchTrk(0),
  fpmdTrk(0),
  chTrk(0),
  pmdTrk(0),
  fHistCent(0){
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

}

//________________________________________________________________________
void AliAnalysisStandaloneMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
 
  fOutputList = new TList();

  fTreeOut = new TTree("fTreeOut", "Reconst ntuple");
  fTreeOut->Branch("fchTrk", &fchTrk,"fchTrk/I");
  fTreeOut->Branch("fpmdTrk", &fpmdTrk,"fpmdTrk/I");
  fTreeOut->Branch("fCentrality", &fCentrality,"fCentrality/I");
  fTreeOut->Branch("fCentrality1", &fCentrality1,"fCentrality1/F");

  Int_t kNBinsEvent = 10; Float_t XminEvent = 0; Float_t XmaxEvent = 10;
  Int_t kNbinsXY = 200; Float_t XminXY = -100.0; Float_t XmaxXY  = 100.0;
  fHistTotEvent = new TH1F("TotEvent","TotEvent",kNBinsEvent,XminEvent,XmaxEvent);  
  fOutputList->Add(fHistTotEvent);
  
  fHistCent = new TH1F("Cent","Cent",101, -1, 100);
  fOutputList->Add(fHistCent);

  fHistEta = new TH1F("fHistEta", "Eta", 100, 2.5, 4.0);
  fOutputList->Add(fHistEta);
}

//________________________________________________________________________
void AliAnalysisStandaloneMC::UserExec(Option_t *) 
{
  // Main loop
  Float_t MipCut1 = 432; // 6 MIP
  Float_t MipCut2 = 648; // 9 MIP
  Int_t Cls = 0,ngammal=0, ngammat=0;
  Int_t nch = 0, PhotonInc=0;
  Float_t pf[20],mf[20];

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
  /*
  mf[0]=1.04088;mf[1]=1.01775;mf[2]=1.01688;mf[3]=1.00541;mf[4]=1.00373;mf[5]=1;mf[6]=1.01624;mf[7]=1.01541;
  mf[8]=1.03513;mf[9]=1.04301;mf[10]=1.06344;mf[11]=1.0722;mf[12]=1.07441;mf[13]=1.05367;
  pf[0]=1.04586;pf[1]=1.04017;pf[2]=1.00709;pf[3]=1.12888;pf[4]=1.08158;pf[5]=1.02348;pf[6]=1.02248;
  pf[7]=1;pf[8]=1.15356;pf[9]=1.08182;pf[10]=1.04051;pf[11]=1.05372;pf[12]=1.17521;pf[13]=1.06774;
  */
  float mPVz,voc,npmdcl,npmdcl1,track1,cl1c,tpcc,voc10;;
  Int_t RunNo,eventno, nvertexCon;


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
  
  //----------------------->Centrality Selection
    AliCentrality *centrality = fESD->GetCentrality();
    Float_t nCentrality = 0,  nCentrality1 = 0, nCentrality2 = 0;
    nCentrality = centrality->GetCentralityPercentile("V0M");
    //Printf("Centrality ============================= %f",nCentrality);       
    fHistCent->Fill(nCentrality);
    fCentrality1 = nCentrality;
    fCentrality = nCentrality/10;

    //---------------------> Vertex Selection 
    const AliESDVertex *vertex = fESD->GetPrimaryVertex();
    if(!vertex) return;
    if(vertex->GetNContributors() < 1) return;
    nvertexCon = vertex->GetNContributors();
    if(vertex->GetZRes() == 0) return;
    if(vertex->GetX()<10e-5 && vertex->GetY()<10e-5 &&  vertex->GetZ()<10e-5) return;
    
    Float_t Vz = vertex->GetZRes();    
    Bool_t zVerStatus = vertex->GetStatus();
    mPVz = vertex->GetZ();
    if(TMath::Abs(mPVz) >= 10.0) return;
    
    
    //----------------->SPD                                                                                                    
    if (fESD->IsPileupFromSPD(3,0.8)) return;
    const AliESDVertex* vtx = fESD->GetPrimaryVertexSPD();
    if (!vtx || !vtx->GetStatus()) return;
    if (vtx->IsFromVertexerZ() && (vtx->GetDispersion() > 0.2 || vtx->GetZRes() > 1.25 * 0.2)) return;
    const AliMultiplicity* mult = fESD->GetMultiplicity();
    if (!mult) return;
    Double_t vz = vtx->GetZ();
    Int_t nTracklets = mult->GetNumberOfTracklets();
    //----------end spd
    
    
    npmdcl   = fESD->GetNumberOfPmdTracks();
    track1   = fESD->GetNumberOfTracks();
    nc++;
    npmdcl1  = npmdcl;
    RunNo    = fESD->GetRunNumber();
    eventno  = fESD->GetEventNumberInFile();
    

    //==================================================================                                                                                               
    //--------->>PMD<<--------------------------------------------------                                                       
    //==================================================================
    
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
      Int_t chTrk=0., pmdTrk=0.;

      for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
	AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
	if (!track) {
	  Printf("ERROR: Could not receive track %d", iTrack);
	  continue;
	}
	if(!stack->IsPhysicalPrimary(iTrack)) continue;
	Int_t mpart  = track->PdgCode();
	Float_t eta   = track->Eta();
	Float_t charge = track->Charge();

	if(eta <2.8 || eta >3.6) continue;

	fHistEta->Fill(eta);

	if(charge !=0) chTrk++;
	if(mpart==22) pmdTrk++;

	//cout<< eta << endl;

      } //mpart

      fchTrk = chTrk;
      fpmdTrk = pmdTrk;
      //cout<< fchTrk <<"     "<< fpmdTrk<< endl;

      fTreeOut->Fill();
      PostData(1, fOutputList);
      PostData(2, fTreeOut);
	
}

/*
//___________________________________________SDM
	void amax(int eta)
	{
	  Float_t f[22];
	  Int_t tbin,bin;
	  Float_t p1[22],c1[22];
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
		} //w==1                
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
		} // w>1    
	      
	    Float_t max=0.0,min=999.0,pmax=0.,cmax=0.,pmin=0.,cmin=0.;
	    Float_t ma=0,mi=0;

	    for(Int_t i =1; i<=tbin;i++)
	      {
		if(p[i] == 0 && c[i]== 0) continue;
		f[i] = p[i]/2.0/(p[i]/2.0+c[i]);
		if(f[i] > max) {
		  max=f[i];  // max fraction           
		  ma=i;    // max phi                                        
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
*/

  //________________________________________________________________________
  void AliAnalysisStandaloneMC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  fTreeOut = dynamic_cast<TTree*> (GetOutputData(2));
  if (!fTreeOut) {
    printf("ERROR: Output list not available\n");
    return;
  }

  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}
  
