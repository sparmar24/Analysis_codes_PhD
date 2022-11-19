#include "TRandom3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
//#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliCentrality.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
/*
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDPmdTrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDVZERO.h" */

#include "AliAnalysisTaskPMDSP.h"
#include "AliAODVertex.h"
#include <TParticle.h>
#include <AliLog.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>
#include <iostream>
#include <fstream>
using namespace std;
//ofstream outfile,outfile1;
//outfile.open("testout", ios::out);
//outfile1.open("testEtaphi", ios::out);


ClassImp(AliAnalysisTaskPMDSP)

AliAnalysisTaskPMDSP::AliAnalysisTaskPMDSP(const char *name)
: AliAnalysisTaskSE(name),
  fAOD(0),
  fOutputList(0),
  RunNo(0),FileNo(0),EventCounter(0),BunchCross(0),Orbit(0),Period(0),
   v0EP(0),v0aEP(0),v0cEP(0),tpcEP(0),tpcEPSub1(0),tpcEPSub2(0),
   CenV0MP(0),CenV0M(0),CenTRK(0),CenCL1(0),Vx(0),Vy(0),
   Vz(0),VxSPD(0),VySPD(0),VzSPD(0),PileSPD(0),mag(0),fcov5(0),
   ntracklet(0),v0amult(0),v0cmult(0),ntrack(0),
   zdcn1en(0),zdcp1en(0),zdcn2en(0),zdcp2en(0),
   
//  fHistTotEvent(0),
  fUseOfflineTrigger(kFALSE),
  mTree(0),

  fCentralityEstimator("V0M"){
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

}

void AliAnalysisTaskPMDSP::UserCreateOutputObjects()
{
//ofstream outfile,outfile1;
  mTree = new TTree("mTree","Reconst ntuple");
  mTree->Branch("RunNo",&RunNo,"RunNo/I");
  mTree->Branch("FileNo",&FileNo,"FileNo/I");
  mTree->Branch("EventCounter",&EventCounter,"EventCounter/I");
  mTree->Branch("BunchCross",&BunchCross,"BunchCross/I");
  mTree->Branch("Orbit",&Orbit,"Orbit/I");
  mTree->Branch("Period",&Period,"Period/I");
  mTree->Branch("mag",&mag,"mag/F");
  mTree->Branch("fcov5",&fcov5,"fcov5/F");
  mTree->Branch("zdcn1en",&zdcn1en,"zdcn1en/F");
  mTree->Branch("zdcp1en",&zdcp1en,"zdcp1en/F");
  mTree->Branch("zdcn2en",&zdcn2en,"zdcn2en/F");
  mTree->Branch("zdcp2en",&zdcp2en,"zdcp2en/F");
  mTree->Branch("v0EP",&v0EP,"v0EP/F");
  mTree->Branch("v0aEP",&v0aEP,"v0aEP/F");
  mTree->Branch("v0cEP",&v0cEP,"v0cEP/F");
  mTree->Branch("tpcEP",&tpcEP,"tpcEP/F");
  mTree->Branch("tpcEPSub1",&tpcEPSub1,"tpcEPSub1/F");
  mTree->Branch("tpcEPSub2",&tpcEPSub2,"tpcEPSub2/F");
  mTree->Branch("CenV0MP",&CenV0MP,"CenV0MP/F");
  mTree->Branch("CenV0M",&CenV0M,"CenV0M/I");
  mTree->Branch("CenTRK",&CenTRK,"CenTRK/F");
  mTree->Branch("CenCL1",&CenCL1,"CenCL1/F");
  mTree->Branch("Vx",&Vx,"Vx/F");
  mTree->Branch("Vy",&Vy,"Vy/F");
  mTree->Branch("Vz",&Vz,"Vz/F");
  mTree->Branch("VxSPD",&VxSPD,"VxSPD/F");
  mTree->Branch("VySPD",&VySPD,"VySPD/F");
  mTree->Branch("VzSPD",&VzSPD,"VzSPD/F");
  mTree->Branch("PileSPD",&PileSPD,"PileSPD/F");
  mTree->Branch("ntracklet",&ntracklet,"ntracklet/F");
  mTree->Branch("v0amult",&v0amult,"v0amult/F");
  mTree->Branch("v0cmult",&v0cmult,"v0cmult/F");
  mTree->Branch("v0a",&v0a,"v0a[32]/F");
  mTree->Branch("v0c",&v0c,"v0c[32]/F");
  mTree->Branch("ntrack",&ntrack,"ntrack/I");
  mTree->Branch("EtaPtFlag",&EtaPtFlag,"EtaPtFlag[ntrack]/I");
  mTree->Branch("PhiITSFbChi",&PhiITSFbChi,"PhiITSFbChi[ntrack]/I");
  mTree->Branch("ncldcaz",&ncldcaz,"ncldcaz[ntrack]/I");
   
//   outfile.open("testout", ios::out);
//   outfile1.open("testetaphi", ios::out);
                                                                         
  fOutputList = new TList();
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
//  fHistTotEvent = new TH1F("fHistTotEvent
//                             "Event statistics;;N_{events}",
//			   4,0.5,4.5);
//  for(Int_t i = 1; i <= 4; i++)
//    fHistTotEvent->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
//  fOutputList->Add(fHistTotEvent);
  
  
}

void AliAnalysisTaskPMDSP::UserExec(Option_t *)
{
  Int_t nmult=0;

fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;}
EventCounter++;
//  fHistTotEvent->Fill(1);
 // Input handler
  AliAODInputHandler *handler = (AliAODInputHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!handler)
    {
      cout<<"FATAL: could not get Input Handler."<<endl;
      return;
    }


Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if(!isSelected) return;
 
//    fHistTotEvent->Fill(2);

//  Int_t nVertex = fAOD->GetNumberOfVertices();
//  if( nVertex > 0 ) {
//Float_t Vx,Vy,Vz;
  AliAODVertex* vertex = (AliAODVertex*)fAOD->GetPrimaryVertex();
 AliAODVertex* vertexspd = (AliAODVertex*)fAOD->GetPrimaryVertexSPD();
  Int_t nTracksPrim = vertex->GetNContributors();
//if(nTracksPrim>1 && vertex->GetZRes() != 0)fHistTotEvent->Fill(3);
   BunchCross = fAOD->GetBunchCrossNumber();
   Orbit      =fAOD->GetOrbitNumber();
   Period     =fAOD->GetPeriodNumber();
  ULong64_t evt = BunchCross +  Orbit*3564 +Period*16777215*3564;
//  outfile<<BunchCross<<"  "<<Orbit<<"  "<<Period<<"  "<<evt<<endl;
  Vz = vertex->GetZ();
 Vx = vertex->GetX();
 Vy = vertex->GetY();
VxSPD = vertexspd->GetX();
VySPD = vertexspd->GetY();
VzSPD = vertexspd->GetZ();
 Double_t fCov[6];
    vertex->GetCovarianceMatrix(fCov);
// outfile<<"fCov[5]"<<fCov[5]<<endl;
fcov5 = fCov[5];
//if(fCov[5] = 0)continue;
if( nTracksPrim < 1 || TMath::Abs(Vz) >= 12.0) {
    return;
  }

  if(TMath::Abs(Vx) >= 0.3) {
    return ;
  }
if(TMath::Abs(Vy) >= 0.3) {
    return;
  }
if(TMath::Abs(Vz) >12.0)  {
    return;
  }
if(vertex->GetType()==AliAODVertex::kPileupSPD){
return;}

//  fHistTotEvent->Fill(4); //Analyzed events
  RunNo =fAOD->GetRunNumber();

/*
string  filename = CurrentFileName(); 
unsigned pos = filename.find("000");
unsigned pos2 = filename.find("/E");
string filename4 = filename.substr(pos,pos2-pos);
unsigned pos3 = filename.find("AOD");
unsigned pos4 = filename.find("/");
string filename5 = filename.substr(pos3,pos4-pos3);
unsigned pos5 = filename5.find("/");
unsigned pos6 = filename5.find("/A");
string filename6 = filename5.substr(pos5,pos6-pos5);
string filename7 =filename6.substr(1,10);
Int_t dub= atoi(filename4.c_str());
Int_t dub1= atoi(filename7.c_str());
//cout<<pos<<"  "<<pos2<<" "<<pos3<<" "<<pos4<<" "<<pos5<<" "<<pos6<<endl;

//cout<<filename<<"  "<<filename4<<"  "<<filename5<<endl;
//cout<<filename6<<" "<<filename7<<endl;
//cout<<RunNo<<'\t'<<dub<<'\t'<<dub1<<endl; 
 if(FileNo != dub1)
{
  FileNo =dub1;
  EventCounter=0;
// outfile<<RunNo<<'\t'<<FileNo<<'\t'<<EventCounter<<'\t'<<dub<<'\t'<<dub1<<endl; 
} 
*/
// **********
//cout<<" **** "<<Vx<<"  "<<VxSPD<<"  "<<Vy<<"  "<<VySPD<<"  "<<Vz<<"  "<<VzSPD<<endl;
// *********** Centrality
 AliCentrality *centrality = fAOD->GetCentrality();
   CenV0MP = centrality->GetCentralityPercentile("V0M");
    CenV0M  = centrality->GetCentralityClass10("V0M");
CenTRK = centrality->GetCentralityPercentile("TRK");
 CenCL1 = centrality->GetCentralityPercentile("CL1");
// outfile<<"Centrality="<<CenV0MP<<"  "<<CenV0M<<"  "<<CenTRK<<"  "<<CenCL1<<endl;
// ***********
PileSPD = fAOD->IsPileupFromSPD();
// outfile<<"PileSPD="<<PileSPD<<endl;


// outfile<<" **** "<<Vx<<"  "<<VxSPD<<"  "<<Vy<<"  "<<VySPD<<"  "<<Vz<<"  "<<VzSPD<<endl;

//   if( CenV0M >8)return;
   if( CenV0M <0)return;
// outfile1<<RunNo<<'\t'<<FileNo<<'\t'<<EventCounter<<'\t'<<dub<<'\t'<<dub1<<endl; 
// outfile<<RunNo<<'\t'<<FileNo<<'\t'<<EventCounter<<'\t'<<dub<<'\t'<<dub1<<endl; 
  if(fAOD->GetNumberOfTracks()<10){
   return;}
Float_t Qx,Qy,EPv0,EPv0a,EPv0c,EPtrk,EPZDC;
Float_t Q1x,Q1y,Q2x,Q2y;
tpcEP=999;
     TVector2 * QvecEP = fAOD->GetEventplane()->GetQVector();
    if(!QvecEP)
     {
      tpcEP=999;
     }
    else
    {
    Qx = QvecEP->X();
     Qy = QvecEP->Y();
tpcEP=TMath::ATan2(Qy,Qx)/2.; 
 }  
if(!fAOD->GetEventplane()->GetQsub1())
  {
  tpcEPSub1 =999.0;
  }
else 
{

Q1x =  fAOD->GetEventplane()->GetQsub1()->X();
     Q1y =  fAOD->GetEventplane()->GetQsub1()->Y();
tpcEPSub1=TMath::ATan2(Q1y,Q1x)/2.; 
}
if(!fAOD->GetEventplane()->GetQsub2())
  {
  tpcEPSub2 =999.0;
  }
  else
{
     Q2x =  fAOD->GetEventplane()->GetQsub2()->X();
     Q2y =  fAOD->GetEventplane()->GetQsub2()->Y();
tpcEPSub2=TMath::ATan2(Q2y,Q2x)/2.; 
}
//cout<<tpcEP<<"  "<<tpcEPSub1<<"  "<<tpcEPSub2<<endl;
//if(tpcEP <0)tpcEP +=3.141592;
 v0EP = fAOD->GetEventplane()->GetEventplane("V0",fAOD,2);
 v0aEP = fAOD->GetEventplane()->GetEventplane("V0A",fAOD,2);
 v0cEP = fAOD->GetEventplane()->GetEventplane("V0C",fAOD,2);
//if(v0EP <0)v0EP +=3.141592;
//if(v0aEP <0)v0aEP +=3.141592;
//if(v0cEP <0)v0cEP +=3.141592;

// EPtrk = fAOD->GetEventplane()->GetEventplane("TPC",fAOD,2);
// EPZDC = fAOD->GetEventplane()->GetEventplane("ZDC",fAOD,2);
//cout<<Qx<<" "<<Qy<<"  "<<"  "<<tpcEP<<"  "<<v0EP<<"  "<<v0aEP<<"  "<<v0cEP<<"  "<<EPtrk<<"  "<<EPZDC<<endl;
//outfile<<tpcEP<<'\t'<<v0EP<<endl;
 mag = fAOD->GetMagneticField();
//outfile<<"mag ="<<mag<<endl;
//****** Zero degree
  zdcn1en =  (Float_t)fAOD->GetZDCN1Energy();
   zdcp1en=  (Float_t)fAOD->GetZDCP1Energy();
   zdcn2en=  (Float_t)fAOD->GetZDCN2Energy();
  zdcp2en =  (Float_t)fAOD->GetZDCP2Energy();
//outfile<<"zdc="<<zdcn1en<<" "<<zdcp1en<<"  "<<zdcn2en<<"  "<<zdcp2en<<endl;
//******* SPD
// AliMultiplicity* mult = fAOD->GetMultiplicity();
//outfile<<"mult="<<mult<<endl;
//if (!mult) return;
 ntracklet = fAOD->GetTracklets()->GetNumberOfTracklets();
 v0amult=fAOD->GetVZEROData()->GetMTotV0A();
 v0cmult=fAOD->GetVZEROData()->GetMTotV0C();
//AliAODVZRO *v0 =fAOD->GetVZEROData();
   for (Int_t i=0; i<32; i++)
                    {
 v0c[i] = (float)(fAOD->GetVZEROData()->GetMultiplicityV0C(i));
//   cout<<i<<'\t'<<v0c[i]<<endl;
      }

   for (Int_t k=0; k<32; k++)
                    {
v0a[k] = (float)(fAOD->GetVZEROData()->GetMultiplicityV0A(k));
//   cout<<k<<'\t'<<v0a[k]<<endl;
  }

// outfile<<"nTracklets="<<ntracklet<<"  "<<v0amult<<"  "<<v0cmult<<endl;
 for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTracks);
            if (!track) {
                  printf("ERROR: Could not receive track %d\n", iTracks);
                        continue;
                            }
if(track->GetType() != AliAODTrack::kPrimary)continue;
		    float p,Pt, phi1, dcaxy,dcaz, eta;
		    double dz[2], covar[3],nclus;
		    float chi2ndf,dcax,dcay,dxy;
		    TBits cluster;
        if (!(track->TestFilterBit(128) || track->TestFilterBit(272))) continue;
 int fbittpc=0,fbithib=0;
  if(track->TestFilterBit(128))fbittpc=1;
  if(track->TestFilterBit(272))fbithib=2;
  int  NITShits = 0;
        if (track->HasPointOnITSLayer(0)) NITShits++;
        if (track->HasPointOnITSLayer(1)) NITShits++;
        if (track->HasPointOnITSLayer(2)) NITShits++;
        if (track->HasPointOnITSLayer(3)) NITShits++;
        if (track->HasPointOnITSLayer(4)) NITShits++;
        if (track->HasPointOnITSLayer(5)) NITShits++;
   double r[3];
  int flag =0;
      Bool_t dcaflag = track->GetXYZ(r);
      if(dcaflag) flag=1;

                    Pt = track->Pt();
                    eta = track->Eta();
                    phi1 = track->Phi();
                      p = track->P();
		    Short_t Charge = track->Charge();
                    chi2ndf = track->Chi2perNDF();
int		        cluster1 = track->GetTPCNcls(); 

                 dcaxy = track->DCA();
                    dcaz  = track->ZAtDCA();
//cout<<"******************"<<"dcaxy="<<dcaxy<<"dcaz="<<dcaz<<endl;

		    if(fabs(eta) > 0.8) continue;
                    if(Pt <0.15)continue;
                    if(Pt >5.0)continue;
                    if(cluster1 < 60 )continue;
 nmult++;
int  phi2 = phi1*1000;
int  eta1 = fabs(eta)*100;
int z100, xy100;
int  pt1 =Pt*100;
int chi =chi2ndf*10;
// outfile<<phi1<<"  "<<phi2<<"  "<<NITShits<<"  "<<fbithib<<"  "<<chi<<endl;
PhiITSFbChi[nmult-1] = phi2*10000+NITShits*1000+fbithib*100+chi;
// outfile<<"PhiITSFbChi[nmult-1]="<<PhiITSFbChi[nmult-1]<<endl;
if(Charge <0)PhiITSFbChi[nmult-1]=-PhiITSFbChi[nmult-1];
//outfile<<"PhiITSFbChi[nmult-1]="<<PhiITSFbChi[nmult-1]<<endl;
//outfile<<eta<<" "<<eta1<<"  "<<Pt<<"  "<<pt1<<endl;
EtaPtFlag[nmult-1] = eta1*10000+10*pt1+flag;
if(eta <0)EtaPtFlag[nmult-1]=-EtaPtFlag[nmult-1];
//outfile<<"EtaPtFlag[nmult-1]="<<EtaPtFlag[nmult-1]<<endl;
  if(flag ==1)
{
z100=fabs(dcaz)*100;
xy100=fabs(dcaxy)*100;
//cout<<dcaz<<"  "<<z100<<"  "<<dcaxy<<"  "<<xy100<<endl;
 }
else
{
          dcax = r[0] - Vx;
          dcay = r[1] - Vy;
          dcaz = r[2] - Vz;
           dcaxy = sqrt(dcax*dcax+dcay*dcay);
z100=fabs(dcaz)*100;
xy100=fabs(dcaxy)*100;
//cout<<dcaz<<"  "<<z100<<"  "<<dcaxy<<"  "<<xy100<<endl;
} 
// outfile<<"chi="<<chi<<"z100="<<z100<<"xy100="<<xy100<<endl;
            ncldcaz[nmult-1] = 1000000*cluster1+1000*z100+xy100;
//outfile<<"ncldcaz[nmult-1]="<<ncldcaz[nmult-1]<<endl;
		  }
  
   
ntrack = nmult;


mTree->Fill();
	     
  PostData(1, fOutputList);
  PostData(2, mTree);
  
}

void AliAnalysisTaskPMDSP::Terminate(Option_t *)
{
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    mTree = dynamic_cast<TTree*> (GetOutputData(2));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }

}
