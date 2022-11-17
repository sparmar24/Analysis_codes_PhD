
//+++++++++++Event-by-event charge separation in pb-Pb @ 2.76 TeV
//(for rec MC)  ++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

#include <TRandom.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliCentrality.h"
#include "AliEventPoolManager.h"
#include "AliEventplane.h"

#include "AliAODVertex.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include <AliLog.h>
#include <AliStack.h>

#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>

#include "TProof.h"
#include "AliPID.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliVHeader.h"
#include "TFormula.h"
#include "AliStack.h"
#include <iostream>
#include <cerrno>
#include <memory>
#include <fstream>
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"

#include "AliAnalysisTaskEbyEChargeSep.h"
using namespace std;

ofstream outfile;

ClassImp(AliAnalysisTaskEbyEChargeSep)
ClassImp(AliMixedEvtTracks)

//-------------------------->Defualt constructor
AliAnalysisTaskEbyEChargeSep::AliAnalysisTaskEbyEChargeSep():
AliAnalysisTaskSE(),
fAOD(0),
fOutputList(0),
farrayMC(0x0),
fVtrxZ(0),
fEvtCentrality(0),
fEvtCentOrMult(0),
fCentEstimator("V0M"),
fEvtCentCL1(0),
fEvtCentTRK(0),
fV0AMult(0),
fV0CMult(0),
fV0ACMult(0),
fBunchCrossNo(0),
fOrbitNo(0),
fPeriodNo(0),
fEventCounter(0),
IsEvtSctPlotReq(kFALSE),
fCompareFileNo(0),
fTreeEventCounter(0),
fRunNumPlotReq(138275),
fFileNumPlotReq(0273),
fEvtNumPlotReq(41),
fAODProdNumber(0),
fCAODNumber(0),
fRunNumber(0),
fhTotEvent(0x0),
fhTotTrack(0x0),
fhTotPart(0x0),
fmcArray(0x0),
fDBEvtTree(0),
fBit(128),
IsDBMethodOn(kTRUE),
fReadMC(kFALSE),
fMCKine(kFALSE),
fDBWindow(90),
fEvtabsZvtx(7.),
fEventZvtx(0),
fTrkPtMin(0.2),
fTrkPtMax(5.0),
fTrkEtaMin(-0.8),
fTrkEtaMax(0.8),
fTrkabsDCAxy(3.0),
fTrkabsDCAz(3.0),
fTrkTPCclus(80),
fhEffValues(0x0),
fMixedEvent(kFALSE),
fFlowQnVectorMgr(),
flowQnVectorTask(),
fPhiRecoDgr(0),
fPhiKineDgr(0),
fIsEffCorr(0),
flMultiAllChrg(0),
flMultPosChrg(0),
flMultNegChrg(0),
V0AEvtAngleRad_latest(0),
V0CEvtAngleRad_latest(0),
TPCEvtAngleRad_latest(0),
fExpectedCorrectionPass("latest"),
fAlternativeCorrectionPass("latest")
{
    //Default Constructor
}

//--------------------------> Named constructor
AliAnalysisTaskEbyEChargeSep::AliAnalysisTaskEbyEChargeSep(const char *name):
AliAnalysisTaskSE(name),
fAOD(0),
fOutputList(0),
farrayMC(0x0),
fVtrxZ(0),
fEvtCentrality(0),
fEvtCentOrMult(0),
fCentEstimator("V0M"),
fEvtCentCL1(0),
fEvtCentTRK(0),
fV0AMult(0),
fV0CMult(0),
fV0ACMult(0),
fBunchCrossNo(0),
fOrbitNo(0),
fPeriodNo(0),
IsEvtSctPlotReq(kFALSE),
fCompareFileNo(0),
fEventCounter(0),
fTreeEventCounter(0),
fRunNumPlotReq(138275),
fFileNumPlotReq(0273),
fEvtNumPlotReq(41),
fAODProdNumber(0),
fCAODNumber(0),
fRunNumber(0),
fhTotEvent(0x0),
fhTotTrack(0x0),
fhTotPart(0x0),
fmcArray(0x0),
fDBEvtTree(0),
fBit(128),
IsDBMethodOn(kTRUE),
fReadMC(kFALSE),
fMCKine(kFALSE),
fDBWindow(90),
fEvtabsZvtx(7.),
fEventZvtx(0),
fTrkPtMin(0.2),
fTrkPtMax(5.0),
fTrkEtaMin(-0.8),
fTrkEtaMax(0.8),
fTrkabsDCAxy(3.0),
fTrkabsDCAz(3.0),
fTrkTPCclus(80),
fhEffValues(0x0),
fMixedEvent(kFALSE),
fFlowQnVectorMgr(),
flowQnVectorTask(),
fPhiRecoDgr(0),
fPhiKineDgr(0),
fIsEffCorr(0),
flMultiAllChrg(0),
flMultPosChrg(0),
flMultNegChrg(0),
V0AEvtAngleRad_latest(0),
V0CEvtAngleRad_latest(0),
TPCEvtAngleRad_latest(0),
fExpectedCorrectionPass("latest"),
fAlternativeCorrectionPass("latest")
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
}

//--------------------------> Destructor
AliAnalysisTaskEbyEChargeSep::~AliAnalysisTaskEbyEChargeSep()
{
    Info("~AliAnalysisTaskEbyEChargeSep","Calling Destructor");
    if(fhEffValues){delete fhEffValues; fhEffValues=0;}
}


//--------------------------> User Outputs
void AliAnalysisTaskEbyEChargeSep::UserCreateOutputObjects()
{
    
  outfile.open("testout",ios::out);

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    TString PartOrTrack = "";
    if(fReadMC){
        PartOrTrack = "McReco";
        if(fMCKine)PartOrTrack = "McPart";
    }else PartOrTrack = "TrkReco";
    
    TString AnaType = "";
    if(fMixedEvent)AnaType = "ME";
    else AnaType = "SE";
    
    Float_t Pi = TMath::Pi();
    fDBEvtTree = new TTree(Form("fDBEvtTree_%s",AnaType.Data()),"Store DB Event Infos");
    fDBEvtTree->Branch("fRunNumber",&fRunNumber,"fRunNumber/I");
    fDBEvtTree->Branch("fVtrxZ",     &fVtrxZ,    "fVtrxZ/F");
    fDBEvtTree->Branch("fEvtCentrality",&fEvtCentrality,"fEvtCentrality/F"); 
    fDBEvtTree->Branch("fEvtCentOrMult",&fEvtCentOrMult,"fEvtCentOrMult/F");
    fDBEvtTree->Branch("fEvtCentCL1",&fEvtCentCL1,"fEvtCentCL1/F");
    fDBEvtTree->Branch("fEvtCentTRK",&fEvtCentTRK,"fEvtCentTRK/F");
    fDBEvtTree->Branch("flMultiAllChrg", &flMultiAllChrg, "flMultiAllChrg/I");
    fDBEvtTree->Branch("flMultPosChrg",  &flMultPosChrg,  "flMultPosChrg/I");
    fDBEvtTree->Branch("flMultNegChrg",  &flMultNegChrg,  "flMultNegChrg/I");


    fhTotEvent = new TH1I("fhTotEvent", "Event Stats",4,-0.5,3.5);
    fhTotEvent->GetXaxis()->SetBinLabel(1,"Total");
    fhTotEvent->GetXaxis()->SetBinLabel(2,"OffTrig");
    fhTotEvent->GetXaxis()->SetBinLabel(3,"VtxSel");
    fhTotEvent->GetXaxis()->SetBinLabel(4,"Selected");
    fOutputList->Add(fhTotEvent);
    
    fhTotTrack = new TH1I("fhTotTrack", "Track Stats",8,-0.5,7.5);
    fhTotTrack->GetXaxis()->SetBinLabel(1,"All");
    fhTotTrack->GetXaxis()->SetBinLabel(2,"FbitPass");
    fhTotTrack->GetXaxis()->SetBinLabel(3,"NegLPass");
    fhTotTrack->GetXaxis()->SetBinLabel(4,"PtSel");
    fhTotTrack->GetXaxis()->SetBinLabel(5,"EtaSel");
    fhTotTrack->GetXaxis()->SetBinLabel(6,"TPCClusSel");
    fhTotTrack->GetXaxis()->SetBinLabel(7,"DCAClusSel");
    fhTotTrack->GetXaxis()->SetBinLabel(8,"Accecpted");
    if(!fMCKine)fOutputList->Add(fhTotTrack);
    
    fhTotPart = new TH1I("fhTotPart", "Generated Part Stats",4,-0.5,3.5);
    fhTotPart->GetXaxis()->SetBinLabel(1,"AllPart");
    fhTotPart->GetXaxis()->SetBinLabel(2,"SelecPt");
    fhTotPart->GetXaxis()->SetBinLabel(3,"SelecEta");
    fhTotPart->GetXaxis()->SetBinLabel(4,"SelectedF");
    if(fMCKine)fOutputList->Add(fhTotPart);

    //Ztvx Distribution //240
    TH1F* fhEvtZvtx  = new TH1F("fhEvtZvtx","fhEvtZvtx",600, -30., 30.);
    fhEvtZvtx->Sumw2();
    fOutputList->Add(fhEvtZvtx);

    //Multiplicity Distribution
    TH1F* fhEvtCentrality = new TH1F("fhEvtCentrality","fhCentrality",12,-2,10);
    TH1F* fhEvtCentralityPer = new TH1F("fhEvtCentralityPer","fhEvtCentralityPer",110,-2,108);        
    TH2F* fhMultVsCentAllChrg = new TH2F(Form("fhMultVsCentAllChrg_%s", AnaType.Data()),"fhMultVsCentAllChrg",Int_t(5e2),0,Int_t(2e3), 14, -2, 12);
    TH2F* fhMultVsCentPosChrg = new TH2F(Form("fhMultVsCentPosChrg_%s", AnaType.Data()),"fhMultVsCentPosChrg",Int_t(5e2),0,Int_t(2e3), 14, -2, 12);
    TH2F* fhMultVsCentNegChrg = new TH2F(Form("fhMultVsCentNegChrg_%s", AnaType.Data()),"fhMultVsCentNegChrg",Int_t(5e2),0,Int_t(2e3), 14, -2, 12);
    fhMultVsCentAllChrg->Sumw2();
    fhMultVsCentPosChrg->Sumw2();
    fhMultVsCentNegChrg->Sumw2();
    fhEvtCentrality->Sumw2();
    fhEvtCentralityPer->Sumw2();

    fOutputList->Add(fhEvtCentrality);
    fOutputList->Add(fhEvtCentralityPer);
    fOutputList->Add(fhMultVsCentAllChrg);
    fOutputList->Add(fhMultVsCentPosChrg);
    fOutputList->Add(fhMultVsCentNegChrg);

    TH1F* fhPt = new TH1F(Form("fhPt_%s", AnaType.Data()), "Trk Pt", 500, 0, 5.0);
    TH1F* fhEta = new TH1F(Form("fhEta_%s", AnaType.Data()), "Trk Eta", 100, -1, 1);
    TH1F* fhPhi = new TH1F(Form("fhPhi_%s", AnaType.Data()), "Trk Phi", 360, 0, 2*Pi);
    TH1F* fhPhid = new TH1F(Form("fhPhid_%s", AnaType.Data()), "Trk Phid", 361, 0, 361);

    fhPt->Sumw2();
    fhEta->Sumw2();
    fhPhi->Sumw2();
    fhPhid->Sumw2();
    fOutputList->Add(fhPt);
    fOutputList->Add(fhEta);
    fOutputList->Add(fhPhi);
    fOutputList->Add(fhPhid);

    TH2F* fhTrkDCAxyz = new TH2F(Form("fhTrkDCAxyz_%s", AnaType.Data()), "DCAxyVsDCAz; ", 500, -5.0, 5.0, 500, -5.0, 5.0);
    fhTrkDCAxyz->Sumw2();
    fOutputList->Add(fhTrkDCAxyz);

    if(!fMCKine){
      TH1F* fhTrkDCAxy = new TH1F(Form("fhTrkDCAxy_%s", AnaType.Data()), "Trk DCA-xy", 500, -5.0, 5.0);
      TH1F* fhTrkDCAz = new TH1F(Form("fhTrkDCAz_%s", AnaType.Data()), "Trk DCA_-z", 500, -5.0, 5.0);
      TH1F* fhTrkClusterTPC = new TH1F(Form("fhTrkClusterTPC_%s", AnaType.Data()), "Trk TPC Cluster", 200, 0, 200);
      TH1F* fhTrkChi2NDF = new TH1F(Form("fhTrkChi2NDF_%s", AnaType.Data()), "Trk Chi2NDF", 700, 0, 7.0);
      TH2F* fhTrkChrgVsDCAxy = new TH2F(Form("fhTrkChrgVsDCAxy_%s", AnaType.Data()), "DCAxyVsCharge; Trk DCA-xy; Charge;", 500, -5.0, 5.0, 10, -5, 5);
      TH2F* fhTrkChrgVsDCAz  = new TH2F(Form("fhTrkChrgVsDCAz_%s", AnaType.Data()), "DCAzVsCharge; Trk DCA-z; Charge;", 500, -5.0, 5.0, 10, -5, 5);
      fhTrkDCAxy->Sumw2();
      fhTrkDCAz->Sumw2();
      fhTrkClusterTPC->Sumw2();
      fhTrkChi2NDF->Sumw2();
      fhTrkChrgVsDCAxy->Sumw2();
      fhTrkChrgVsDCAz->Sumw2();
      fOutputList->Add(fhTrkDCAxy);
      fOutputList->Add(fhTrkDCAz);
      fOutputList->Add(fhTrkClusterTPC);
      fOutputList->Add(fhTrkChi2NDF);
      fOutputList->Add(fhTrkChrgVsDCAxy);
      fOutputList->Add(fhTrkChrgVsDCAz);
    }

    TH2F* fhChrgVsPt    = new TH2F(Form("fhChrgVsPt%s_%s",PartOrTrack.Data(), AnaType.Data()), "PtVsCharge; Trk p_{T}; Charge;", 50, 0, 5.0, 10, -5, 5);
    TH2F* fhChrgVsPhi   = new TH2F(Form("fhChrgVsPhi%s_%s",PartOrTrack.Data(), AnaType.Data()), "PhiVsCharge; Trk #varphi; Charge;", 360, 0, 2*Pi, 10, -5, 5);
    TH2F* fhChrgVsPhiD   = new TH2F(Form("fhChrgVsPhiD%s_%s",PartOrTrack.Data(), AnaType.Data()), "PhiDVsCharge; Trk #varphi; Charge;", 360, 0, 360, 10, -5, 5);
    TH2F* fhChrgVsEta   = new TH2F(Form("fhChrgVsEta%s_%s",PartOrTrack.Data(), AnaType.Data()), "EtaVsCharge; Trk #eta; Charge;", 100, -1., 1., 10, -5, 5);
    fhChrgVsPt->Sumw2();
    fhChrgVsPhi->Sumw2();
    fhChrgVsPhiD->Sumw2();
    fhChrgVsEta->Sumw2();
    fOutputList->Add(fhChrgVsPt);
    fOutputList->Add(fhChrgVsPhi);
    fOutputList->Add(fhChrgVsPhiD);
    fOutputList->Add(fhChrgVsEta);

    TH1F* fhMultV0ADet  = new TH1F("fhMultV0ADet","fhMultV0ADet",Int_t(5e2),0,Int_t(10e3));
    TH1F* fhMultV0CDet  = new TH1F("fhMultV0CDet","fhMultV0CDet",Int_t(5e2),0,Int_t(10e3));
    TH1F* fhMultV0ACDet = new TH1F("fhMultV0ACDet","fhMultV0ACDet",Int_t(5e2),0,Int_t(10e3));

    TH3F* fhPtVsChrgVsCent = new TH3F(Form("fhPtVsChrgVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()), Form("fhPtVsChrgVsCent%s",PartOrTrack.Data()), 600, 0, 6, 10, -5, 5, 12, -1, 11);
    TH3F* fhEtaVsChrgVsCent = new TH3F(Form("fhEtaVsChrgVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()), Form("fhEtaVsChrgVsCent%s",PartOrTrack.Data()), 100, -1,1, 10, -5, 5, 12, -1, 11);
    TH3F* fhPhiVsChrgVsCent = new TH3F(Form("fhPhiVsChrgVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()), Form("fhPhiVsChrgVsCent%s",PartOrTrack.Data()), 360, 0, 2*Pi, 10, -5, 5, 12, -1, 11);

    fhMultV0ADet->Sumw2();
    fhMultV0CDet->Sumw2();
    fhMultV0ACDet->Sumw2();
    fhPtVsChrgVsCent->Sumw2();
    fhEtaVsChrgVsCent->Sumw2();
    fhPhiVsChrgVsCent->Sumw2();

    fOutputList->Add(fhMultV0ADet);
    fOutputList->Add(fhMultV0CDet);
    fOutputList->Add(fhMultV0ACDet);
    fOutputList->Add(fhPtVsChrgVsCent);
    fOutputList->Add(fhEtaVsChrgVsCent);
    fOutputList->Add(fhPhiVsChrgVsCent);

    //_____________________ Qn Correction framework                                                                       
    flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    if (flowQnVectorTask != NULL) {
      fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    }
    else {
      AliFatal("Flow Qn vector corrections framework needed but it is not present. ABORTING!!!");
    }

    PostData(1, fOutputList);
    PostData(2, fDBEvtTree);
}

//--------------------------> User Exec
void AliAnalysisTaskEbyEChargeSep::UserExec(Option_t *){
    
    TString PartOrTrack = "";
    if(fReadMC){
        PartOrTrack = "McReco";
        if(fMCKine)PartOrTrack = "McPart";
    }else PartOrTrack = "TrkReco";
    
    TString AnaType = "";
    if(fMixedEvent)AnaType = "ME";
    else AnaType = "SE";   
    
    //Start of analysis
    fAOD = (AliAODEvent*)dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD && AODEvent() && IsStandardAOD()) {
        fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
    } else if(!fAOD)  {
        printf("AOD not found!\n");
        return;
    }
    
    fhTotEvent->Fill(0);
    fEventCounter++;
    
    AliAODInputHandler *handler = (AliAODInputHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!handler)return;
    fhTotEvent->Fill(1);
    
    if(fMCKine){
        fmcArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if (!fmcArray){
            AliError("Could not find Monte-Carlo in AOD");
            return;
        }
        AliAODMCHeader *mcHeader=NULL;
        mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!mcHeader) {
            AliError("Could not find MC Header in AOD");
            return;
        }
    }

    Bool_t EventSel = IsEventSelectionPassed();
    if(!EventSel)return;
    fhTotEvent->Fill(3);

    //-------------->variables stored in tree 
    fPeriodNo     = fAOD->GetPeriodNumber();
    fOrbitNo      = fAOD->GetOrbitNumber();
    fBunchCrossNo = fAOD->GetBunchCrossNumber();
    //cout<< fPeriodNo<<"    "<< fOrbitNo<<"    "<< fBunchCrossNo<< endl;
    fCAODNumber = GetFileNumberOfCurrentAOD();

    if(fCompareFileNo != fCAODNumber){
      fCompareFileNo = fCAODNumber;
      fTreeEventCounter = 0;
    }
    fTreeEventCounter++;
    
    if(IsEvtSctPlotReq){
      if(fRunNumPlotReq != fAOD->GetRunNumber())return;
      if(fFileNumPlotReq  != GetFileNumberOfCurrentAOD())return;
      if(fEvtNumPlotReq != fTreeEventCounter)return;
      // Going Back from here                                                               
    }    
    
    //----------------|| Get corrected Event plane V0A, V0C and TPC -----||
    int myHarmonic = 2;

    const AliQnCorrectionsQnVector *myQnVectorV0_latest;
    double myEventPlaneV0_latest = 0.0;
    myQnVectorV0_latest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","latest","latest");
    if (myQnVectorV0_latest != NULL) {
      myEventPlaneV0_latest = myQnVectorV0_latest->EventPlane(myHarmonic);
      V0AEvtAngleRad_latest = myEventPlaneV0_latest;
    }
    //----------------V0C
    const AliQnCorrectionsQnVector *myQnVectorV0C_latest;
    double myEventPlaneV0C_latest = 0.0;
    myQnVectorV0C_latest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","latest","latest");
    if (myQnVectorV0C_latest != NULL) {
      myEventPlaneV0C_latest = myQnVectorV0C_latest->EventPlane(myHarmonic);
      V0CEvtAngleRad_latest = myEventPlaneV0C_latest;
    }

    //---------->>>TPC Histograms                                                                                       
    const AliQnCorrectionsQnVector *myQnVectorTPC_latest;
    double myEventPlaneTPC_latest = 0.0;
    myQnVectorTPC_latest = fFlowQnVectorMgr->GetDetectorQnVector("TPC","latest","latest");
    if (myQnVectorTPC_latest != NULL) {
      myEventPlaneTPC_latest = myQnVectorTPC_latest->EventPlane(myHarmonic);
      TPCEvtAngleRad_latest = myEventPlaneTPC_latest;
    }
   
    //Filtered TrackLoop:
    Int_t lMultiAllChrg, lMultPosChrg, lMultNegChrg;
    Int_t NTrkPosEta, NTrkNegEta, TorPPhiRadiBin;
    Float_t TorPPhiRad = 0.;
    Int_t nSEMETracks = 0;
    Float_t PtWeightedCos2Phi=0.,PtWeightedSin2Phi=0.;
    Float_t TrkAryAllChrgPhi[Int_t(5e3)], TrkAryPosChrgPhi[Int_t(5e3)], TrkAryNegChrgPhi[Int_t(5e3)];
    Float_t TrkAryPosEtaPhi[Int_t(5e3)], TrkAryNegEtaPhi[Int_t(5e3)];
    Float_t NTrkAryPosChrgPhiRadBin[361], NTrkAryNegChrgPhiRadBin[361], rpp[361];
    Float_t qx,qy,q2x,q2y,qxp,qyp,q2xp,q2yp,q4x,q4y,q4xp,q4yp;
    Float_t qxn,qyn,q2xn,q2yn,q4yn,q4xn,fddb;
    Float_t cos2phi,sin2phi,cos4phi,sin4phi;
    Float_t cosch, sinch;
    Float_t phitr[5000];

    //mix charge
    Int_t chrm[5000], phibin =0, win;
    Float_t phi;

	for(Int_t iPhiDegree=0; iPhiDegree<=360; iPhiDegree++){
	  NTrkAryPosChrgPhiRadBin[iPhiDegree]=0.;
	  NTrkAryNegChrgPhiRadBin[iPhiDegree]=0.;
        }

	win=0;
	lMultiAllChrg = 0, lMultPosChrg = 0, lMultNegChrg = 0;
              
	  Int_t nTrackOrPart = 0; 
	  if(fMCKine)nTrackOrPart = fmcArray->GetEntriesFast(); //kine                               
	  else nTrackOrPart = fAOD->GetNumberOfTracks();        // Data || reco

	  for(int iTorPStr=0; iTorPStr < nTrackOrPart; iTorPStr++){  // trk loop
                Float_t TorPP = 0, TorPPt =0., TorPEta = 0., TorPPhi=0, TorPPhiDeg = 0;
                Short_t TorPChrg =0;
	                  
                if(fMCKine){ // mc kine
		  AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fmcArray->At(iTorPStr));
                    if (!mcPart){
                        AliError(Form("Failed casting particle from MC array!, Skipping particle for %d", iTorPStr));
                        continue;
                    }
		    if(!IsGeneratedPartPassed(mcPart))continue;
		    if(!mcPart->IsPhysicalPrimary())continue;

                    TorPP = mcPart->P();
                    TorPPt = mcPart->Pt();
                    TorPPhi = mcPart->Phi();
                    TorPEta = mcPart->Eta();
                    TorPChrg = mcPart->Charge();
                    TorPPhiDeg = TorPPhi*180.0/TMath::Pi();
                }else{  // data or reco

		  AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTorPStr);
                    if (!track){
                        AliError(Form("Failed casting particle from track array!, Skipping particle for %d", iTorPStr));
                        continue;
                    }
		    if(!IsTrackSelectionPassed(track))continue; 
		    //cout << "track-Pt() >>" << track->Pt() <<",   Charge >> " <<  track->Charge() << endl;
                    TorPP = track->P();
                    TorPPt = track->Pt();
                    TorPPhi = track->Phi();
                    TorPEta = track->Eta();
                    TorPChrg = track->Charge();
                    TorPPhiDeg = TorPPhi*180.0/TMath::Pi();
		
		if(TorPChrg==0)continue;
		//cout<< " Charge >> " << TorPChrg <<endl;
		if(!fMixedEvent){
		  ((TH1F*)fOutputList->FindObject(Form("fhPt_%s", AnaType.Data())))->Fill(TorPPt);
		  ((TH1F*)fOutputList->FindObject(Form("fhEta_%s", AnaType.Data())))->Fill(TorPEta);
		  ((TH1F*)fOutputList->FindObject(Form("fhPhi_%s", AnaType.Data())))->Fill(TorPPhi);
		  ((TH1F*)fOutputList->FindObject(Form("fhPhid_%s", AnaType.Data())))->Fill(TorPPhiDeg);
		  ((TH1F*)fOutputList->FindObject(Form("fhTrkDCAxy_%s", AnaType.Data())))->Fill(track->DCA());
		  ((TH1F*)fOutputList->FindObject(Form("fhTrkDCAz_%s", AnaType.Data())))->Fill(track->ZAtDCA());
		  ((TH1F*)fOutputList->FindObject(Form("fhTrkClusterTPC_%s", AnaType.Data())))->Fill(track->GetTPCNcls());
		  ((TH1F*)fOutputList->FindObject(Form("fhTrkChi2NDF_%s", AnaType.Data())))->Fill(track->Chi2perNDF());
		  ((TH2F*)fOutputList->FindObject(Form("fhTrkChrgVsDCAxy_%s", AnaType.Data())))->Fill(track->DCA(), TorPChrg+0.5);
		  ((TH2F*)fOutputList->FindObject(Form("fhTrkChrgVsDCAz_%s", AnaType.Data())))->Fill(track->ZAtDCA(), TorPChrg+0.5);
		  ((TH2F*)fOutputList->FindObject(Form("fhTrkDCAxyz_%s", AnaType.Data())))->Fill(track->DCA(),track->ZAtDCA()); 
		  fhTotTrack->Fill(7);
                }
		
		//	((TH3F*)fOutputList->FindObject(Form("fhTorPProp%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPt,TorPEta,TorPPhi);
		((TH2F*)fOutputList->FindObject(Form("fhChrgVsPt%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPt, TorPChrg+0.5);
		((TH2F*)fOutputList->FindObject(Form("fhChrgVsPhi%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPhi, TorPChrg+0.5);
		((TH2F*)fOutputList->FindObject(Form("fhChrgVsPhiD%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPhiDeg, TorPChrg+0.5);
		((TH2F*)fOutputList->FindObject(Form("fhChrgVsEta%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPEta, TorPChrg+0.5);
		((TH3F*)fOutputList->FindObject(Form("fhPtVsChrgVsCent%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPt, TorPChrg+0.5, fEvtCentrality+0.5);
		((TH3F*)fOutputList->FindObject(Form("fhEtaVsChrgVsCent%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPEta, TorPChrg+0.5, fEvtCentrality+0.5);
		((TH3F*)fOutputList->FindObject(Form("fhPhiVsChrgVsCent%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPhi, TorPChrg+0.5, fEvtCentrality+0.5);

		}

                lMultiAllChrg++;
		phitr[lMultiAllChrg]=TorPPhi;
		chrm[lMultiAllChrg]=TorPChrg;
		TrkAryAllChrgPhi[lMultiAllChrg]=TorPPhi;

                if(TorPChrg > 0){
		  lMultPosChrg++;
		  TrkAryPosChrgPhi[lMultPosChrg]=0.0;
		  TrkAryPosChrgPhi[lMultPosChrg]=TorPPhi;
                }else if(TorPChrg < 0){
		  lMultNegChrg++;
		  TrkAryNegChrgPhi[lMultNegChrg]=0.0;
		  TrkAryNegChrgPhi[lMultNegChrg]=TorPPhi;
                }
		
		if(IsEvtSctPlotReq){
		  Float_t ScatXValue = 200.0*TorPPt*cos(TorPPhi)/TorPP;
		  Float_t ScatYValue =200.0*TorPPt*sin(TorPPhi)/TorPP;
		  //cout<< "helooo = "<<"    "<< ScatXValue<< "     "<< ScatYValue << endl;
		}

                TorPPhiRad = 0.;
                TorPPhiRad = TorPPhi*180.0/TMath::Pi(); //Radian to degree                        
                TorPPhiRadiBin = Int_t(TorPPhiRad + 1);//Bin 1 to 360 degree                      
                if(TorPPhiRadiBin > 360)TorPPhiRadiBin = 360; //Protection                          

                if(TorPChrg>0)NTrkAryPosChrgPhiRadBin[TorPPhiRadiBin]++;
                else if(TorPChrg<0)NTrkAryNegChrgPhiRadBin[TorPPhiRadiBin]++;
	  
		//========<< Apply efficiency Correction:                                     
                Int_t iPt = int(TorPPt*10);
                if(fIsEffCorr && !fMCKine){
		  Double_t fEffPt = 1;
		  fEffPt = GetTrackWeight(iPt);
		  if(fEffPt >0) {
		    if(TorPChrg>0)NTrkAryPosChrgPhiRadBin[TorPPhiRadiBin] = NTrkAryPosChrgPhiRadBin[TorPPhiRadiBin]+(1/fEffPt)-1;
		    else if(TorPChrg<0)NTrkAryNegChrgPhiRadBin[TorPPhiRadiBin] = NTrkAryNegChrgPhiRadBin[TorPPhiRadiBin]+(1/fEffPt)-1;
		  }
                }
	  } //trk loop  
    
	flMultiAllChrg = lMultiAllChrg;
        flMultPosChrg = lMultPosChrg;
        flMultNegChrg = lMultNegChrg;
	//cout<< lMultiAllChrg <<"    "<<lMultPosChrg<<"    "<<lMultNegChrg<< endl;

	((TH2F*)fOutputList->FindObject(Form("fhMultVsCentAllChrg_%s", AnaType.Data())))->Fill(lMultiAllChrg, fEvtCentrality+0.5);
        ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentPosChrg_%s", AnaType.Data())))->Fill(lMultPosChrg, fEvtCentrality+0.5);
        ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentNegChrg_%s", AnaType.Data())))->Fill(lMultNegChrg, fEvtCentrality+0.5);


	if(IsEvtSctPlotReq)return; // Going back if second run requested 


	 fDBEvtTree->Fill();
	 PostData(1, fOutputList);
	 PostData(2, fDBEvtTree);
	
}


//--------------------------> Get phi efficiency values                          
Double_t AliAnalysisTaskEbyEChargeSep::GetTrackWeight(Int_t PtBin){
  if(!fhEffValues)return 1;
  if(fhEffValues->IsBinUnderflow(PtBin)||fhEffValues->IsBinOverflow(PtBin))return 1.;
  return fhEffValues->GetBinContent(PtBin);
}


//--------------------------> Get which AOD was used   
Int_t AliAnalysisTaskEbyEChargeSep::GetFileNumberOfCurrentAOD(){

  fRunNumber =fAOD->GetRunNumber();
  TString  CrFileName = CurrentFileName() ;
  TString  FileNameType = "AliAOD.root";
  unsigned PosI = CrFileName.Index(fAODProdNumber)+1+fAODProdNumber.Sizeof()-1 + 1; //index=x-1 and sizeof = x+1  
  unsigned PosF = (CrFileName.Sizeof()-1)-(PosI)-(FileNameType.Sizeof()-1);
  TString  AODNumberStr(CrFileName(PosI-1,PosF));
  //cout<<CrFileName<<"     "<< FileNameType<<"     "<<PosI<<"     "<<PosF<<"      "<<AODNumberStr <<endl;
  return Int_t(AODNumberStr.Atoi());
}

//--------------------------> Selection of Events
Bool_t AliAnalysisTaskEbyEChargeSep::IsEventSelectionPassed(){
    
    //Trigger selection
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if(!isSelected) return kFALSE;
    else fhTotEvent->Fill(2);
    
    //Vertex related selection
    Int_t nVertex = fAOD->GetNumberOfVertices();
    if(nVertex <=0)return kFALSE;
    
    AliAODVertex* Vtx = (AliAODVertex*)fAOD->GetPrimaryVertex();
    
    Int_t nTracksPrim = Vtx->GetNContributors();
    if( nTracksPrim < 1)return kFALSE;
    
    Float_t VtrxX = Vtx->GetX();
    if(TMath::Abs(VtrxX) >= 0.5)return kFALSE; 
    
    Float_t VtrxY = Vtx->GetY();
    if(TMath::Abs(VtrxY) >= 0.5)return kFALSE;
    
    Float_t VtrxZ = Vtx->GetZ();
    //    if(TMath::Abs(VtrxZ) > fEvtabsZvtx)return kFALSE;
    
    fVtrxZ = VtrxZ;
    fEventZvtx = Vtx->GetZ(); //for event Pool
    
    if(Vtx->GetType()==AliAODVertex::kPileupSPD)return kFALSE;
    
    //VZERO-AMPLITUDE
    AliAODVZERO *fAODV0 = fAOD->GetVZEROData();
    
    //Centrality Selection
    AliCentrality *EvtCent = fAOD->GetCentrality();
    fEvtCentrality = EvtCent->GetCentralityClass10("V0M");
    //if(fEvtCentrality < 0 || fEvtCentrality > 8)return kFALSE;
    if(fEvtCentrality < 0)return kFALSE;

    fEvtCentOrMult = Float_t(EvtCent->GetCentralityPercentile("V0M"));
    fEvtCentCL1 = Float_t(EvtCent->GetCentralityPercentile("CL1"));
    fEvtCentTRK = Float_t(EvtCent->GetCentralityPercentile("TRK"));

    if(fEvtCentOrMult < 0)return kFALSE;
    if(fEvtCentCL1 < 0)return kFALSE;
    if(fEvtCentTRK < 0)return kFALSE;
    fV0AMult = Float_t(fAODV0->GetMTotV0A());
    fV0CMult = Float_t(fAODV0->GetMTotV0C());
    fV0ACMult = Float_t(fAODV0->GetMTotV0A()+fAODV0->GetMTotV0C());

    //QA Plots after above cuts are applied     
    ((TH1F*)fOutputList->FindObject("fhMultV0ADet"))->Fill(fAODV0->GetMTotV0A());
    ((TH1F*)fOutputList->FindObject("fhMultV0CDet"))->Fill(fAODV0->GetMTotV0C());
    ((TH1F*)fOutputList->FindObject("fhMultV0ACDet"))->Fill(fAODV0->GetMTotV0A() + fAODV0->GetMTotV0C());
    ((TH1F*)fOutputList->FindObject("fhEvtZvtx"))->Fill(VtrxZ);
    ((TH1F*)fOutputList->FindObject("fhEvtCentrality"))->Fill(fEvtCentrality+0.5);
    ((TH1F*)fOutputList->FindObject("fhEvtCentralityPer"))->Fill(Float_t(EvtCent->GetCentralityPercentile("V0M")));
    
    return kTRUE;
}


//--------------------------> Selection of Tracks
Bool_t AliAnalysisTaskEbyEChargeSep::IsTrackSelectionPassed(AliAODTrack *fTrack){
    
  //Filterbit
  if(!fTrack->TestFilterBit(128)) return kFALSE;
  else fhTotTrack->Fill(1);
    
  //if (fTrack->GetLabel()<0)return kFALSE;   
  //else fhTotTrack->Fill(2);
   
    //Kinematic Selection
    Float_t Pt = fTrack->Pt();
    if(Pt < fTrkPtMin || Pt > fTrkPtMax)return kFALSE;
    else fhTotTrack->Fill(3);
    
    Float_t Eta = fTrack->Eta();
    if(Eta < fTrkEtaMin || Eta > fTrkEtaMax)return kFALSE;
    else fhTotTrack->Fill(4);
    
    TBits ClusterTPC = fTrack->GetTPCClusterMap();
    if(ClusterTPC.CountBits() < fTrkTPCclus)return kFALSE;
    else fhTotTrack->Fill(5);

    Float_t chi2ndf = fTrack->Chi2perNDF();
    if(chi2ndf >= 4.0)return kFALSE;
    
    Float_t DCAxy = fTrack->DCA();
    Float_t DCAz  = fTrack->ZAtDCA();
    if(fabs(DCAxy) > fTrkabsDCAxy || fabs(DCAz) > fTrkabsDCAz) return kFALSE;
    else fhTotTrack->Fill(6);
    
    return kTRUE;
}

//--------------------------> Selection of Tracks
Bool_t AliAnalysisTaskEbyEChargeSep::IsGeneratedPartPassed(AliAODMCParticle* TorP){
    
    //Kinematic Selection
    Float_t Pt = TorP->Pt();
    if(Pt < fTrkPtMin || Pt > fTrkPtMax)return kFALSE;
    else fhTotPart->Fill(1);
    
    Float_t Eta = TorP->Eta();
    if(Eta < fTrkEtaMin || Eta > fTrkEtaMax)return kFALSE;
    else fhTotPart->Fill(2);
    
    return kTRUE;
}

//--------------------------> Terminate Class
void AliAnalysisTaskEbyEChargeSep::Terminate(Option_t *)
{
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available \n");
        return;
    }
    fDBEvtTree = dynamic_cast<TTree*> (GetOutputData(2));
}


//End

