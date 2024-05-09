// This code includes variabels calculated online within this code. for exp: 2-, 3- particle calculators calculated in this code
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

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliCentrality.h"
#include "AliEventPoolManager.h"

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

#include "AliAnalysisTaskEbyEChargeSep.h"
using namespace std;
ofstream of;

ClassImp(AliAnalysisTaskEbyEChargeSep)
ClassImp(AliMixedEvtTracks)
ClassImp(AliMEOfflineTreeParam)

//-------------------------->Defualt constructor
AliAnalysisTaskEbyEChargeSep::AliAnalysisTaskEbyEChargeSep():
AliAnalysisTaskSE(),
fAOD(0),
fOutputList(0),
farrayMC(0x0),
fEvtCentrality(0),
fEvtCentOrMult(0),
fCentEstimator("V0M"),
fEventCounter(0),
fhTotEvent(0x0),
fhTotTrack(0x0),
fhTotPart(0x0),
fmcArray(0x0),
IsEvtSctPlotReq(kFALSE),
fRunNumPlotReq(135322),
fFileNumPlotReq(720),
fEvtNumPlotReq(112),
fAODProdNumber(0),
fCAODNumber(0),
fRunNumber(0),
fDBEvtTree(0),
fTreeTr(0),
fBit(128),
IsDBMethodOn(kTRUE),
IsCorrltrMethodOn(kTRUE),
fCorrltrLS(kFALSE),
fCorrltrULS(kFALSE),
fReadMC(kFALSE),
fMCKine(kFALSE),
fEvtq2Cor(0),
fEvtq2CorPos(0),
fEvtq2CorNeg(0),
fEvtq3Cor(0),
fEvtq3CorPos(0),
fEvtq3CorNeg(0),
fEvtq3CorUnNPP(0),
fEvtq3CorUnPNN(0),
fEvtq4Cor(0),
fEvtq4CorPos(0),
fEvtq4CorNeg(0),
fDBWindow(90),
fEvtDBValue(0.),
fEvtRndmValue(0.),
DBPMValueMethod(0),
DBPMValueRndm(0),
iPhiMax(0),
fMaxPhi(0),
fASymtry(0),
ASymtryPar(0),
fCompareFileNo(0),
fTreeEventCounter(0),
fEvtabsZvtx(10.),
fEventZvtx(0),
fTrkPtMin(0.2),
fTrkPtMax(5.0),
fTrkEtaMin(-0.8),
fTrkEtaMax(0.8),
fTrkabsDCAxy(2.4),
fTrkabsDCAz(3.2),
fCorrCos3LSPos(0.),
fCorrCos3LSNeg(0.),
fCorrCos3ULSPos(0.),
fCorrCos3ULSNeg(0.),
fTrkTPCclus(70),
fhEffValues(0x0),
fIsEffCorr(kFALSE),
fPoolMngr(0x0),
fEvtPool(0x0),
fMixedEvent(kFALSE),
fMaxTrkDepth(5000),
fMaxEvtMixed(5000),
fMEnTrkPerEvt(10),
fFillMETree(kFALSE),
ftPt_Tr(0.),
ftPhi_Tr(0.),
ftEta_Tr(0.),
ftCharge_Tr(0.),
ftMltOrCnt_Tr(0.),
ftZvtx_Tr(0.),
ftPeriod_Tr(0),
ftOrbit_Tr(0),
ftBC_Tr(0)
{
    //Default Constructor
}

//--------------------------> Named constructor
AliAnalysisTaskEbyEChargeSep::AliAnalysisTaskEbyEChargeSep(const char *name):
AliAnalysisTaskSE(name),
fAOD(0),
fOutputList(0),
farrayMC(0x0),
fEvtCentrality(0),
fEvtCentOrMult(0),
fCentEstimator("V0M"),
fEventCounter(0),
fhTotEvent(0x0),
fhTotTrack(0x0),
fhTotPart(0x0),
fmcArray(0x0),
IsEvtSctPlotReq(kFALSE),
fRunNumPlotReq(135322),
fFileNumPlotReq(720),
fEvtNumPlotReq(112),
fAODProdNumber(0),
fCAODNumber(0),
fRunNumber(0),
fDBEvtTree(0),
fTreeTr(0),
fBit(128),
IsDBMethodOn(kTRUE),
IsCorrltrMethodOn(kTRUE),
fCorrltrLS(kFALSE),
fCorrltrULS(kFALSE),
fReadMC(kFALSE),
fMCKine(kFALSE),
fEvtq2Cor(0),
fEvtq2CorPos(0),
fEvtq2CorNeg(0),
fEvtq3Cor(0),
fEvtq3CorPos(0),
fEvtq3CorNeg(0),
fEvtq3CorUnNPP(0),
fEvtq3CorUnPNN(0),
fEvtq4Cor(0),
fEvtq4CorPos(0),
fEvtq4CorNeg(0),
fDBWindow(90),
fEvtDBValue(0.),
fEvtRndmValue(0.),
DBPMValueMethod(0),
DBPMValueRndm(0),
iPhiMax(0),
fMaxPhi(0),
fASymtry(0),
ASymtryPar(0),
fCompareFileNo(0),
fTreeEventCounter(0),
fEvtabsZvtx(10.),
fEventZvtx(0),
fTrkPtMin(0.2),
fTrkPtMax(5.0),
fTrkEtaMin(-0.8),
fTrkEtaMax(0.8),
fTrkabsDCAxy(2.4),
fTrkabsDCAz(3.2),
fCorrCos3LSPos(0.),
fCorrCos3LSNeg(0.),
fCorrCos3ULSPos(0.),
fCorrCos3ULSNeg(0.),
fTrkTPCclus(70),
fhEffValues(0x0),
fIsEffCorr(kFALSE),
fPoolMngr(0x0),
fEvtPool(0x0),
fMixedEvent(kFALSE),
fMaxTrkDepth(5000),
fMaxEvtMixed(5000),
fMEnTrkPerEvt(10),
fFillMETree(kFALSE),
ftPt_Tr(0.),
ftPhi_Tr(0.),
ftEta_Tr(0.),
ftCharge_Tr(0.),
ftMltOrCnt_Tr(0.),
ftZvtx_Tr(0.),
ftPeriod_Tr(0),
ftOrbit_Tr(0),
ftBC_Tr(0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
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

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    of.open("testoutputfile", ios::out);
    
    TString PartOrTrack = "";
    if(fReadMC){
        PartOrTrack = "McReco";
        if(fMCKine)PartOrTrack = "McPart";
    }else PartOrTrack = "TrkReco";
    
    TString AnaType = "";
    if(fMixedEvent)AnaType = "ME";
    else AnaType = "SE";
    
    Float_t Pi = TMath::Pi();
    fDBEvtTree = new TTree(Form("fDBEvtTree_%s",PartOrTrack.Data()),"Store DB Event Infos");
    fDBEvtTree->Branch("RunNo",&fRunNumber,"fRunNumber/I");
    fDBEvtTree->Branch("AODNumber",&fCAODNumber,"fCAODNumber/I");
    fDBEvtTree->Branch("nEvent",&fTreeEventCounter,"fTreeEventCounter/I");
    fDBEvtTree->Branch("EvtCentrality",&fEvtCentrality,"fEvtCentrality/F");
    fDBEvtTree->Branch("Q2PosPart",   &fEvtq2CorPos, "fEvtq2CorPos/F");
    fDBEvtTree->Branch("Q2NegPart",   &fEvtq2CorNeg, "fEvtq2CorNeg/F");
    fDBEvtTree->Branch("Q3Part",      &fEvtq3Cor,    "fEvtq3Cor/F");
    fDBEvtTree->Branch("Q3PosPart",   &fEvtq3CorPos, "fEvtq3CorPos/F");
    fDBEvtTree->Branch("Q3NegPart",   &fEvtq3CorNeg, "fEvtq3CorNeg/F");
    fDBEvtTree->Branch("Q3PartUnPNN", &fEvtq3CorUnPNN,"fEvtq3CorUnPNN/F");
    fDBEvtTree->Branch("Q3PartUnNPP", &fEvtq3CorUnNPP,"fEvtq3CorUnNPP/F");
    fDBEvtTree->Branch("Q4Part",      &fEvtq4Cor,    "fEvtq4Cor/F");
    fDBEvtTree->Branch("Q4PartPos",   &fEvtq4CorPos,  "fEvtq4CorPos/F");
    fDBEvtTree->Branch("Q4PartNeg",   &fEvtq4CorNeg,  "fEvtq4CorNeg/F");
    fDBEvtTree->Branch("MaxPhi",   &fMaxPhi,        "fMaxPhi/I");
    fDBEvtTree->Branch("Asy",      &fASymtry,       "fASymtry/F");
    fDBEvtTree->Branch(Form("DBValue%d", fDBWindow),&fEvtDBValue,"fEvtDBValue/F");
    fDBEvtTree->Branch(Form("RndmValue%d", fDBWindow),&fEvtRndmValue,"fEvtRndmValue/F");
    if(IsCorrltrMethodOn){
        if(fCorrltrLS){
            fDBEvtTree->Branch("CorrCosLSPos",&fCorrCos3LSPos,"fCorrCos3LSPos/F");
            fDBEvtTree->Branch("CorrCosLSNeg",&fCorrCos3LSNeg,"fCorrCos3LSNeg/F");
        }
        if(fCorrltrULS){
            fDBEvtTree->Branch("CorrCosULSPos",&fCorrCos3ULSPos,"fCorrCos3ULSPos/F");
            fDBEvtTree->Branch("CorrCosULSNeg",&fCorrCos3ULSNeg,"fCorrCos3ULSNeg/F");
        }
    }
    
    
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
    
    
    //Trk Level
    TH3F* fhTorPProp = new TH3F(Form("fhTorPProp%s_%s", PartOrTrack.Data(), AnaType.Data()),  "Properties; Trk p_{T}; Trk #eta; Trk #varphi;", 100, 0, 6.0, 100, -1., 1., 360, 0, 2*Pi);
    fhTorPProp->Sumw2();
    fOutputList->Add(fhTorPProp);
    
    
    if(!fMCKine && !fMixedEvent){
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
    
    //Ztvx Distribution
    TH1F* fhEvtZvtx  = new TH1F("fhEvtZvtx","fhEvtZvtx",240, -12., 12.);
    fhEvtZvtx->Sumw2();
    fOutputList->Add(fhEvtZvtx);
    
    //Multiplicity Distribution
    TH1F* fhMultV0ADet  = new TH1F("fhMultV0ADet","fhMultV0ADet",Int_t(5e2),0,Int_t(10e3));
    TH1F* fhMultV0CDet  = new TH1F("fhMultV0CDet","fhMultV0CDet",Int_t(5e2),0,Int_t(10e3));
    TH1F* fhMultV0ACDet = new TH1F("fhMultV0ACDet","fhMultV0ACDet",Int_t(5e2),0,Int_t(10e3));
    TH1F* fhEvtCentrality = new TH1F("fhEvtCentrality","fhCentrality",12,-2,10);
    TH1F* fhEvtCentralityPer = new TH1F("fhEvtCentralityPer","fhEvtCentralityPer",110,-2,108);
    
    TH2F* fhMultVsCentAllChrg = new TH2F(Form("fhMultVsCentAllChrg_%s", AnaType.Data()),"fhMultVsCentAllChrg",Int_t(2e2),0,Int_t(2e3), 14, -2, 12);
    TH2F* fhMultVsCentPosChrg = new TH2F(Form("fhMultVsCentPosChrg_%s", AnaType.Data()),"fhMultVsCentPosChrg",Int_t(2e2),0,Int_t(2e3), 14, -2, 12);
    TH2F* fhMultVsCentNegChrg = new TH2F(Form("fhMultVsCentNegChrg_%s", AnaType.Data()),"fhMultVsCentNegChrg",Int_t(2e2),0,Int_t(2e3), 14, -2, 12);
    TH3F* fhPtVsChrgVsCent = new TH3F(Form("fhPtVsChrgVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()), Form("fhPtVsChrgVsCent%s",PartOrTrack.Data()), 600, 0, 6, 10, -5, 5, 12, -1, 11);
    TH3F* fhEtaVsChrgVsCent = new TH3F(Form("fhEtaVsChrgVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()), Form("fhEtaVsChrgVsCent%s",PartOrTrack.Data()), 100, -1, 1, 10, -5, 5, 12, -1, 11);
    TH3F* fhPhiVsChrgVsCent = new TH3F(Form("fhPhiVsChrgVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()), Form("fhPhiVsChrgVsCent%s",PartOrTrack.Data()), 360, 0, 2*Pi, 10, -5, 5, 12, -1, 11);
    
    fhMultV0ADet->Sumw2();
    fhMultV0CDet->Sumw2();
    fhMultV0ACDet->Sumw2();
    fhMultVsCentAllChrg->Sumw2();
    fhMultVsCentPosChrg->Sumw2();
    fhMultVsCentNegChrg->Sumw2();
    fhEvtCentrality->Sumw2();
    fhEvtCentralityPer->Sumw2();
    fhPtVsChrgVsCent->Sumw2();
    fhEtaVsChrgVsCent->Sumw2();
    fhPhiVsChrgVsCent->Sumw2();
    
    fOutputList->Add(fhMultV0ADet);
    fOutputList->Add(fhMultV0CDet);
    fOutputList->Add(fhMultV0ACDet);
    fOutputList->Add(fhEvtCentrality);
    fOutputList->Add(fhEvtCentralityPer);
    fOutputList->Add(fhMultVsCentAllChrg);
    fOutputList->Add(fhMultVsCentPosChrg);
    fOutputList->Add(fhMultVsCentNegChrg);
    fOutputList->Add(fhPtVsChrgVsCent);
    fOutputList->Add(fhEtaVsChrgVsCent);
    fOutputList->Add(fhPhiVsChrgVsCent);
    
    
    //1. DB Distribution
    //Method Wise
    if(IsDBMethodOn){
        TH2F* fhDBwWindowVsCent  = new TH2F(Form("fhDBwW%dVsCent%s_%s",fDBWindow,PartOrTrack.Data(), AnaType.Data()),Form("fhDBwW%dVsCent%s",fDBWindow,PartOrTrack.Data()), 100, 0, 5, 12, -1, 11);
        fhDBwWindowVsCent->Sumw2();
        fOutputList->Add(fhDBwWindowVsCent);
        
        
        TH2F* fhRndmwWindowVsCent  = new TH2F(Form("fhRndwW%dVsCent%s_%s",fDBWindow, PartOrTrack.Data(), AnaType.Data()),Form("fhRndwW%dVsCent%s",fDBWindow,PartOrTrack.Data()), 100, 0, 5, 12, -1, 11);
        fhRndmwWindowVsCent->Sumw2();
        fOutputList->Add(fhRndmwWindowVsCent);
    }
    
    //2. Correlator method
    if(IsCorrltrMethodOn){
      /*
        //V0 Event Plane Distribution
        TH1F* fhV0Plane   = new TH1F("fhV0Plane","fhV0Plane", 240, -20, 200);
        TH1F* fhV0APlane  = new TH1F("fhV0APlane","fhV0APlane", 240, -20, 200);
        TH1F* fhV0CPlane  = new TH1F("fhV0CPlane","fhV0CPlane", 240, -20, 200);
        fhV0Plane->Sumw2();
        fhV0APlane->Sumw2();
        fhV0CPlane->Sumw2();
        fOutputList->Add(fhV0Plane);
        fOutputList->Add(fhV0APlane);
        fOutputList->Add(fhV0CPlane);
        
        
        if(fCorrltrLS){ //like sign
            TH2F* fhCorrCos3SumlsPP  = new TH2F(Form("fhCorrCos3pLSVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()),Form("fhCorrCos3pLSVsCent%s",PartOrTrack.Data()), 12, -1, 11, 100, -0.001, 0.001);
            fhCorrCos3SumlsPP->Sumw2();
            fOutputList->Add(fhCorrCos3SumlsPP);
            
            TH2F* fhCorrCos3SumlsNN  = new TH2F(Form("fhCorrCos3nLSVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()),Form("fhCorrCos3nLSVsCent%s",PartOrTrack.Data()), 12, -1, 11, 100, -0.001, 0.001);
            fhCorrCos3SumlsNN->Sumw2();
            fOutputList->Add(fhCorrCos3SumlsNN);
        }
        
        if(fCorrltrULS){//unlike sign
            TH2F* fhCorrCos3SumUlsPP  = new TH2F(Form("fhCorrCos3pULSVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()),Form("fhCorrCos3pULSVsCent%s",PartOrTrack.Data()), 12, -1, 11, 100, -0.001, 0.001);
            fhCorrCos3SumUlsPP->Sumw2();
            fOutputList->Add(fhCorrCos3SumUlsPP);
            
            TH2F* fhCorrCos3SumUlsNN  = new TH2F(Form("fhCorrCos3nULSVsCent%s_%s",PartOrTrack.Data(), AnaType.Data()),Form("fhCorrCos3nULSVsCent%s",PartOrTrack.Data()), 12, -1, 11, 100, -0.001, 0.001);
            fhCorrCos3SumUlsNN->Sumw2();
            fOutputList->Add(fhCorrCos3SumUlsNN);
        }*/
    }
    
    //Step2. Analysis Plot
    if(IsEvtSctPlotReq){
        TH2F* fHXYScatterVsCentPlotPos = new TH2F(Form("fhXYSctPlotPos%s_%s",PartOrTrack.Data(), AnaType.Data()), "fHXYScatterPlotPos", 600, -300, 300, 600, -300, 300);
        fHXYScatterVsCentPlotPos->Sumw2();
        fOutputList->Add(fHXYScatterVsCentPlotPos);
        
        TH2F* fHXYScatterVsCentPlotNeg = new TH2F(Form("fhXYSctPlotNeg%s_%s",PartOrTrack.Data(), AnaType.Data()), "fHXYScatterPlotNeg", 600, -300, 300, 600, -300, 300);
        fHXYScatterVsCentPlotNeg->Sumw2();
        fOutputList->Add(fHXYScatterVsCentPlotNeg);
    }
    
    //Event Mixing attempt
    
    //Defination
    Int_t nBinsCentorMult = 2;
    Double_t centralityBins[] = {20, 40, 70};
    Int_t nZvtxBins = 4;
    Double_t vertexBins[] = {-10, -2.5, 0., 2.5, 10.};
    
    fPoolMngr = new AliEventPoolManager(fMaxEvtMixed,fMaxTrkDepth,nBinsCentorMult,(Double_t*)centralityBins,nZvtxBins,(Double_t*) vertexBins);
    
    //Offline TREE for the first time
    if(fFillMETree && fMixedEvent) {
        fTreeTr = new TTree(Form("fTreeTr_%s", AnaType.Data()),"TTree for Event Tracks");
        fTreeTr->Branch("Pt_Tr",      &ftPt_Tr,      "ftPt_Tr/F");
        fTreeTr->Branch("Phi_Tr",     &ftPhi_Tr,     "ftPhi_Tr/F");
        fTreeTr->Branch("Eta_Tr",     &ftEta_Tr,     "ftEta_Tr/F");
        fTreeTr->Branch("Charge_Tr",  &ftCharge_Tr,  "ftCharge_Tr/F");
        fTreeTr->Branch("MltOrCnt_Tr",&ftMltOrCnt_Tr,"ftMltOrCnt_Tr/F");
        fTreeTr->Branch("Zvtx_Tr",    &ftZvtx_Tr,    "ftZvtx_Tr/F");
        fTreeTr->Branch("Period_Tr",  &ftPeriod_Tr,  "ftPeriod_Tr/I");
        fTreeTr->Branch("Orbit_Tr",   &ftOrbit_Tr,   "ftOrbit_Tr/I");
        fTreeTr->Branch("BC_Tr",      &ftBC_Tr,      "ftBC_Tr/I");
        PostData(3, fTreeTr);
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
    
    fCAODNumber = GetFileNumberOfCurrentAOD();
    
    if(fCompareFileNo != fCAODNumber){
        fCompareFileNo = fCAODNumber;
        fTreeEventCounter = 0;
    }
    
    fTreeEventCounter++;
    
    if(IsEvtSctPlotReq){
        if(fRunNumPlotReq != fAOD->GetRunNumber())return;
        if(fFileNumPlotReq != GetFileNumberOfCurrentAOD())return;
        if(fEvtNumPlotReq != fTreeEventCounter)return;
        // Going Back if it's not mached here
    }
    
    TObjArray* fTrackArray = new TObjArray;
    fTrackArray->SetOwner(kTRUE);
    
    //Storing Track or KinePart in Array !
    Int_t nTrackOrPart = 0;
    if(fMCKine)nTrackOrPart = fmcArray->GetEntriesFast();
    else nTrackOrPart = fAOD->GetNumberOfTracks();
    
    for (Int_t iTrack = 0; iTrack < nTrackOrPart; iTrack++){
        
        if(fMCKine){ // Generated Particles
            AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fmcArray->At(iTrack));
            if (!mcPart) continue;
            fhTotPart->Fill(0);
            if(!IsGeneratedPartPassed(mcPart))continue;
            if(!mcPart->IsPhysicalPrimary())continue;
            fhTotPart->Fill(3);
            fTrackArray->Add(mcPart);
            
        } else{ // Reco Particles
            AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrack);
            if (!track)continue;
            
            fhTotTrack->Fill(0);
            if(!IsTrackSelectionPassed(track)) continue;
            if (track->Charge()==0)continue;
            if(fReadMC){
                AliAODMCParticle *TorP  = (AliAODMCParticle*)fMCEvent->GetTrack(track->GetLabel());
                if (!TorP->IsPhysicalPrimary()) continue;
            }
            fTrackArray->Add(track);
        }
    }
    
    if(fFillMETree && fMixedEvent)OfflineMETree(fAOD, fTrackArray);
    
    //for Mixed Event Pool
    TObjArray* fTrackMEArrayTemp = new TObjArray;
    fTrackMEArrayTemp->SetOwner(kTRUE);
    
    if(fMixedEvent){
        if(!fPoolMngr)return;
        fEvtPool = fPoolMngr->GetEventPool(fEvtCentOrMult, fEventZvtx);
        if (!fEvtPool){
            //AliInfo(Form("No-1: pool found for multiplicity = %f, zVtx = %f cm", fEvtCentOrMult, fEventZvtx));
            return;
        }
        
        if(!fEvtPool->IsReady()){
            //AliInfo(Form("No-2: pool not ready for multiplicity = %f, zVtx = %f cm", fEvtCentOrMult, fEventZvtx));
        }else cout <<"Pool is ready" << endl;
        
        fTrackMEArrayTemp = (TObjArray*)StoreMEEventsandTracks(fTrackArray);
    }
    
    
    //Filtered TrackLoop:
    Int_t lMultiAllChrg, lMultPosChrg, lMultNegChrg;
    Int_t NTrkPosEta=0, NTrkNegEta=0, TorPPhiRadiBin=0;
    Float_t TorPPhiRad = 0.;
    Float_t Cos2Phi=0.,Sin2Phi=0.;
    Float_t PtWeightedCos2Phi=0.,PtWeightedSin2Phi=0.;
    Float_t TrkAryAllChrgPhi[Int_t(5e3)], TrkAryPosChrgPhi[Int_t(5e3)],TrkAryNegChrgPhi[Int_t(5e3)];
    Float_t TrkAryPosEtaPhi[Int_t(5e3)],TrkAryNegEtaPhi[Int_t(5e3)];
    Float_t NTrkAryPosChrgPhiRadBin[361],NTrkAryNegChrgPhiRadBin[361];
    Float_t qx,qy,q2x,q2y,qxp,qyp,q2xp,q2yp,q4x,q4y,q4xp,q4yp;
    Float_t qxn,qyn,q2xn,q2yn,q4yn,q4xn;
    
    for(Int_t iPhiDegree=0; iPhiDegree<=360; iPhiDegree++){
        NTrkAryPosChrgPhiRadBin[iPhiDegree]=0.;
        NTrkAryNegChrgPhiRadBin[iPhiDegree]=0.;
    }

    qx =0.,qy=0.,q2x=0.0,q2y=0.0, q4x=0.0,q4y=0.0;
    qxp =0.,qyp=0.,q2xp=0.0,q2yp=0.0,q4xn=0.0,q4yn=0.0;
    qxn =0.,qyn=0.,q2xn=0.0,q2yn=0.0,q4xp=0.0,q4yp=0.0;
    lMultiAllChrg = 0, lMultPosChrg = 0, lMultNegChrg = 0;
    
    Int_t NEvents = 0;// SE or ME
    if(fMixedEvent){
        if(!fEvtPool->IsReady()) NEvents =0;
        else {
            NEvents = fEvtPool->GetCurrentNEvents();
            Int_t diffEvt = NEvents*fMEnTrkPerEvt - fTrackArray->GetEntriesFast();
            if(TMath::Abs(diffEvt)>50)NEvents = 0; // 50 = mult difference
            
            TObjArray*  fSEMEEvtTrackstemp = new TObjArray;
            fSEMEEvtTrackstemp->SetOwner(kTRUE);
            Bool_t IsEventPoolReady =  kFALSE;
            
            Int_t fMEMultTemp = 0.;
            for (Int_t jMix =0; jMix < NEvents; jMix++){
                
                fSEMEEvtTrackstemp = fEvtPool->GetEvent(jMix); // replacing from pool
                if(!fSEMEEvtTrackstemp)continue;
                
                for(int i=0; i < fSEMEEvtTrackstemp->GetEntriesFast(); i++)fMEMultTemp++;
            }
            
            if(fMEMultTemp >= fTrackArray->GetEntriesFast())IsEventPoolReady = kTRUE;
            if(!IsEventPoolReady)return;
        }
    }else NEvents = 1;// SE
    
    
    
    TObjArray*  fSEMEEvtTracks = new TObjArray;
    fSEMEEvtTracks->SetOwner(kTRUE);
    
    Int_t fMEMultTemp2 = 0.;
    
    //cout << "Total Events Number = " << NEvents << endl;
    for (Int_t jMix =0; jMix < NEvents; jMix++){
        
        //cout << " -- Event Number = " << jMix << endl;
        if(!fMixedEvent)fSEMEEvtTracks = (TObjArray*)fTrackArray; //Stored Track or Part
        else if(fMixedEvent)fSEMEEvtTracks = fEvtPool->GetEvent(jMix); // replacing from pool
        
        if(!fSEMEEvtTracks)continue;
        
        //cout << " ---- Entries = " << fSEMEEvtTracks->GetEntriesFast() << endl;
        for(int iTorPStr=0; iTorPStr < fSEMEEvtTracks->GetEntriesFast(); iTorPStr++){
            
            //constaint of MEmult=SEmult
            fMEMultTemp2++;
            if(fMEMultTemp2 > fTrackArray->GetEntriesFast())break;
            //cout << fMEMultTemp2 << ") " << fTrackArray->GetEntriesFast() << endl;
            
            Float_t TorPP = 0, TorPPt =0., TorPEta = 0., TorPPhi=0, TorPPhiDeg = 0;
            Short_t TorPChrg =0;
	       
            if(fMCKine){
                AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fSEMEEvtTracks->At(iTorPStr));
                if (!mcPart){
                    AliError(Form("Failed casting particle from MC array!, Skipping particle for %d", iTorPStr));
                    continue;
                }
                TorPP = mcPart->P();
                TorPPt = mcPart->Pt();
                TorPPhi = mcPart->Phi();
                TorPEta = mcPart->Eta();
                TorPChrg = mcPart->Charge();
                TorPPhiDeg = TorPPhi*180.0/TMath::Pi();
                
            }else {
                AliAODTrack* track = (AliAODTrack*)fSEMEEvtTracks->At(iTorPStr);
                if (!track){
                    AliError(Form("Failed casting particle from track array!, Skipping particle for %d", iTorPStr));
                    continue;
                }
                
                //cout << "track-P() >>" << track->P() <<",   Charge >> " <<  track->Charge() << endl;
                TorPP = track->P();
                TorPPt = track->Pt();
                TorPPhi = track->Phi();
                TorPEta = track->Eta();
                TorPChrg = track->Charge();
                TorPPhiDeg = TorPPhi*180.0/TMath::Pi();
	            
                if(!fMixedEvent){
                    ((TH1F*)fOutputList->FindObject(Form("fhTrkDCAxy_%s", AnaType.Data())))->Fill(track->DCA());
                    ((TH1F*)fOutputList->FindObject(Form("fhTrkDCAz_%s", AnaType.Data())))->Fill(track->ZAtDCA());
                    ((TH1F*)fOutputList->FindObject(Form("fhTrkClusterTPC_%s", AnaType.Data())))->Fill(track->GetTPCNcls());
                    ((TH1F*)fOutputList->FindObject(Form("fhTrkChi2NDF_%s", AnaType.Data())))->Fill(track->Chi2perNDF());
                    ((TH2F*)fOutputList->FindObject(Form("fhTrkChrgVsDCAxy_%s", AnaType.Data())))->Fill(track->DCA(), TorPChrg+0.5);
                    ((TH2F*)fOutputList->FindObject(Form("fhTrkChrgVsDCAz_%s", AnaType.Data())))->Fill(track->ZAtDCA(), TorPChrg+0.5);
                    fhTotTrack->Fill(7);
                }
            }
            
            //cout << iTorPStr << ") TEST Pt = " << TorPPt << endl;
	    if(TorPChrg==0)continue;
            ((TH3F*)fOutputList->FindObject(Form("fhTorPProp%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPt,TorPEta,TorPPhi);
            ((TH2F*)fOutputList->FindObject(Form("fhChrgVsPt%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPt, TorPChrg+0.5);
            ((TH2F*)fOutputList->FindObject(Form("fhChrgVsPhi%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPhi, TorPChrg+0.5);
            ((TH2F*)fOutputList->FindObject(Form("fhChrgVsPhiD%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPhiDeg, TorPChrg+0.5);
	    ((TH2F*)fOutputList->FindObject(Form("fhChrgVsEta%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPEta, TorPChrg+0.5);
            ((TH3F*)fOutputList->FindObject(Form("fhPtVsChrgVsCent%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPt, TorPChrg+0.5, fEvtCentrality+0.5);
            ((TH3F*)fOutputList->FindObject(Form("fhEtaVsChrgVsCent%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPEta, TorPChrg+0.5, fEvtCentrality+0.5);
            ((TH3F*)fOutputList->FindObject(Form("fhPhiVsChrgVsCent%s_%s",PartOrTrack.Data(),AnaType.Data())))->Fill(TorPPhi, TorPChrg+0.5, fEvtCentrality+0.5);
            
            lMultiAllChrg++;
            TrkAryAllChrgPhi[lMultiAllChrg]=TorPPhi;
	    qx +=cos(TorPPhi);
	    q2x +=cos(2.*TorPPhi);
	    q4x +=cos(4.*TorPPhi);
	    qy +=sin(TorPPhi);
	    q2y +=sin(2.*TorPPhi);
	    q4y +=sin(4.*TorPPhi);

            if(TorPChrg > 0){
                lMultPosChrg++;
                TrkAryPosChrgPhi[lMultPosChrg]=0.0;
                TrkAryPosChrgPhi[lMultPosChrg]=TorPPhi;
		qxp  +=cos(TorPPhi);
		q2xp +=cos(2.*TorPPhi);
		q4xp +=cos(4.*TorPPhi);
		qyp  +=sin(TorPPhi);
		q2yp +=sin(2.*TorPPhi);
		q4yp +=sin(4.*TorPPhi);

            }else if(TorPChrg < 0){
                lMultNegChrg++;
                TrkAryNegChrgPhi[lMultNegChrg]=0.0;
                TrkAryNegChrgPhi[lMultNegChrg]=TorPPhi;
		qxn +=cos(TorPPhi);
		q2xn +=cos(2.*TorPPhi);
		q4xn +=cos(4.*TorPPhi);
		qyn +=sin(TorPPhi);
		q2yn +=sin(2.*TorPPhi);
		q4yn +=sin(4.*TorPPhi);
            }
            
            //Eta based QA and calculations
            if(TorPEta > 0){
                NTrkPosEta++;
                TrkAryPosEtaPhi[NTrkPosEta]=0.0;
                TrkAryPosEtaPhi[NTrkPosEta]=TorPPhi;
                
            }else if(TorPEta <= 0){
                NTrkNegEta++;
                TrkAryNegEtaPhi[NTrkNegEta]=0.0;
                TrkAryNegEtaPhi[NTrkNegEta]=TorPPhi;
            }
            
            if(IsEvtSctPlotReq && TorPP!=0){
                Float_t ScatXValue = 200.0*TorPPt*cos(TorPPhi)/TorPP;
                Float_t ScatYValue =200.0*TorPPt*sin(TorPPhi)/TorPP;
                if(TorPChrg > 0)((TH2F*)fOutputList->FindObject(Form("fhXYSctPlotPos%s_%s",PartOrTrack.Data(), AnaType.Data())))->Fill(ScatXValue,ScatYValue);
                else if (TorPChrg < 0)((TH2F*)fOutputList->FindObject(Form("fhXYSctPlotNeg%s_%s",PartOrTrack.Data(), AnaType.Data())))->Fill(ScatXValue,ScatYValue);
            }
            
            Cos2Phi += cos(2.0*TorPPhi);
            Sin2Phi += sin(2.0*TorPPhi);
            
            Float_t PtWeight = TorPPt;
            if(PtWeight > 2.)PtWeight= 2;
            PtWeightedCos2Phi += PtWeight*cos(2.0*TorPPhi);
            PtWeightedSin2Phi += PtWeight*sin(2.0*TorPPhi);
            
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
	  
	}
    }
    
    if(fMixedEvent){
        TObjArray* tracksMEClone = CloneAcceptAndReduceTracks(fTrackMEArrayTemp);
        fEvtPool->UpdatePool(tracksMEClone);
        if(!fEvtPool->IsReady()) return;
    }

    cout<<lMultiAllChrg<<"    "<<lMultPosChrg<<"     "<<lMultNegChrg<<endl;

    ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentAllChrg_%s", AnaType.Data())))->Fill(lMultiAllChrg, fEvtCentrality+0.5);
    ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentPosChrg_%s", AnaType.Data())))->Fill(lMultPosChrg, fEvtCentrality+0.5);
    ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentNegChrg_%s", AnaType.Data())))->Fill(lMultNegChrg, fEvtCentrality+0.5);
    
    // Going back if second step of RUN is requested
    if(IsEvtSctPlotReq)return;
     
    Float_t q2Cor    = Cal2PartCor(q2x,q2y,lMultiAllChrg);  // same charge                                                                 
    Float_t q2CorPos = Cal2PartCor(q2xp,q2yp,lMultPosChrg); //                                                                             
    Float_t q2CorNeg = Cal2PartCor(q2xn,q2yn,lMultNegChrg); //                                                                             
    fEvtq2Cor    = q2Cor;
    fEvtq2CorPos = q2CorPos;
    fEvtq2CorNeg = q2CorNeg;

    //------------->3-particle (same charge)                                                                                               
    Float_t q3Cor    = Cal3PartCor(qx,qy,q2x,q2y,lMultiAllChrg);
    Float_t q3CorPos = Cal3PartCor(qxp,qyp,q2xp,q2yp,lMultPosChrg); // same ch                                                             
    Float_t q3CorNeg = Cal3PartCor(qxn,qyn,q2xn,q2yn,lMultNegChrg); // same ch                                                             
    fEvtq3Cor    = q3Cor;
    fEvtq3CorPos = q3CorPos;
    fEvtq3CorNeg = q3CorNeg;

    //................>3-particle (unlike charge)                                                                                          
    Float_t q3CorPNN = Cal3PartCorUn(qxp,qyp,q2xp,q2yp,qxn,qyn,q2xn,q2yn,lMultPosChrg,lMultNegChrg); //+--                                 
    Float_t q3CorNPP = Cal3PartCorUn(qxn,qyn,q2xn,q2yn,qxp,qyp,q2xp,q2yp,lMultNegChrg,lMultPosChrg); //-++                                 
    fEvtq3CorUnPNN = q3CorPNN; //q3unp                                                                                                     
    fEvtq3CorUnNPP = q3CorNPP; //q3unn                                                                                                     

    //------------->4-particle                                                                                                             
    Float_t q4Cor    = Cal4PartCor(q2x,q2y,q4x,q4y,lMultiAllChrg);
    Float_t q4CorPos = Cal4PartCor(q2xp,q2yp,q4xp,q4yp,lMultPosChrg);
    Float_t q4CorNeg = Cal4PartCor(q2xn,q2yn,q4xn,q4yn,lMultNegChrg);
    fEvtq4Cor    = q4Cor;
    fEvtq4CorPos = q4CorPos;
    fEvtq4CorNeg = q4CorNeg;

    //1) DB Extraction with phi window: +90Max
    if(IsDBMethodOn){
        Float_t TempTotPosChargArray[451],TempTotNegChargArray[451];
        for(Int_t iPhi =1; iPhi<=360; iPhi++){
            TempTotPosChargArray[iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
            TempTotNegChargArray[iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
            if(iPhi <=90){
                TempTotPosChargArray[360+iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
                TempTotNegChargArray[360+iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
            }
        }
        
        //DBMethod
        //Float_t fOrgDBValue = CalEvtPhiDumbellMethod(fDBWindow,TempTotPosChargArray,TempTotNegChargArray, kTRUE);
	CalEvtPhiDumbellMethod(90,TempTotPosChargArray,TempTotNegChargArray);
	fEvtDBValue = DBPMValueMethod;
	fMaxPhi = iPhiMax;
	fASymtry = ASymtryPar;
        ((TH2F*)fOutputList->FindObject(Form("fhDBwW%dVsCent%s_%s",fDBWindow, PartOrTrack.Data(), AnaType.Data())))->Fill(fEvtDBValue, fEvtCentrality+0.5);

        //RandomMethod
        //Float_t fRndmDBValue = CalEvtPhiDumbellMethod(fDBWindow,TempTotPosChargArray,TempTotNegChargArray, kFALSE);
        ((TH2F*)fOutputList->FindObject(Form("fhRndwW%dVsCent%s_%s",fDBWindow,PartOrTrack.Data(), AnaType.Data())))->Fill(fEvtRndmValue, fEvtCentrality+0.5);
        fEvtRndmValue = DBPMValueRndm;
    }
    
    
    //2) Correlator Method
    if(IsCorrltrMethodOn){
      /*
        //2a) Two Particle Correlator Like Sign with Event Plance
        Float_t Corr2PCandEvtPAllTrk     = CalCorlVia2PCLikeSignwEvtP(lMultiAllChrg,TrkAryAllChrgPhi);
        Float_t Corr2PCandEvtPLikeSign   = CalCorlVia2PCLikeSignwEvtP(NTrkPosEta,TrkAryPosEtaPhi);
        Float_t Corr2PCandEvtPUnlikeSign = CalCorlVia2PCLikeSignwEvtP(NTrkNegEta,TrkAryNegEtaPhi);
        
        //2b) Two Particle Correlator via Like Sign
        Float_t EvtPhiRad = ExtractEvtPhiInProperRange(Sin2Phi, Cos2Phi);
        Float_t EvtPhiDegree = EvtPhiRad*180.0/TMath::Pi();
        
        Float_t EvtPhiWeightRad    = ExtractEvtPhiInProperRange(PtWeightedSin2Phi, PtWeightedCos2Phi);
        Float_t EvtPhiWeightDegree = EvtPhiWeightRad*180.0/TMath::Pi();
        
        AliEventplane* Eventplane =  fAOD->GetEventplane();
        //Float_t V0Plane  = Eventplane->GetEventplane("V0", fAOD, 2)*180.0/TMath::Pi();
        Float_t V0Plane  = Eventplane->GetEventplane("V0", fAOD, 2);
        Float_t V0APlane = Eventplane->GetEventplane("V0A",fAOD, 2);
        Float_t V0CPlane = Eventplane->GetEventplane("V0C",fAOD, 2);
        
        if(V0Plane  < 0)V0Plane  =  180. + V0Plane;
        if(V0APlane < 0)V0APlane =  180. + V0APlane;
        if(V0CPlane < 0)V0CPlane =  180. + V0CPlane;
        
        ((TH1F*)fOutputList->FindObject("fhV0Plane"))->Fill(V0Plane);
        ((TH1F*)fOutputList->FindObject("fhV0APlane"))->Fill(V0APlane);
        ((TH1F*)fOutputList->FindObject("fhV0CPlane"))->Fill(V0CPlane);
        
	if(fCorrltrLS){
            CalCorlVia2PCLikeSign(lMultPosChrg, EvtPhiRad, EvtPhiWeightRad, V0Plane, V0APlane, V0CPlane,TrkAryPosChrgPhi,kTRUE);
            CalCorlVia2PCLikeSign(lMultNegChrg, EvtPhiRad, EvtPhiWeightRad, V0Plane, V0APlane, V0CPlane, TrkAryNegChrgPhi,kFALSE);
        }
        
        //2c) Two Particle Correlator via Unlike Sign
        if(fCorrltrULS){
            CalCorlVia2PCUnlikeSign(lMultPosChrg, lMultNegChrg, EvtPhiRad, EvtPhiWeightRad, V0Plane, V0APlane, V0CPlane, TrkAryPosChrgPhi, TrkAryNegChrgPhi, kTRUE);
            CalCorlVia2PCUnlikeSign(lMultNegChrg, lMultPosChrg, EvtPhiRad, EvtPhiWeightRad, V0Plane, V0APlane, V0CPlane, TrkAryNegChrgPhi, TrkAryPosChrgPhi, kFALSE);
        }*/
    }
    
    fDBEvtTree->Fill();
    PostData(1, fOutputList);
    PostData(2, fDBEvtTree);
    
}

Float_t  AliAnalysisTaskEbyEChargeSep::Cal2PartCor(Float_t cos2x,Float_t sin2y, Float_t nt)   // same charge (++, --)                          
{
  Float_t qc =0.;
  qc =cos2x*cos2x+sin2y*sin2y;
  //  if((nt*(nt-1.0)) > 0) qc = (qc-nt)/(nt*(nt-1.0));                                                                                        
  qc = qc-nt;
  return qc;
}

Float_t  AliAnalysisTaskEbyEChargeSep::Cal3PartCor(float cosx,float siny,float cos2x,float sin2y, Int_t nt)  //same charge                     
{
  Float_t qc =0., q2c =0., q3c =0., q3 =0.;
  qc  = cosx*cosx+siny*siny;
  q2c = cos2x*cos2x+sin2y*sin2y;
  q3c = cos2x*(cosx*cosx-siny*siny)+2.0*sin2y*cosx*siny;
  q3 = q3c-2.0*qc-q2c+2.0*nt;
  return q3;
}

Float_t  AliAnalysisTaskEbyEChargeSep::Cal3PartCorUn(float cosxp,float sinyp,float cos2xp,float sin2yp,float cosxn,float sinyn,float cos2xn,float sin2yn, Int_t fpos, Int_t fneg) // opp charge                                                                                              
{
  float qqq = cosxp*cosxn*cos2xn+cosxp*sinyn*sin2yn+sinyp*cosxn*sin2yn-sinyp*sinyn*cos2xn;
  float qq  = (cos2xp*cos2xn+sin2yp*sin2yn);
  float re  = cosxp*cosxn+sinyp*sinyn;
  float  q3 = qqq-re;
  //float  q3 = qqq-re/(fpos*(fpos-1.0)*fneg);                                                                                                 
  return q3;
}

Float_t  AliAnalysisTaskEbyEChargeSep::Cal4PartCor(float cos2x,float sin2y,float cos4x,float sin4y, Int_t nt)
{
  Float_t q4c;
  float  qc  = cos2x*cos2x+sin2y*sin2y;
  float  q2c = cos4x*cos4x+sin4y*sin4y;
  float q3c  = cos4x*(cos2x*cos2x-sin2y*sin2y)+2.0*sin4y*cos2x*sin2y;
  //   float q3c = cos4x*(cos2x**2.-sin2y**2.);                                                                                                
  float wt1 = (nt*(nt-1.)*(nt-2.)*(nt-3.));
  if(wt1>0) q4c = (qc*qc+q2c-2.*q3c-4.*(nt-2.)*qc)/float(wt1);
  //  if(((nt-1.)*(nt-2.)) >0) {                                                                                                               
  float wt = (nt-1.)*(nt-2.);
  if(wt>0) q4c +=2./float(wt);
  return q4c;
}

//--------------------------> Calculated Eventt Phi via Dumbell Method
//Float_t AliAnalysisTaskEbyEChargeSep::CalEvtPhiDumbellMethod(Int_t fPhiWindow, Float_t *TempTotChrgPos,Float_t *TempTotChrgNeg, Bool_t isMethod){
void AliAnalysisTaskEbyEChargeSep::CalEvtPhiDumbellMethod(Int_t fPhiWindow, Float_t *TempTotChrgPos,Float_t *TempTotChrgNeg)
{  
    Float_t SumPCInWindow[361],SumNCInWindow[361];
    Int_t tempPhi =0;
    
    for(Int_t iPhi=1; iPhi<=360; iPhi++){
        tempPhi = 0;
        SumPCInWindow[iPhi]=0.;
        SumPCInWindow[iPhi] = TempTotChrgPos[iPhi];
        
        SumNCInWindow[iPhi]=0.;
        SumNCInWindow[iPhi] = TempTotChrgNeg[iPhi];
        
        for(Int_t jPhi=iPhi+1; jPhi<=fPhiWindow-1+iPhi; jPhi++){
            tempPhi =jPhi;
            if (jPhi>360)tempPhi=tempPhi-360;
            SumPCInWindow[iPhi] = SumPCInWindow[iPhi]+TempTotChrgPos[tempPhi];
            SumNCInWindow[iPhi] = SumNCInWindow[iPhi]+TempTotChrgNeg[tempPhi];
        }
    }

    Float_t RefFracAB=0.,FracAB=0.,FracAPosFront=0.,FracBNegBack=0, ASymtry=999., maxpn=0., minnp=0.;
    //    Float_t DBPMValueMethod = 0.0, DBPMValueRndm = -9.0;
    Int_t iPhiPosFront180=0, iPhiNegBack180=0;
    
    //METHOD1: Per Dumb-bell DB values
    for(Int_t iPhi=1; iPhi<=360; iPhi++){
        if(SumPCInWindow[iPhi]<0 || SumNCInWindow[iPhi]<0)continue; //no case but protection
        if((SumPCInWindow[iPhi] + SumNCInWindow[iPhi])<=0) continue; //0 possible
        
        iPhiPosFront180 = iPhi;
        if(iPhi<=180)iPhiNegBack180 = iPhi+180;
        else if(iPhi>180)iPhiNegBack180 = iPhi-180;
        
        if((SumPCInWindow[iPhiPosFront180]+SumNCInWindow[iPhiPosFront180]) <= 0)continue;
        FracAPosFront = SumPCInWindow[iPhiPosFront180]/(SumPCInWindow[iPhiPosFront180]+SumNCInWindow[iPhiPosFront180]);
        
        if((SumPCInWindow[iPhiNegBack180]+SumNCInWindow[iPhiNegBack180]) <= 0)continue;
        FracBNegBack  = SumNCInWindow[iPhiNegBack180]/(SumPCInWindow[iPhiNegBack180]+SumNCInWindow[iPhiNegBack180]);

	if(FracAPosFront < 0.5 || FracBNegBack < 0.5) continue;
	if((FracAPosFront+FracBNegBack)<=0)continue; //protection

	maxpn =SumPCInWindow[iPhiPosFront180]-SumNCInWindow[iPhiPosFront180];    
	minnp =SumNCInWindow[iPhiNegBack180]-SumPCInWindow[iPhiNegBack180];

	if((maxpn+minnp) ==0) ASymtry=999;
	if((maxpn+minnp)>0)  ASymtry = (maxpn-minnp)/(maxpn+minnp);
	if(ASymtry <0) ASymtry = -ASymtry;
	if(ASymtry >0.25) continue;

	FracAB = FracAPosFront + FracBNegBack;

        if(FracAB > RefFracAB){ //leading sum
            RefFracAB=FracAB;
            DBPMValueMethod=FracAB;
            iPhiMax=iPhi;
	    ASymtryPar=ASymtry;
	}
    }

    //METHOD2: Random DB Distribution (Validity KindOf)
    Float_t FracABRndm=0.,FracAPosFrontRndm=0.,FracBNegBackRndm=0;
    Int_t iPhiRndmPosFront=0, iPhiRndmNegBack=0;
    Int_t iPhiRndm = gRandom->Uniform(1,360);
    
    iPhiRndmPosFront = iPhiRndm;
    if(iPhiRndm<=180)iPhiRndmNegBack = iPhiRndm+180;
    else if(iPhiRndm>180)iPhiRndmNegBack = iPhiRndm-180;
    
    if((SumPCInWindow[iPhiRndmPosFront]+SumNCInWindow[iPhiRndmPosFront]) > 0){
        FracAPosFrontRndm = SumPCInWindow[iPhiRndmPosFront]/(SumPCInWindow[iPhiRndmPosFront]+SumNCInWindow[iPhiRndmPosFront]);
    }
    
    if((SumPCInWindow[iPhiRndmNegBack]+SumNCInWindow[iPhiRndmNegBack]) > 0 ){
        FracBNegBackRndm  = SumNCInWindow[iPhiRndmNegBack]/(SumPCInWindow[iPhiRndmNegBack]+SumNCInWindow[iPhiRndmNegBack]);
    }
    
    FracABRndm = FracAPosFrontRndm + FracBNegBackRndm;
    if(FracABRndm>0)DBPMValueRndm = FracABRndm;
    
    /*    if(isMethod)return DBPMValueMethod;
	  else return DBPMValueRndm;*/
}

/*
//--------------------------> Correlator 2 Particle
Float_t AliAnalysisTaskEbyEChargeSep::CalCorlVia2PCLikeSignwEvtP(Int_t TrkLength,Float_t *TrkAry){
    Float_t SumCos2Phi=0.0, nCosPair2 = 0.0;
    for(Int_t iTrk =1; iTrk<TrkLength; iTrk++){
        for (Int_t jTrk=iTrk+1; jTrk<=TrkLength; jTrk++){
            SumCos2Phi +=cos(2.0*(TrkAry[iTrk]-TrkAry[jTrk]));
            nCosPair2++;
        }
    }
    return SumCos2Phi /=nCosPair2;
    
}

//--------------------------> Correlator 2 Particle via Like Sign
void AliAnalysisTaskEbyEChargeSep::CalCorlVia2PCLikeSign(Int_t nTrk, Float_t EvtPhiRad, Float_t EvtPhiWeightRad, Float_t V0Plane,Float_t V0APlane,Float_t V0CPlane, Float_t *TrkAry, Bool_t isCorrNorPLS){
    
    Double_t CorrCosSumls=0., CorrCosDiffls=0., CorrCosDiff2ls=0., CorrCosCrossls=0.,CorrCos3Sumls=0.;
    Double_t CorrSinSqls=0., CorrCosSumw2PCls=0., CorrCosSumw2PCWeightls=0.;
    Double_t CorrCosSumw2V0ls=0., CorrCosSumw2V0Als=0., CorrCosSumw2V0Cls=0.;
    Float_t nCrossPair2ls = 0.0, nCrossPair3ls = 0.0;
    
    
    TString PartOrTrack = "";
    if(fReadMC){
        PartOrTrack = "McReco";
        if(fMCKine)PartOrTrack = "McPart";
    }else PartOrTrack = "TrkReco";
    
    TString nOrp = "";
    if(isCorrNorPLS)nOrp += "p";
    else if(!isCorrNorPLS)nOrp += "n";
    
    TString AnaType = "";
    if(fMixedEvent)AnaType = "ME";
    else AnaType = "SE";
    
    
    for (Int_t iTr =1; iTr<= nTrk; iTr++){
        for (Int_t jTr=1; jTr< nTrk; jTr++){
            if(jTr==iTr)continue;
            if(iTr > jTr){
                CorrCosSumls     += TMath::Cos(TrkAry[iTr]+TrkAry[jTr]);
                CorrCosDiffls    += TMath::Cos(TrkAry[iTr]-TrkAry[jTr]);
                CorrCosDiff2ls   += TMath::Cos(2.0*(TrkAry[iTr]-TrkAry[jTr]));
                CorrCosCrossls   += TMath::Cos(TrkAry[iTr])*cos(TrkAry[jTr]);
                CorrSinSqls      += TMath::Sin(TrkAry[iTr])*TMath::Sin(TrkAry[jTr]);
                
                CorrCosSumw2PCls += TMath::Cos(TrkAry[iTr]+TrkAry[jTr]-2.0*EvtPhiRad);
                CorrCosSumw2PCWeightls += TMath::Cos(TrkAry[iTr]+TrkAry[jTr]-2.0*EvtPhiWeightRad);
                
                CorrCosSumw2V0ls  += TMath::Cos(TrkAry[iTr]+TrkAry[jTr]-2.0*V0Plane);
                CorrCosSumw2V0Als += TMath::Cos(TrkAry[iTr]+TrkAry[jTr]-2.0*V0APlane);
                CorrCosSumw2V0Cls += TMath::Cos(TrkAry[iTr]+TrkAry[jTr]-2.0*V0CPlane);
                nCrossPair2ls++;
            }
            
            for (Int_t kTr=jTr+1;  kTr<=nTrk; kTr++){
                if(kTr == iTr)continue;
                CorrCos3Sumls  += TMath::Cos(TrkAry[jTr]+TrkAry[kTr]-2.0*TrkAry[iTr]);
                nCrossPair3ls++;
            }
        }
    }
    
    if(nCrossPair2ls>0){
        CorrCosSumls      /=nCrossPair2ls;
        CorrCosDiffls     /=nCrossPair2ls;
        CorrCosDiff2ls    /=nCrossPair2ls;
        CorrCosCrossls    /=nCrossPair2ls;
        CorrSinSqls       /=nCrossPair2ls;
        CorrCosSumw2PCls  /=nCrossPair2ls;
        CorrCosSumw2PCWeightls /=nCrossPair2ls;
        CorrCosSumw2V0ls  /=nCrossPair2ls;
        CorrCosSumw2V0Als /=nCrossPair2ls;
        CorrCosSumw2V0Cls /=nCrossPair2ls;
    }
    
    if(nCrossPair3ls>0) {
        CorrCos3Sumls     /=nCrossPair3ls;
    }
    
    if(isCorrNorPLS)fCorrCos3LSPos = CorrCos3Sumls;
    else if(!isCorrNorPLS)fCorrCos3LSNeg = CorrCos3Sumls;
    
    ((TH2F*)fOutputList->FindObject(Form("fhCorrCos3%sLSVsCent%s_%s", nOrp.Data(), PartOrTrack.Data(), AnaType.Data())))->Fill(fEvtCentrality+0.5, CorrCos3Sumls);
    return;
}


//--------------------------> Correlator 2 Particle via Unlike Sign
void AliAnalysisTaskEbyEChargeSep::CalCorlVia2PCUnlikeSign(Int_t nTrkP, Int_t nTrkN, Float_t EvtPhiRad, Float_t EvtPhiWeightRad, Float_t V0Plane, Float_t V0APlane, Float_t V0CPlane, Float_t *TrkPArray, Float_t *TrkNArray, Bool_t isCorrNorPULS){
    
    Double_t CorrCosSumUls=0., CorrCosDiffUls=0., CorrCosCrossUls=0.,CorrCos3SumUls=0.;
    Double_t CorrSinSqUls=0., CorrCosSumw2PCUls=0., CorrCosSumw2PCWeightUls=0.;
    Double_t CorrCosSumw2V0Uls=0., CorrCosSumw2V0AUls=0., CorrCosSumw2V0CUls=0.;
    
    Int_t nCrossPair2Uls = 0, nCrossPair3Uls;
    
    TString PartOrTrack = "";
    if(fReadMC){
        PartOrTrack = "McReco";
        if(fMCKine)PartOrTrack = "McPart";
    }else PartOrTrack = "TrkReco";
    
    TString nOrp = "";
    if(isCorrNorPULS)nOrp += "p";
    else if(!isCorrNorPULS)nOrp += "n";
    
    TString AnaType = "";
    if(fMixedEvent)AnaType = "ME";
    else AnaType = "SE";
    
    for (Int_t iTr =1; iTr<=nTrkP; iTr++){
        for (Int_t jTr =1; jTr<=nTrkP; jTr++){
            for (Int_t kTr =1; kTr<=nTrkN; kTr++){
                if(iTr==1){
                    CorrCosSumUls        += TMath::Cos(TrkNArray[kTr]+TrkPArray[jTr]);
                    CorrCosDiffUls       += TMath::Cos(TrkNArray[kTr]-TrkPArray[jTr]);
                    CorrCosCrossUls      += TMath::Cos(TrkNArray[kTr])*cos(TrkPArray[jTr]);
                    
                    CorrSinSqUls         += TMath::Sin(TrkNArray[kTr])*TMath::Sin(TrkPArray[jTr]);
                    
                    CorrCosSumw2PCUls    += TMath::Cos(TrkNArray[kTr]+TrkPArray[jTr]-2.0*EvtPhiRad);
                    CorrCosSumw2PCWeightUls += TMath::Cos(TrkNArray[kTr]+TrkPArray[jTr]-2.0*EvtPhiWeightRad);
                    CorrCosSumw2V0Uls    += TMath::Cos(TrkNArray[kTr]+TrkPArray[jTr]-2.0*V0Plane);
                    CorrCosSumw2V0CUls   += TMath::Cos(TrkNArray[kTr]+TrkPArray[jTr]-2.0*V0APlane);
                    CorrCosSumw2V0CUls   += TMath::Cos(TrkNArray[kTr]+TrkPArray[jTr]-2.0*V0CPlane);
                    nCrossPair2Uls++;
                }
                if(jTr==iTr)continue;
                
                CorrCos3SumUls  += TMath::Cos(TrkPArray[jTr]+TrkNArray[kTr]-2.0*TrkPArray[iTr]);
                nCrossPair3Uls++;
            }
        }
    }
    
    if(nCrossPair2Uls>0){
        CorrCosSumUls      /=nCrossPair2Uls;
        CorrCosDiffUls     /=nCrossPair2Uls;
        CorrCosCrossUls    /=nCrossPair2Uls;
        CorrSinSqUls       /=nCrossPair2Uls;
        CorrCosSumw2PCUls  /=nCrossPair2Uls;
        CorrCosSumw2PCWeightUls /=nCrossPair2Uls;
        CorrCosSumw2V0Uls  /=nCrossPair2Uls;
        CorrCosSumw2V0AUls /=nCrossPair2Uls;
        CorrCosSumw2V0CUls /=nCrossPair2Uls;
    }
    if(CorrCos3SumUls>0){
        CorrCos3SumUls     /=nCrossPair3Uls;
    }
    
    if(isCorrNorPULS)fCorrCos3ULSPos = CorrCos3SumUls;
    else if(!isCorrNorPULS)fCorrCos3ULSNeg = CorrCos3SumUls;
    
    ((TH2F*)fOutputList->FindObject(Form("fhCorrCos3%sULSVsCent%s_%s", nOrp.Data(), PartOrTrack.Data(), AnaType.Data())))->Fill(fEvtCentrality+0.5, CorrCos3SumUls);
    return;
}
*/

//--------------------------> Get phi efficiency values
Double_t AliAnalysisTaskEbyEChargeSep::GetTrackWeight(Int_t PhiBin){
    if(!fhEffValues)return 1;
    if(fhEffValues->IsBinUnderflow(PhiBin)||fhEffValues->IsBinOverflow(PhiBin))return 1.;
    return fhEffValues->GetBinContent(PhiBin);
} 



//--------------------------> Get which AOD was used
Int_t AliAnalysisTaskEbyEChargeSep::GetFileNumberOfCurrentAOD(){
    
    fRunNumber =fAOD->GetRunNumber();
    TString  CrFileName = CurrentFileName() ;
    TString  FileNameType = "AliAOD.root";
    unsigned PosI = CrFileName.Index(fAODProdNumber)+1+fAODProdNumber.Sizeof()-1 + 1; //index=x-1 and sizeof = x+1
    unsigned PosF = (CrFileName.Sizeof()-1)-(PosI)-(FileNameType.Sizeof()-1);
    TString  AODNumberStr(CrFileName(PosI-1,PosF));
    
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
    if(TMath::Abs(VtrxX) >= 0.3)return kFALSE;
    
    Float_t VtrxY = Vtx->GetY();
    if(TMath::Abs(VtrxY) >= 0.3)return kFALSE;
    
    Float_t VtrxZ = Vtx->GetZ();
    if(TMath::Abs(VtrxZ) > fEvtabsZvtx)return kFALSE;
    
    fEventZvtx = Vtx->GetZ(); //for event Pool
    
    if(Vtx->GetType()==AliAODVertex::kPileupSPD)return kFALSE;
    
    //VZERO-AMPLITUDE
    AliAODVZERO *fAODV0 = fAOD->GetVZEROData();
    
    //Centrality Selection
    AliCentrality *EvtCent = fAOD->GetCentrality();
    fEvtCentrality  = EvtCent->GetCentralityClass10("V0M");
    if(fEvtCentrality < 1 || fEvtCentrality > 8)return kFALSE;
    
    fEvtCentOrMult = Float_t(EvtCent->GetCentralityPercentile("V0M"));
    
    float mult = fAODV0->GetMTotV0A();
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
    if (!fTrack->TestFilterBit(fBit)) return kFALSE;
    else fhTotTrack->Fill(1);
    
    if (fTrack->GetLabel()<0)return kFALSE;
    else fhTotTrack->Fill(2);
    
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

//--------------------------> Getting Correct Phi Range Tan2 like
Float_t AliAnalysisTaskEbyEChargeSep::ExtractEvtPhiInProperRange(Float_t X, Float_t Y){
    
    Float_t Temp = atan(X/Y);
    if(Y < 0 && X > 0)Temp = TMath::Pi()   + Temp;
    if(Y < 0 && X < 0)Temp = TMath::Pi()   + Temp;
    if(Y > 0 && X < 0)Temp = TMath::Pi()*2 + Temp;
    Temp /=2.;
    
    return Temp;
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



//--------------------------> Store Tracks and Events for ME: Under development
TObjArray* AliAnalysisTaskEbyEChargeSep::StoreMEEventsandTracks(TObjArray*  fMEEvtTracksIn){
    
    
    TObjArray* fMEEvtTracksOut = new TObjArray;
    fMEEvtTracksOut->SetOwner(kTRUE);
    
    //One can take >1 track per event : Large stats is reqire for ME
    Bool_t IsAlreadyFilled = kFALSE;
    Int_t iRndmPos = 0.;
    
    if(fMEEvtTracksIn->GetEntriesFast()<fMEnTrkPerEvt)fMEnTrkPerEvt=fMEEvtTracksIn->GetEntriesFast()-1;
    
    for(int iTorPStr=0; iTorPStr < fMEnTrkPerEvt; iTorPStr++){
        
        iRndmPos = gRandom->Uniform(0, fMEEvtTracksIn->GetEntriesFast()-1);
        
        if(fMCKine){
            AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fMEEvtTracksIn->At(iTorPStr));
            if (!mcPart)continue;
            //cout<<iTorPStr << ") Rndm Position --> " << iRndmPos << endl;
            fMEEvtTracksOut->Add(mcPart); // //For Mixed Event
        }else {
            AliAODTrack* track = (AliAODTrack*)fMEEvtTracksIn->At(iTorPStr);
            if (!track)continue;
            fMEEvtTracksOut->Add(track); //For Mixed Event
        }
    }
    return fMEEvtTracksOut;
}


//________________ Required track list for Mixed Events
TObjArray* AliAnalysisTaskEbyEChargeSep::CloneAcceptAndReduceTracks(TObjArray* fMEtracks)
{
    TObjArray* tracksME = new TObjArray;
    tracksME->SetOwner(kTRUE);
    
    for (Int_t i=0; i<fMEtracks->GetEntriesFast(); i++){
        AliVParticle* Part = (AliVParticle*)fMEtracks->UncheckedAt(i);
        AliMixedEvtTracks* MEtrackarray = new AliMixedEvtTracks(Part->Eta(), Part->Phi(), Part->Pt(), Part->P(), Part->Charge());
        MEtrackarray->SetUniqueID(Part->GetUniqueID());
        tracksME->Add(MEtrackarray);
    }
    
    return tracksME;
}

//Offline Mixing and filling of of Tree Event and Tracks

//____________________| Filling Associated track tree: offline
void AliAnalysisTaskEbyEChargeSep::OfflineMETree(AliAODEvent* aod, TObjArray* fTrackClnArray){
    
    
    
    for(Int_t iTrack = 0; iTrack < fTrackClnArray->GetEntriesFast(); iTrack++){
        
        ftPt_Tr = 0.;
        ftPhi_Tr = 0.;
        ftEta_Tr = 0.;
        ftCharge_Tr = 0.;
        ftMltOrCnt_Tr = 0;
        ftZvtx_Tr = 0.;
        ftPeriod_Tr = 0;
        ftOrbit_Tr = 0;
        ftBC_Tr = 0;
        
        ftMltOrCnt_Tr = fEvtCentOrMult;
        ftZvtx_Tr = fEventZvtx;
        ftPeriod_Tr = (UInt_t)aod->GetPeriodNumber();
        ftOrbit_Tr = (UInt_t)aod->GetOrbitNumber();
        ftBC_Tr = (UShort_t)aod->GetBunchCrossNumber();
        
        //cout << "ftPeriod_Tr Number = " << aod->GetPeriodNumber() << endl;
        //cout << "ftOrbit_Tr Number = " << aod->GetOrbitNumber() << endl;
        //cout << "ftBC_Tr Number = " << aod->GetBunchCrossNumber() << endl;
        
        if(fMCKine){
            AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fTrackClnArray->At(iTrack));
            if (!mcPart){
                AliError(Form("Failed casting particle from MC array!, Skipping particle for %d", iTrack));
                continue;
            }
            
            ftPt_Tr  = mcPart->Pt();
            ftPhi_Tr = mcPart->Phi();
            ftEta_Tr = mcPart->Eta();
            ftCharge_Tr = mcPart->Charge();
        }else {
            AliAODTrack* track = (AliAODTrack*)fTrackClnArray->At(iTrack);
            if (!track){
                AliError(Form("Failed casting particle from track array!, Skipping particle for %d", iTrack));
                continue;
            }
            ftPt_Tr  = track->Pt();
            ftPhi_Tr = track->Phi();
            ftEta_Tr = track->Eta();
            ftCharge_Tr = track->Charge();
        }
        fTreeTr->Fill();
    } //end of track loop
    
    PostData(3,fTreeTr);
    
    return;
}

//End

