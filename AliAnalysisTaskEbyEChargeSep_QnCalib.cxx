

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
fPileSPD(0),
fEventCounter(0),
fhTotEvent(0x0),
fhTotTrack(0x0),
fhTotPart(0x0),
fmcArray(0x0),
fDBEvtTree(0),
fBit(128),
IsDBMethodOn(kTRUE),
fReadMC(kFALSE),
fMCKine(kFALSE),
fdbBin(0),
fDBWindow(90),
fEvtDBValue(0.),
fEvtmDBValue(0.),
fEvtRndmValue(0.),
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
fddb(0),
dbi(0),
fPhiRecoDgr(0),
fPhiKineDgr(0),
fIsEffCorr(0),
flMultiAllChrg(0),
flMultPosChrg(0),
flMultNegChrg(0),
fRunNumber(0),
DBPMValueMethod(0),
DBPMValueRndm(0),
iPhiMax(0),                  
fMaxPhi(0),
fASymtry(0),
ASymtry(0),
rpmax(0),
V0AEvtAngleRad_raw(0),
V0AEvtAngleRad_plain(0),
V0AEvtAngleRad_rec(0),
V0AEvtAngleRad_align(0),
V0AEvtAngleRad_twist(0),
V0AEvtAngleRad_rescale(0),
V0AEvtAngleRad_latest(0),
V0CEvtAngleRad_raw(0),
V0CEvtAngleRad_plain(0),
V0CEvtAngleRad_rec(0),
V0CEvtAngleRad_align(0),
V0CEvtAngleRad_twist(0),
V0CEvtAngleRad_rescale(0),
V0CEvtAngleRad_latest(0),
TPCEvtAngleRad_raw(0),
TPCEvtAngleRad_plain(0),
TPCEvtAngleRad_rec(0),
TPCEvtAngleRad_align(0),
TPCEvtAngleRad_twist(0),
TPCEvtAngleRad_rescale(0),
TPCEvtAngleRad_latest(0),
fExpectedCorrectionPass("latest"),
fAlternativeCorrectionPass("latest")
{
    //Default Constructor
  for (int i = 0; i < 7; i++) {
    fEvtPlaneFrameworkVZeroVsCentStep[i]=NULL;
    fEvtPlaneFrameworkVZeroCVsCentStep[i]=NULL;
    fEvtPlaneFrameworkTPCVsCentStep[i]=NULL;
  }
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
fPileSPD(0),
fEventCounter(0),
fhTotEvent(0x0),
fhTotTrack(0x0),
fhTotPart(0x0),
fmcArray(0x0),
fDBEvtTree(0),
fBit(128),
IsDBMethodOn(kTRUE),
fReadMC(kFALSE),
fMCKine(kFALSE),
fdbBin(0),
fDBWindow(90),
fEvtDBValue(0.),
fEvtmDBValue(0.),
fEvtRndmValue(0.),
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
fddb(0),
dbi(0),
fPhiRecoDgr(0),
fPhiKineDgr(0),
fIsEffCorr(0),
flMultiAllChrg(0),
flMultPosChrg(0),
flMultNegChrg(0),
fRunNumber(0),
DBPMValueMethod(0),
DBPMValueRndm(0),
iPhiMax(0),
fMaxPhi(0),
fASymtry(0),
ASymtry(0),
rpmax(0),
V0AEvtAngleRad_raw(0),
V0AEvtAngleRad_plain(0),
V0AEvtAngleRad_rec(0),
V0AEvtAngleRad_align(0),
V0AEvtAngleRad_twist(0),
V0AEvtAngleRad_rescale(0),
V0AEvtAngleRad_latest(0),
V0CEvtAngleRad_raw(0),
V0CEvtAngleRad_plain(0),
V0CEvtAngleRad_rec(0),
V0CEvtAngleRad_align(0),
V0CEvtAngleRad_twist(0),
V0CEvtAngleRad_rescale(0),
V0CEvtAngleRad_latest(0),
TPCEvtAngleRad_raw(0),
TPCEvtAngleRad_plain(0),
TPCEvtAngleRad_rec(0),
TPCEvtAngleRad_align(0),
TPCEvtAngleRad_twist(0),
TPCEvtAngleRad_rescale(0),
TPCEvtAngleRad_latest(0),
fExpectedCorrectionPass("latest"),
fAlternativeCorrectionPass("latest")
{
  for (int i = 0; i < 7; i++) {
    fEvtPlaneFrameworkVZeroVsCentStep[i]=NULL;
    fEvtPlaneFrameworkVZeroCVsCentStep[i]=NULL;
    fEvtPlaneFrameworkTPCVsCentStep[i]=NULL;
  }
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
    fDBEvtTree->Branch("RunNo",&fRunNumber,"fRunNumber/I");
    fDBEvtTree->Branch("VertexZ",     &fVtrxZ,    "fVtrxZ/F");
    fDBEvtTree->Branch("EvtCentrality",&fEvtCentrality,"fEvtCentrality/F"); 
    fDBEvtTree->Branch("V0AEvtAngleRad_raw",       &V0AEvtAngleRad_raw,       "V0AEvtAngleRad_raw/F");
    fDBEvtTree->Branch("V0AEvtAngleRad_plain",     &V0AEvtAngleRad_plain,     "V0AEvtAngleRad_plain/F");
    fDBEvtTree->Branch("V0AEvtAngleRad_rec",       &V0AEvtAngleRad_rec,       "V0AEvtAngleRad_rec/F");
    fDBEvtTree->Branch("V0AEvtAngleRad_align",     &V0AEvtAngleRad_align,     "V0AEvtAngleRad_align/F");
    fDBEvtTree->Branch("V0AEvtAngleRad_twist",     &V0AEvtAngleRad_twist,     "V0AEvtAngleRad_twist/F");
    fDBEvtTree->Branch("V0AEvtAngleRad_rescale",   &V0AEvtAngleRad_rescale,   "V0AEvtAngleRad_rescale/F");
    fDBEvtTree->Branch("V0AEvtAngleRad_latest",    &V0AEvtAngleRad_latest,    "V0AEvtAngleRad_latest/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_raw",       &V0CEvtAngleRad_raw,       "V0CEvtAngleRad_raw/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_plain",     &V0CEvtAngleRad_plain,     "V0CEvtAngleRad_plain/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_rec",       &V0CEvtAngleRad_rec,       "V0CEvtAngleRad_rec/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_align",     &V0CEvtAngleRad_align,     "V0CEvtAngleRad_align/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_twist",     &V0CEvtAngleRad_twist,     "V0CEvtAngleRad_twist/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_rescale",   &V0CEvtAngleRad_rescale,   "V0CEvtAngleRad_rescale/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_latest",    &V0CEvtAngleRad_latest,    "V0CEvtAngleRad_latest/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_raw",       &TPCEvtAngleRad_raw,       "TPCEvtAngleRad_raw/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_plain",     &TPCEvtAngleRad_plain,     "TPCEvtAngleRad_plain/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_rec",       &TPCEvtAngleRad_rec,       "TPCEvtAngleRad_rec/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_align",     &TPCEvtAngleRad_align,     "TPCEvtAngleRad_align/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_twist",     &TPCEvtAngleRad_twist,     "TPCEvtAngleRad_twist/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_rescale",   &TPCEvtAngleRad_rescale,   "TPCEvtAngleRad_rescale/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_latest",    &TPCEvtAngleRad_latest,    "TPCEvtAngleRad_latest/F");
    fDBEvtTree->Branch("MultiAllChrg", &flMultiAllChrg, "flMultiAllChrg/I");
    fDBEvtTree->Branch("MultPosChrg",  &flMultPosChrg,  "flMultPosChrg/I");
    fDBEvtTree->Branch("MultNegChrg",  &flMultNegChrg,  "flMultNegChrg/I");
    fDBEvtTree->Branch("MaxPhi",       &fMaxPhi,        "fMaxPhi/I");
    fDBEvtTree->Branch("Asy",      &fASymtry,       "fASymtry/I");
    fDBEvtTree->Branch(Form("DBValue%d", fDBWindow),&fEvtDBValue,"fEvtDBValue/F");
    fDBEvtTree->Branch(Form("RndmValue%d", fDBWindow),&fEvtRndmValue,"fEvtRndmValue/F");
    

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

    //Ztvx Distribution
    TH1F* fhEvtZvtx  = new TH1F("fhEvtZvtx","fhEvtZvtx",240, -12., 12.);
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

    //1. DB Distribution
    //Method Wise            
    if(IsDBMethodOn){
      TH2F* fhDBwWindowVsCent  = new TH2F(Form("fhDBwW%dVsCent_%s",fDBWindow, AnaType.Data()),Form("fhDBwW%dVsCent%s",fDBWindow,AnaType.Data()), 100, 0, 5, 12, -1, 11);
      fhDBwWindowVsCent->Sumw2();
      fOutputList->Add(fhDBwWindowVsCent);
      TH2F* fhRndmwWindowVsCent  = new TH2F(Form("fhRndwW%dVsCent_%s",fDBWindow, AnaType.Data()),Form("fhRndwW%dVsCent%s",fDBWindow,AnaType.Data()), 100, 0, 5, 12, -1, 11);
      fhRndmwWindowVsCent->Sumw2();
      fOutputList->Add(fhRndmwWindowVsCent);
    }

    fEvtPlaneFrameworkVZeroVsCentStep[0] = new TH2F("fEvtPlaneV0_Framework_rawCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroVsCentStep[1] = new TH2F("fEvtPlaneV0_Framework_plainCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroVsCentStep[2] = new TH2F("fEvtPlaneV0_Framework_recCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroVsCentStep[3] = new TH2F("fEvtPlaneV0_Framework_alignCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroVsCentStep[4] = new TH2F("fEvtPlaneV0_Framework_twistCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroVsCentStep[5] = new TH2F("fEvtPlaneV0_Framework_rescaleCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroVsCentStep[6] = new TH2F("fEvtPlaneV0_Framework_latestCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fOutputList->Add(fEvtPlaneFrameworkVZeroVsCentStep[0]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroVsCentStep[1]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroVsCentStep[2]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroVsCentStep[3]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroVsCentStep[4]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroVsCentStep[5]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroVsCentStep[6]);
    //----------------- V0C
    fEvtPlaneFrameworkVZeroCVsCentStep[0] = new TH2F("fEvtPlaneV0C_Framework_rawCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);  //2-d
    fEvtPlaneFrameworkVZeroCVsCentStep[1] = new TH2F("fEvtPlaneV0C_Framework_plainCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroCVsCentStep[2] = new TH2F("fEvtPlaneV0C_Framework_recCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroCVsCentStep[3] = new TH2F("fEvtPlaneV0C_Framework_alignCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroCVsCentStep[4] = new TH2F("fEvtPlaneV0C_Framework_twistCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroCVsCentStep[5] = new TH2F("fEvtPlaneV0C_Framework_rescaleCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkVZeroCVsCentStep[6] = new TH2F("fEvtPlaneV0C_Framework_latestCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fOutputList->Add(fEvtPlaneFrameworkVZeroCVsCentStep[0]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroCVsCentStep[1]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroCVsCentStep[2]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroCVsCentStep[3]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroCVsCentStep[4]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroCVsCentStep[5]);
    fOutputList->Add(fEvtPlaneFrameworkVZeroCVsCentStep[6]);
    // TPC step wise Correction Histograms                                                                              
    fEvtPlaneFrameworkTPCVsCentStep[0] = new TH2F("fEvtPlaneTPC_Framework_rawCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);  //2-d
    fEvtPlaneFrameworkTPCVsCentStep[1] = new TH2F("fEvtPlaneTPC_Framework_plainCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkTPCVsCentStep[2] = new TH2F("fEvtPlaneTPC_Framework_recCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkTPCVsCentStep[3] = new TH2F("fEvtPlaneTPC_Framework_alignCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkTPCVsCentStep[4] = new TH2F("fEvtPlaneTPC_Framework_twistCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkTPCVsCentStep[5] = new TH2F("fEvtPlaneTPC_Framework_rescaleCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fEvtPlaneFrameworkTPCVsCentStep[6] = new TH2F("fEvtPlaneTPC_Framework_latestCent", "", 400, -2*TMath::Pi(), 2*TMath::Pi(), 12, 0, 12);
    fOutputList->Add(fEvtPlaneFrameworkTPCVsCentStep[0]);
    fOutputList->Add(fEvtPlaneFrameworkTPCVsCentStep[1]);
    fOutputList->Add(fEvtPlaneFrameworkTPCVsCentStep[2]);
    fOutputList->Add(fEvtPlaneFrameworkTPCVsCentStep[3]);
    fOutputList->Add(fEvtPlaneFrameworkTPCVsCentStep[4]);
    fOutputList->Add(fEvtPlaneFrameworkTPCVsCentStep[5]);
    fOutputList->Add(fEvtPlaneFrameworkTPCVsCentStep[6]);

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
    fRunNumber = 0;
    fRunNumber = (Int_t)fAOD->GetRunNumber();

    int myHarmonic = 2;

    //------------- V0A
    const AliQnCorrectionsQnVector *myQnVectorV0_raw;
    double myEventPlaneV0_raw = 0.0;
    myQnVectorV0_raw = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","raw","raw");
    if (myQnVectorV0_raw != NULL) {
      myEventPlaneV0_raw = myQnVectorV0_raw->EventPlane(myHarmonic);
      V0AEvtAngleRad_raw = myEventPlaneV0_raw;
      fEvtPlaneFrameworkVZeroVsCentStep[0] ->Fill(myEventPlaneV0_raw, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_plain;
    double myEventPlaneV0_plain = 0.0;
    myQnVectorV0_plain = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","plain","plain");
    if (myQnVectorV0_plain != NULL) {
      myEventPlaneV0_plain = myQnVectorV0_plain->EventPlane(myHarmonic);
      V0AEvtAngleRad_plain = myEventPlaneV0_plain;
      fEvtPlaneFrameworkVZeroVsCentStep[1] ->Fill(myEventPlaneV0_plain, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_rec;
    double myEventPlaneV0_rec = 0.0;
    myQnVectorV0_rec = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rec","rec");
    if (myQnVectorV0_rec != NULL) {
      myEventPlaneV0_rec = myQnVectorV0_rec->EventPlane(myHarmonic);
      V0AEvtAngleRad_rec = myEventPlaneV0_rec;
      fEvtPlaneFrameworkVZeroVsCentStep[2] ->Fill(myEventPlaneV0_rec, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_align;
    double myEventPlaneV0_align = 0.0;
    myQnVectorV0_align = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","align","align");
    if (myQnVectorV0_align != NULL) {
      myEventPlaneV0_align = myQnVectorV0_align->EventPlane(myHarmonic);
      V0AEvtAngleRad_align = myEventPlaneV0_align;
      fEvtPlaneFrameworkVZeroVsCentStep[3] ->Fill(myEventPlaneV0_align, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_twist;
    double myEventPlaneV0_twist = 0.0;
    myQnVectorV0_twist = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","twist","twist");
    if (myQnVectorV0_twist != NULL) {
      myEventPlaneV0_twist = myQnVectorV0_twist->EventPlane(myHarmonic);
      V0AEvtAngleRad_twist = myEventPlaneV0_twist;
      fEvtPlaneFrameworkVZeroVsCentStep[4] ->Fill(myEventPlaneV0_twist, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_rescale;
    double myEventPlaneV0_rescale = 0.0;
    myQnVectorV0_rescale = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rescale","rescale");
    if (myQnVectorV0_rescale != NULL) {
      myEventPlaneV0_rescale = myQnVectorV0_rescale->EventPlane(myHarmonic);
      V0AEvtAngleRad_rescale = myEventPlaneV0_rescale;
      fEvtPlaneFrameworkVZeroVsCentStep[5] ->Fill(myEventPlaneV0_rescale, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_latest;
    double myEventPlaneV0_latest = 0.0;
    myQnVectorV0_latest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","latest","latest");
    if (myQnVectorV0_latest != NULL) {
      myEventPlaneV0_latest = myQnVectorV0_latest->EventPlane(myHarmonic);
      V0AEvtAngleRad_latest = myEventPlaneV0_latest;
      fEvtPlaneFrameworkVZeroVsCentStep[6] ->Fill(myEventPlaneV0_latest, fEvtCentrality+0.5);
    }
    //----------------V0C
    const AliQnCorrectionsQnVector *myQnVectorV0C_raw;
    double myEventPlaneV0C_raw = 0.0;
    myQnVectorV0C_raw = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","raw","raw");
    if (myQnVectorV0C_raw != NULL) {
      myEventPlaneV0C_raw = myQnVectorV0C_raw->EventPlane(myHarmonic);
      V0CEvtAngleRad_raw = myEventPlaneV0C_raw;
      fEvtPlaneFrameworkVZeroCVsCentStep[0]->Fill(myEventPlaneV0C_raw, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_plain;
    double myEventPlaneV0C_plain = 0.0;
    myQnVectorV0C_plain = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","plain","plain");
    if (myQnVectorV0C_plain != NULL) {
      myEventPlaneV0C_plain = myQnVectorV0C_plain->EventPlane(myHarmonic);
      V0CEvtAngleRad_plain = myEventPlaneV0C_plain;
      fEvtPlaneFrameworkVZeroCVsCentStep[1] ->Fill(myEventPlaneV0C_plain, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_rec;
    double myEventPlaneV0C_rec = 0.0;
    myQnVectorV0C_rec = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rec","rec");
    if (myQnVectorV0C_rec != NULL) {
      myEventPlaneV0C_rec = myQnVectorV0C_rec->EventPlane(myHarmonic);
      V0CEvtAngleRad_rec = myEventPlaneV0C_rec;
      fEvtPlaneFrameworkVZeroCVsCentStep[2] ->Fill(myEventPlaneV0C_rec, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_align;
    double myEventPlaneV0C_align = 0.0;
    myQnVectorV0C_align = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","align","align");
    if (myQnVectorV0C_align != NULL) {
      myEventPlaneV0C_align = myQnVectorV0C_align->EventPlane(myHarmonic);
      V0CEvtAngleRad_align = myEventPlaneV0C_align;
      fEvtPlaneFrameworkVZeroCVsCentStep[3] ->Fill(myEventPlaneV0C_align, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_twist;
    double myEventPlaneV0C_twist = 0.0;
    myQnVectorV0C_twist = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","twist","twist");
    if (myQnVectorV0C_twist != NULL) {
      myEventPlaneV0C_twist = myQnVectorV0C_twist->EventPlane(myHarmonic);
      V0CEvtAngleRad_twist = myEventPlaneV0C_twist;
      fEvtPlaneFrameworkVZeroCVsCentStep[4] ->Fill(myEventPlaneV0C_twist, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_rescale;
    double myEventPlaneV0C_rescale = 0.0;
    myQnVectorV0C_rescale = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rescale","rescale");
    if (myQnVectorV0C_rescale != NULL) {
      myEventPlaneV0C_rescale = myQnVectorV0C_rescale->EventPlane(myHarmonic);
      V0CEvtAngleRad_rescale = myEventPlaneV0C_rescale;
      fEvtPlaneFrameworkVZeroCVsCentStep[5] ->Fill(myEventPlaneV0C_rescale, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_latest;
    double myEventPlaneV0C_latest = 0.0;
    myQnVectorV0C_latest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","latest","latest");
    if (myQnVectorV0C_latest != NULL) {
      myEventPlaneV0C_latest = myQnVectorV0C_latest->EventPlane(myHarmonic);
      V0CEvtAngleRad_latest = myEventPlaneV0C_latest;
      fEvtPlaneFrameworkVZeroCVsCentStep[6] ->Fill(myEventPlaneV0C_latest, fEvtCentrality+0.5);
    }

    //---------->>>TPC Histograms                                                                                       
    const AliQnCorrectionsQnVector *myQnVectorTPC_raw;
    double myEventPlaneTPC_raw = 0.0;
    myQnVectorTPC_raw = fFlowQnVectorMgr->GetDetectorQnVector("TPC","raw","raw"); // Null                               
    if (myQnVectorTPC_raw != NULL) {
      myEventPlaneTPC_raw = myQnVectorTPC_raw->EventPlane(myHarmonic);
      TPCEvtAngleRad_raw = myEventPlaneTPC_raw;
      fEvtPlaneFrameworkTPCVsCentStep[0]->Fill(myEventPlaneTPC_raw, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_plain;
    double myEventPlaneTPC_plain = 0.0;
    myQnVectorTPC_plain = fFlowQnVectorMgr->GetDetectorQnVector("TPC","plain","plain");
    if (myQnVectorTPC_plain != NULL) {
      myEventPlaneTPC_plain = myQnVectorTPC_plain->EventPlane(myHarmonic);
      TPCEvtAngleRad_plain = myEventPlaneTPC_plain;
      fEvtPlaneFrameworkTPCVsCentStep[1]->Fill(myEventPlaneTPC_plain, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_rec;
    double myEventPlaneTPC_rec = 0.0;
    myQnVectorTPC_rec = fFlowQnVectorMgr->GetDetectorQnVector("TPC","rec","rec");
    if (myQnVectorTPC_rec != NULL) {
      myEventPlaneTPC_rec = myQnVectorTPC_rec->EventPlane(myHarmonic);
      TPCEvtAngleRad_rec = myEventPlaneTPC_rec;
      fEvtPlaneFrameworkTPCVsCentStep[2]->Fill(myEventPlaneTPC_rec, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_align;
    double myEventPlaneTPC_align = 0.0;
    myQnVectorTPC_align = fFlowQnVectorMgr->GetDetectorQnVector("TPC","align","align"); // Null                         
    if (myQnVectorTPC_align != NULL) {
      myEventPlaneTPC_align = myQnVectorTPC_align->EventPlane(myHarmonic);
      TPCEvtAngleRad_align = myEventPlaneTPC_align;
      fEvtPlaneFrameworkTPCVsCentStep[3]->Fill(myEventPlaneTPC_align, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_twist;
    double myEventPlaneTPC_twist = 0.0;
    myQnVectorTPC_twist = fFlowQnVectorMgr->GetDetectorQnVector("TPC","twist","twist");
    if (myQnVectorTPC_twist != NULL) {
      myEventPlaneTPC_twist = myQnVectorTPC_twist->EventPlane(myHarmonic);
      TPCEvtAngleRad_twist = myEventPlaneTPC_twist;
      fEvtPlaneFrameworkTPCVsCentStep[4]->Fill(myEventPlaneTPC_twist, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_rescale;
    double myEventPlaneTPC_rescale = 0.0;
    myQnVectorTPC_rescale = fFlowQnVectorMgr->GetDetectorQnVector("TPC","rescale","rescale"); // Null                   
    if (myQnVectorTPC_rescale != NULL) {
      myEventPlaneTPC_rescale = myQnVectorTPC_rescale->EventPlane(myHarmonic);
      TPCEvtAngleRad_rescale = myEventPlaneTPC_rescale;
      fEvtPlaneFrameworkTPCVsCentStep[5]->Fill(myEventPlaneTPC_rescale, fEvtCentrality+0.5);
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_latest;
    double myEventPlaneTPC_latest = 0.0;
    myQnVectorTPC_latest = fFlowQnVectorMgr->GetDetectorQnVector("TPC","latest","latest");
    if (myQnVectorTPC_latest != NULL) {
      myEventPlaneTPC_latest = myQnVectorTPC_latest->EventPlane(myHarmonic);
      TPCEvtAngleRad_latest = myEventPlaneTPC_latest;
      fEvtPlaneFrameworkTPCVsCentStep[6]->Fill(myEventPlaneTPC_latest, fEvtCentrality+0.5);
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


	for(Int_t iPhiDegree=0; iPhiDegree<=360; iPhiDegree++){
	  NTrkAryPosChrgPhiRadBin[iPhiDegree]=0.;
	  NTrkAryNegChrgPhiRadBin[iPhiDegree]=0.;
        }

        lMultiAllChrg = 0, lMultPosChrg = 0, lMultNegChrg = 0;
              
	  Int_t nTrackOrPart = 0; 
	  if(fMCKine)nTrackOrPart = fmcArray->GetEntriesFast(); // kine                                                             
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
		}

		if(TorPChrg==0)continue;
		//cout<< " Charge >> " << TorPChrg <<endl;

                lMultiAllChrg++;
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

	flMultiAllChrg =lMultiAllChrg;
        flMultPosChrg = lMultPosChrg;
        flMultNegChrg = lMultNegChrg;
	//cout<< lMultiAllChrg <<"    "<<lMultPosChrg<<"    "<<lMultNegChrg<< endl;


        ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentAllChrg_%s", AnaType.Data())))->Fill(lMultiAllChrg, fEvtCentrality+0.5);
        ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentPosChrg_%s", AnaType.Data())))->Fill(lMultPosChrg, fEvtCentrality+0.5);
        ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentNegChrg_%s", AnaType.Data())))->Fill(lMultNegChrg, fEvtCentrality+0.5);



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
	}
	  

	 fDBEvtTree->Fill();
	 PostData(1, fOutputList);
	 PostData(2, fDBEvtTree);
	
}


//--------------------------> Calculated Eventt Phi via Dumbell Method    
void AliAnalysisTaskEbyEChargeSep::CalEvtPhiDumbellMethod(Int_t fPhiWindow, Float_t *TempTotChrgPos,Float_t *TempTotChrgNeg){

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

  Float_t RefFracAB=0.,FracAB=0.,FracAPosFront=0.,FracBNegBack=0; //, ASymtry=999.;
  //  Float_t DBPMValueMethod = 0.0,
  //Float_t DBPMValueRndm = -9.0;
  //  Int_t iPhiMax = 0, iPhiPosFront180=0, iPhiNegBack180=0;
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
    FracBNegBack  = SumPCInWindow[iPhiNegBack180]/(SumPCInWindow[iPhiNegBack180]+SumNCInWindow[iPhiNegBack180]);

    if(FracAPosFront < 0.5 || FracBNegBack < 0.5) continue;
    if((FracAPosFront+FracBNegBack)<=0)continue; //protection   
                                                           
    ASymtry = (FracAPosFront-FracBNegBack)/(FracAPosFront+FracBNegBack);
    if(ASymtry <0) ASymtry = -ASymtry;
    if(ASymtry >0.25) continue;

    // NEW July 15 2017                                                                           
    //        if(TMath::Abs(ASymtry > 0.25)) continue;

    FracAB = FracAPosFront + FracBNegBack;
    if(FracAB > RefFracAB){ //leading sum                                                                                                                                              
      RefFracAB=FracAB;
      DBPMValueMethod=FracAB;
      iPhiMax=iPhi;
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
    FracBNegBackRndm  = SumPCInWindow[iPhiRndmNegBack]/(SumPCInWindow[iPhiRndmNegBack]+SumNCInWindow[iPhiRndmNegBack]);
  }

  FracABRndm = FracAPosFrontRndm + FracBNegBackRndm;
  if(FracABRndm>0)DBPMValueRndm = FracABRndm;

  //  DBPMValueRndm=FracABRndm;

  /*  if(isMethod)return DBPMValueMethod;
      else return DBPMValueRndm;*/
}


//--------------------------> Get phi efficiency values                          
Double_t AliAnalysisTaskEbyEChargeSep::GetTrackWeight(Int_t PtBin){
  if(!fhEffValues)return 1;
  if(fhEffValues->IsBinUnderflow(PtBin)||fhEffValues->IsBinOverflow(PtBin))return 1.;
  return fhEffValues->GetBinContent(PtBin);
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
    if(TMath::Abs(VtrxX) >= 1.0)return kFALSE;  // 0.3
    
    Float_t VtrxY = Vtx->GetY();
    if(TMath::Abs(VtrxY) >= 1.0)return kFALSE;
    
    Float_t VtrxZ = Vtx->GetZ();
    if(TMath::Abs(VtrxZ) > fEvtabsZvtx)return kFALSE;
    

    fVtrxZ = VtrxZ;
    fEventZvtx = Vtx->GetZ(); //for event Pool
    
    if(Vtx->GetType()==AliAODVertex::kPileupSPD)return kFALSE;
    
    //VZERO-AMPLITUDE
    AliAODVZERO *fAODV0 = fAOD->GetVZEROData();
    
    //Centrality Selection
    AliCentrality *EvtCent = fAOD->GetCentrality();
    fEvtCentrality  = EvtCent->GetCentralityClass10("V0M");
    if(fEvtCentrality < 1 || fEvtCentrality > 8)return kFALSE;
    
    fEvtCentOrMult = Float_t(EvtCent->GetCentralityPercentile("V0M"));


    //QA Plots after above cuts are applied
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
    if(chi2ndf > 4.0)return kFALSE;
    
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

