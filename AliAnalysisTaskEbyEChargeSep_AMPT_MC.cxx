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
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"

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
fCentrality(0),
fMCImpactParameter(0),
fVtrxX(0),
fVtrxY(0),
fVtrxZ(0),
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
flMultiAllChrg(0),
flMultPosChrg(0),
flMultNegChrg(0),
fEvtq1Cor(0),
fEvtq1CorPos(0),
fEvtq1CorNeg(0),
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
  maxpos(0),
  fmaxpos(0),
  rpmax(0),
  maxp(0),
  maxn(0),
  minp(0),
  minn(0),
  asyOld(0),
  dbOld(0),
  ddb20(0),
  ddb30(0),
  ddb40(0),
  ddb60(0),
  ddb90(0),
  dmaxpos20(0),
  dmaxpos30(0),
  dmaxpos40(0),
  dmaxpos60(0),
  dmaxpos90(0),
  fDBValue90Old(0),
  fasy(0),
  mdb20(0),
  mdb30(0),
  mdb40(0),
  mdb60(0),
  mdb90(0),
  mmaxpos20(0),
  mmaxpos30(0),
  mmaxpos40(0),
  mmaxpos60(0),
  mmaxpos90(0),
  cosC(0),
  sinC(0),
  fmEvtq1Cor(0),
  fmEvtq1CorPos(0),
  fmEvtq1CorNeg(0),
  fmEvtq2Cor(0),
  fmEvtq2CorPos(0),
  fmEvtq2CorNeg(0),
  fmEvtq3Cor(0),
  fmEvtq3CorPos(0),
  fmEvtq3CorNeg(0),
  fmEvtq3CorUnNPP(0),
  fmEvtq3CorUnPNN(0),
  fmEvtq4Cor(0),
  fmEvtq4CorPos(0),
  fmEvtq4CorNeg(0),
  mcosC(0),
  msinC(0),
  fmDBValue90Old(0),
  fmmaxpos(0),
  cos1(0),
  sin1(0),
  cos2(0),
  sin2(0),
  cos4(0),
  sin4(0),
  cos1p(0),
  cos2p(0),
  cos4p(0),
  sin1p(0),
  sin2p(0),
  sin4p(0),
  cos1n(0),
  cos2n(0),
  cos4n(0),
  sin1n(0),
  sin2n(0),
  sin4n(0),
  mcos1p(0),
  mcos2p(0),
  mcos4p(0),
  msin1p(0),
  msin2p(0),
  msin4p(0),
  mcos1n(0),
  mcos2n(0),
  mcos4n(0),
  msin1n(0),
  msin2n(0),
  msin4n(0),
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
fCentrality(0),
fMCImpactParameter(0),
fVtrxX(0),
fVtrxY(0),
fVtrxZ(0),
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
flMultiAllChrg(0),
flMultPosChrg(0),
flMultNegChrg(0),
fEvtq1Cor(0),
fEvtq1CorPos(0),
fEvtq1CorNeg(0),
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
maxpos(0),
fmaxpos(0),
rpmax(0),
maxp(0),
maxn(0),
minp(0),
minn(0),
asyOld(0),
dbOld(0),
ddb20(0),
ddb30(0),
ddb40(0),
ddb60(0),
ddb90(0),
dmaxpos20(0),
dmaxpos30(0),
dmaxpos40(0),
dmaxpos60(0),
dmaxpos90(0),
fDBValue90Old(0),
fasy(0),
mdb20(0),
mdb30(0),
mdb40(0),
mdb60(0),
mdb90(0),
mmaxpos20(0),
mmaxpos30(0),
mmaxpos40(0),
mmaxpos60(0),
mmaxpos90(0),
cosC(0),
sinC(0),
fmEvtq1Cor(0),
fmEvtq1CorPos(0),
fmEvtq1CorNeg(0),
fmEvtq2Cor(0),
fmEvtq2CorPos(0),
fmEvtq2CorNeg(0),
fmEvtq3Cor(0),
fmEvtq3CorPos(0),
fmEvtq3CorNeg(0),
fmEvtq3CorUnNPP(0),
fmEvtq3CorUnPNN(0),
fmEvtq4Cor(0),
fmEvtq4CorPos(0),
fmEvtq4CorNeg(0),
mcosC(0),
msinC(0),
fmDBValue90Old(0),
fmmaxpos(0),
cos1(0),
sin1(0),
cos2(0),
sin2(0),
cos4(0),
sin4(0),
cos1p(0),
cos2p(0),
cos4p(0),
sin1p(0),
sin2p(0),
sin4p(0),
cos1n(0),
cos2n(0),
cos4n(0),
sin1n(0),
sin2n(0),
sin4n(0),
mcos1p(0),
mcos2p(0),
mcos4p(0),
msin1p(0),
msin2p(0),
msin4p(0),
mcos1n(0),
mcos2n(0),
mcos4n(0),
msin1n(0),
msin2n(0),
msin4n(0),
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
    fDBEvtTree->Branch("VertexX",      &fVtrxX,       "fVtrxX/F");
    fDBEvtTree->Branch("VertexY",      &fVtrxY,       "fVtrxY/F");
    fDBEvtTree->Branch("VertexZ",      &fVtrxZ,       "fVtrxZ/F");
    fDBEvtTree->Branch("EvtCentrality",&fEvtCentrality,"fEvtCentrality/F");
    fDBEvtTree->Branch("MultiAllChrg", &flMultiAllChrg, "flMultiAllChrg/I");
    fDBEvtTree->Branch("MultPosChrg",  &flMultPosChrg,  "flMultPosChrg/I");
    fDBEvtTree->Branch("MultNegChrg",  &flMultNegChrg,  "flMultNegChrg/I");
    fDBEvtTree->Branch("Q1Part",      &fEvtq1Cor,    "fEvtq1Cor/F");
    fDBEvtTree->Branch("Q1PosPart",   &fEvtq1CorPos, "fEvtq1CorPos/F");
    fDBEvtTree->Branch("Q1NegPart",   &fEvtq1CorNeg, "fEvtq1CorNeg/F");
    fDBEvtTree->Branch("Q2Part",      &fEvtq2Cor,    "fEvtq2Cor/F");
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
    fDBEvtTree->Branch("mQ1Part",     &fmEvtq1Cor,    "fmEvtq1Cor/F");
    fDBEvtTree->Branch("mQ1PosPart",  &fmEvtq1CorPos, "fmEvtq1CorPos/F");
    fDBEvtTree->Branch("mQ1NegPart",  &fmEvtq1CorNeg, "fmEvtq1CorNeg/F");
    fDBEvtTree->Branch("mQ2Part",     &fmEvtq2Cor,    "fmEvtq2Cor/F");
    fDBEvtTree->Branch("mQ2PosPart",  &fmEvtq2CorPos, "fmEvtq2CorPos/F");
    fDBEvtTree->Branch("mQ2NegPart",  &fmEvtq2CorNeg, "fmEvtq2CorNeg/F");
    fDBEvtTree->Branch("mQ3Part",     &fmEvtq3Cor,    "fmEvtq3Cor/F");
    fDBEvtTree->Branch("mQ3PosPart",  &fmEvtq3CorPos, "fmEvtq3CorPos/F");
    fDBEvtTree->Branch("mQ3NegPart",  &fmEvtq3CorNeg, "fmEvtq3CorNeg/F");
    fDBEvtTree->Branch("mQ3PartUnPNN",&fmEvtq3CorUnPNN,"fmEvtq3CorUnPNN/F");
    fDBEvtTree->Branch("mQ3PartUnNPP",&fmEvtq3CorUnNPP,"fmEvtq3CorUnNPP/F");
    fDBEvtTree->Branch("mQ4Part",     &fmEvtq4Cor,    "fmEvtq4Cor/F");
    fDBEvtTree->Branch("mQ4PartPos",  &fmEvtq4CorPos,    "fmEvtq4CorPos/F");
    fDBEvtTree->Branch("mQ4PartNeg",  &fmEvtq4CorNeg,    "fmEvtq4CorNeg/F");
    fDBEvtTree->Branch("mcosC",       &mcosC,       "mcosC/F");
    fDBEvtTree->Branch("msinC",       &msinC,       "msinC/F");
    fDBEvtTree->Branch("cos1",&cos1,"cos1/F");
    fDBEvtTree->Branch("cos2",&cos2,"cos2/F");
    fDBEvtTree->Branch("cos4",&cos4,"cos4/F");
    fDBEvtTree->Branch("sin1",&sin1,"sin1/F");
    fDBEvtTree->Branch("sin2",&sin2,"sin2/F");
    fDBEvtTree->Branch("sin4",&sin4,"sin4/F");
    fDBEvtTree->Branch("cos1p",&cos1p,"cos1p/F");
    fDBEvtTree->Branch("cos2p",&cos2p,"cos2p/F");
    fDBEvtTree->Branch("cos4p",&cos4p,"cos4p/F");
    fDBEvtTree->Branch("sin1p",&sin1p,"sin1p/F");
    fDBEvtTree->Branch("sin2p",&sin2p,"sin2p/F");
    fDBEvtTree->Branch("sin4p",&sin4p,"sin4p/F");
    fDBEvtTree->Branch("cos1n",&cos1n,"cos1n/F");
    fDBEvtTree->Branch("cos2n",&cos2n,"cos2n/F");
    fDBEvtTree->Branch("cos4n",&cos4n,"cos4n/F");
    fDBEvtTree->Branch("sin1n",&sin1n,"sin1n/F");
    fDBEvtTree->Branch("sin2n",&sin2n,"sin2n/F");
    fDBEvtTree->Branch("sin4n",&sin4n,"sin4n/F");
    fDBEvtTree->Branch("cosC",&cosC,"cosC/F");
    fDBEvtTree->Branch("sinC",&sinC,"sinC/F");
    fDBEvtTree->Branch("mcos1p",&mcos1p,"mcos1p/F");
    fDBEvtTree->Branch("mcos2p",&mcos2p,"mcos2p/F");
    fDBEvtTree->Branch("mcos4p",&mcos4p,"mcos4p/F");
    fDBEvtTree->Branch("msin1p",&msin1p,"msin1p/F");
    fDBEvtTree->Branch("msin2p",&msin2p,"msin2p/F");
    fDBEvtTree->Branch("msin4p",&msin4p,"msin4p/F");
    fDBEvtTree->Branch("mcos1n",&mcos1n,"mcos1n/F");
    fDBEvtTree->Branch("mcos2n",&mcos2n,"mcos2n/F");
    fDBEvtTree->Branch("mcos4n",&mcos4n,"mcos4n/F");
    fDBEvtTree->Branch("msin1n",&msin1n,"msin1n/F");
    fDBEvtTree->Branch("msin2n",&msin2n,"msin2n/F");
    fDBEvtTree->Branch("msin4n",&msin4n,"msin4n/F");
    fDBEvtTree->Branch("MultiAllChrg", &flMultiAllChrg, "flMultiAllChrg/I");
    fDBEvtTree->Branch("MultPosChrg", &flMultPosChrg, "flMultPosChrg/I");
    fDBEvtTree->Branch("MultNegChrg", &flMultNegChrg, "flMultNegChrg/I");
    fDBEvtTree->Branch("MaxPhi",      &fMaxPhi,       "fMaxPhi/I");
    fDBEvtTree->Branch(Form("DBValue%d", fDBWindow),&fEvtDBValue,"fEvtDBValue/F");
    fDBEvtTree->Branch(Form("RndmValue%d", fDBWindow),&fEvtRndmValue,"fEvtRndmValue/F");
    //fDBEvtTree->Branch("ddb",&fDBValue90Old,"fDBValue90Old/F");
    //fDBEvtTree->Branch("dmaxpos",&fmaxpos,"fmaxpos/I");
    fDBEvtTree->Branch("ddb90",&ddb90,"ddb90/F");
    fDBEvtTree->Branch("ddb60",&ddb60,"ddb60/F");
    fDBEvtTree->Branch("ddb40",&ddb40,"ddb40/F");
    fDBEvtTree->Branch("dmaxpos90",&dmaxpos90,"dmaxpos90/I");
    fDBEvtTree->Branch("dmaxpos60",&dmaxpos60,"dmaxpos60/I");
    fDBEvtTree->Branch("dmaxpos40",&dmaxpos40,"dmaxpos40/I");
    fDBEvtTree->Branch("mdb90",&mdb90,"mdb90/F");
    fDBEvtTree->Branch("mdb60",&mdb60,"mdb60/F");
    fDBEvtTree->Branch("mdb40",&mdb40,"mdb40/F");
    fDBEvtTree->Branch("mmaxpos90",&mmaxpos90,"mmaxpos90/I");
    fDBEvtTree->Branch("mmaxpos60",&mmaxpos60,"mmaxpos60/I");
    fDBEvtTree->Branch("mmaxpos40",&mmaxpos40,"mmaxpos40/I");

    fDBEvtTree->Branch("MaxPhi",   &fMaxPhi,        "fMaxPhi/I");
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
    if(IsCorrltrMethodOn){}
    
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
    Float_t qxn,qyn,q2xn,q2yn,q4yn,q4xn, cosch, sinch;
    Float_t phitr[5000];
    //mix charge
    Int_t chrm[5000], phibin =0, win;
    Float_t mqx =0.,mqy=0.,mq2x=0.0,mq2y=0.,mq4x=0.,mq4y=0.;
    Float_t mqxp =0.,mqyp=0.,mq2xp=0.0,mq2yp=0.0,mq4xp=0.,mq4yp=0.;
    Float_t mqxn =0.,mqyn=0.,mq2xn=0.0,mq2yn=0.0,mq4xn=0.,mq4yn=0.;
    Float_t mcosch=0.,msinch=0.;
    Float_t phi;
    
    for(Int_t iPhiDegree=0; iPhiDegree<=360; iPhiDegree++){
        NTrkAryPosChrgPhiRadBin[iPhiDegree]=0.;
        NTrkAryNegChrgPhiRadBin[iPhiDegree]=0.;
    }

    qx =0.,qy=0.,q2x=0.0,q2y=0.0, q4x=0.0,q4y=0.0;
    qxp =0.,qyp=0.,q2xp=0.0,q2yp=0.0,q4xn=0.0,q4yn=0.0;
    qxn =0.,qyn=0.,q2xn=0.0,q2yn=0.0,q4xp=0.0,q4yp=0.0;
    lMultiAllChrg = 0, lMultPosChrg = 0, lMultNegChrg = 0;
    //    cos2phi=0.0, sin2phi =0.0, cos4phi=0.0,sin4phi=0.0;
    cosch=0.,sinch=0.;
    
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
	    phitr[lMultiAllChrg]= TorPPhi;
	    chrm[lMultiAllChrg]=TorPChrg;
	    cosch += TorPChrg*cos(TorPPhi);
	    sinch += TorPChrg*sin(TorPPhi);
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
            
            
            if(IsEvtSctPlotReq && TorPP!=0){
                Float_t ScatXValue = 200.0*TorPPt*cos(TorPPhi)/TorPP;
                Float_t ScatYValue =200.0*TorPPt*sin(TorPPhi)/TorPP;
                if(TorPChrg > 0)((TH2F*)fOutputList->FindObject(Form("fhXYSctPlotPos%s_%s",PartOrTrack.Data(), AnaType.Data())))->Fill(ScatXValue,ScatYValue);
                else if (TorPChrg < 0)((TH2F*)fOutputList->FindObject(Form("fhXYSctPlotNeg%s_%s",PartOrTrack.Data(), AnaType.Data())))->Fill(ScatXValue,ScatYValue);
            }
            
	    //            Cos2Phi += cos(2.0*TorPPhi);
	    //            Sin2Phi += sin(2.0*TorPPhi);
            
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
	  
	}//trk loop
    }
    
    if(fMixedEvent){
        TObjArray* tracksMEClone = CloneAcceptAndReduceTracks(fTrackMEArrayTemp);
        fEvtPool->UpdatePool(tracksMEClone);
        if(!fEvtPool->IsReady()) return;
    }

    flMultiAllChrg =lMultiAllChrg;
    flMultPosChrg = lMultPosChrg;
    flMultNegChrg = lMultNegChrg;
    //cout<<lMultiAllChrg<<"    "<<lMultPosChrg<<"     "<<lMultNegChrg<<endl;

    ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentAllChrg_%s", AnaType.Data())))->Fill(lMultiAllChrg, fEvtCentrality+0.5);
    ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentPosChrg_%s", AnaType.Data())))->Fill(lMultPosChrg, fEvtCentrality+0.5);
    ((TH2F*)fOutputList->FindObject(Form("fhMultVsCentNegChrg_%s", AnaType.Data())))->Fill(lMultNegChrg, fEvtCentrality+0.5);

    cosC=cosch;
    sinC=sinch;

    cos1=qx;
    sin1=qy;
    cos2=q2x;
    sin2=q2y;
    cos4=q4x;
    sin4=q4y;
    cos1p=qxp;
    sin1p=qyp;
    cos2p=q2xp;
    sin2p=q2yp;
    cos4p=q4xp;
    sin4p=q4yp;
    cos1n=qxn;
    sin1n=qyn;
    cos2n=q2xn;
    sin2n=q2yn;
    cos4n=q4xn;
    sin4n=q4yn;

    // Going back if second step of RUN is requested
    if(IsEvtSctPlotReq)return;

    Float_t q1Cor    = qxp*qxn+qyp*qyn;  // opp charge  q1 <2>     
    Float_t q1CorPos = qxp*qxp+qyp*qyp-lMultPosChrg; // same q1p <2>              
    Float_t q1CorNeg = qxn*qxn+qyn*qyn-lMultNegChrg; // same q1n<2>
    fEvtq1Cor    = q1Cor;
    fEvtq1CorPos = q1CorPos;
    fEvtq1CorNeg = q1CorNeg;

    Float_t q2Cor    = Cal2PartCor(q2x,q2y,lMultiAllChrg);  // same charge 
    Float_t q2CorPos = Cal2PartCor(q2xp,q2yp,lMultPosChrg); //  
    Float_t q2CorNeg = Cal2PartCor(q2xn,q2yn,lMultNegChrg); //                 
    fEvtq2Cor    = q2Cor;
    fEvtq2CorPos = q2CorPos;
    fEvtq2CorNeg = q2CorNeg;

    //------------->3-particle (same charge)                              
    Float_t q3Cor    = Cal3PartCor(qx,qy,q2x,q2y,lMultiAllChrg);
    Float_t q3CorPos = Cal3PartCor(qxp,qyp,q2xp,q2yp,lMultPosChrg); //same ch
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
    //if(IsDBMethodOn){
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
	CalEvtPhiDumbellMethod(40,TempTotPosChargArray,TempTotNegChargArray);
	fEvtDBValue = DBPMValueMethod;
	fMaxPhi = iPhiMax;
	//fASymtry = ASymtryPar;
        ((TH2F*)fOutputList->FindObject(Form("fhDBwW%dVsCent%s_%s",fDBWindow, PartOrTrack.Data(), AnaType.Data())))->Fill(fEvtDBValue, fEvtCentrality+0.5);

        //RandomMethod
        //Float_t fRndmDBValue = CalEvtPhiDumbellMethod(fDBWindow,TempTotPosChargArray,TempTotNegChargArray, kFALSE);
        ((TH2F*)fOutputList->FindObject(Form("fhRndwW%dVsCent%s_%s",fDBWindow,PartOrTrack.Data(), AnaType.Data())))->Fill(fEvtRndmValue, fEvtCentrality+0.5);
        fEvtRndmValue = DBPMValueRndm;
	//}
    
	for(int winl=1;winl<4;winl++)
	  {
	    if(winl ==1) {
	      testmax(40,TempTotPosChargArray,TempTotNegChargArray);
	      ddb40=dbOld;
	      dmaxpos40=maxpos;
	      //cout<< "dumbell 1= "<< ddb40<<"     "<< dmaxpos40 << endl;
	    }
	    if(winl ==2) {
	      testmax(60,TempTotPosChargArray,TempTotNegChargArray);
	      ddb60=dbOld;
	      dmaxpos60=maxpos;
	      //cout<< "dumbell 2= "<< ddb60<<"     "<< dmaxpos60 << endl;
	    }
	    if(winl ==3) {
	      testmax(90,TempTotPosChargArray,TempTotNegChargArray);
	      ddb90=dbOld;
	      dmaxpos90=maxpos;
	      //cout<< "dumbell 3= "<< ddb90<<"     "<< dmaxpos90 << endl;
	    }
	  }

	    //________________random shuffling of charge___________||
	    for(Int_t iPhiDegree=0; iPhiDegree<=360; iPhiDegree++){     
	      NTrkAryPosChrgPhiRadBin[iPhiDegree]=0.;
	      NTrkAryNegChrgPhiRadBin[iPhiDegree]=0.;
	    }

	    gRandom->SetSeed();
	    srand ( time(NULL) );
	    for (int ii = 1; ii<=lMultiAllChrg; ii++)
	      {
		int j = rand() % (lMultiAllChrg+1-ii);
		j=ii+j;
		//outfile<<ii<<"  "<<j<<endl;           
		int temp = chrm[ii]; chrm[ii] = chrm[j]; chrm[j] = temp;
	      }


	    Float_t phimix;
	    for (Int_t m=1;m<=lMultiAllChrg;m++)
	      {
		phimix =phitr[m];
		if(phimix <0)phimix=-phimix;

		Float_t phi1 =phimix*180.0/TMath::Pi();
		if(phi1 <0)phi1=-phi1;
		phibin = phi1+1.0;
		if(phibin > 360)phibin=360;
		if(chrm[m] >  0)
		  {
		    mcosch +=cos(phimix);
		    msinch +=sin(phimix);
		    //    fpp[phibin] = fpp[phibin]+1;
		    mqxp +=cos(phimix);
		    mq2xp +=cos(2.*phimix);
		    mq4xp +=cos(4.*phimix);
		    mqyp +=sin(phimix);
		    mq2yp +=sin(2.*phimix);
		    mq4yp +=sin(4.*phimix);
		    NTrkAryPosChrgPhiRadBin[phibin]++;
		  }
		if(chrm[m]<0 )
		  {
		    mcosch +=-cos(phimix);
		    msinch +=-sin(phimix);
		    //    fnn[phibin] = fnn[phibin]+1;
		    mqxn +=cos(phimix);
		    mq2xn +=cos(2.*phimix);
		    mq4xn +=cos(4.*phimix);
		    mqyn +=sin(phimix);
		    mq2yn +=sin(2.*phimix);
		    mq4yn +=sin(4.*phimix);
		    NTrkAryNegChrgPhiRadBin[phibin]++;
		  }
		//cout<< m<<"   "<<chrm[m]<< endl;
	      } // trk loop

	    mcosC=mcosch;
	    msinC=msinch;

	    mcos1p=mqxp;
	    msin1p=mqyp;
	    mcos2p=mq2xp;
	    msin2p=mq2yp;
	    mcos4p=mq4xp;
	    msin4p=mq4yp;
	    mcos1n=mqxn;
	    msin1n=mqyn;
	    mcos2n=mq2xn;
	    msin2n=mq2yn;
	    mcos4n=mq4xn;
	    msin4n=mq4yn;

	    Float_t mq1Cor    = mqxp*mqxn+mqyp*mqyn; 
	    Float_t mq1CorPos = mqxp*mqxp+mqyp*mqyp-lMultPosChrg;
	    Float_t mq1CorNeg = mqxn*mqxn+mqyn*mqyn-lMultNegChrg; 
	    fmEvtq1Cor = mq1Cor;
	    fmEvtq1CorPos = mq1CorPos;
	    fmEvtq1CorNeg = mq1CorNeg;
	    Float_t mq2Cor    = Cal2PartCor(q2x,q2y,lMultiAllChrg);   
	    Float_t mq2CorPos = Cal2PartCor(q2xp,q2yp,lMultPosChrg); 
	    Float_t mq2CorNeg = Cal2PartCor(q2xn,q2yn,lMultNegChrg); 
	    fmEvtq2Cor    = mq2Cor;
	    fmEvtq2CorPos = mq2CorPos;
	    fmEvtq2CorNeg = mq2CorNeg;
	    Float_t mq3Cor    = Cal3PartCor(mqx,mqy,mq2x,mq2y,lMultiAllChrg);
	    Float_t mq3CorPos = Cal3PartCor(mqxp,mqyp,mq2xp,mq2yp,lMultPosChrg);
	    Float_t mq3CorNeg = Cal3PartCor(mqxn,mqyn,mq2xn,mq2yn,lMultNegChrg);
	    fmEvtq3Cor    = mq3Cor;
	    fmEvtq3CorPos = mq3CorPos;
	    fmEvtq3CorNeg = mq3CorNeg;
	    Float_t mq3CorPNN = Cal3PartCorUn(mqxp,mqyp,mq2xp,mq2yp,mqxn,mqyn,mq2xn,mq2yn,lMultPosChrg,lMultNegChrg); 
	    Float_t mq3CorNPP = Cal3PartCorUn(mqxn,mqyn,mq2xn,mq2yn,mqxp,mqyp,mq2xp,mq2yp,lMultNegChrg,lMultPosChrg);
	    fmEvtq3CorUnPNN = mq3CorPNN;
	    fmEvtq3CorUnNPP = mq3CorNPP;
	    Float_t mq4Cor    = Cal4PartCor(q2x,q2y,q4x,q4y,lMultiAllChrg);
	    Float_t mq4CorPos = Cal4PartCor(q2xp,q2yp,q4xp,q4yp,lMultPosChrg);
	    Float_t mq4CorNeg = Cal4PartCor(q2xn,q2yn,q4xn,q4yn,lMultNegChrg);
	    fmEvtq4Cor    = mq4Cor;
	    fmEvtq4CorPos = mq4CorPos;
	    fmEvtq4CorNeg = mq4CorNeg;

	    for(Int_t iPhi =1; iPhi<=360; iPhi++){
	      TempTotPosChargArray[iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
	      TempTotNegChargArray[iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
	      if(iPhi <=90){
		TempTotPosChargArray[360+iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
		TempTotNegChargArray[360+iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
	      }
	    }

	    for(int winl=1;winl<4;winl++)
	      {	   
		if(winl ==1) {
	      testmax(40,TempTotPosChargArray,TempTotNegChargArray);
	      mdb40=dbOld;
	      mmaxpos40=maxpos;
	      //cout<< "mix dumbell 1= "<< mdb40<<"     "<< mmaxpos40 << endl;                                                     
	    }
	    if(winl ==2)  {
	      testmax(60,TempTotPosChargArray,TempTotNegChargArray);
	      mdb60=dbOld;
	      mmaxpos60=maxpos;
	      //cout<< "mix dumbell 2= "<< mdb60<<"     "<< mmaxpos60 << endl;      
	    }
	    if(winl ==3)  {
	      testmax(90,TempTotPosChargArray,TempTotNegChargArray);
	      mdb90=dbOld;
	      mmaxpos90=maxpos;
	      //cout<< "mix dumbell 3= "<< mdb90<<"     "<< mmaxpos90 << endl;         
	    }
	  }

	// Post output data.
	fDBEvtTree->Fill();
	PostData(1, fOutputList);
	PostData(2, fDBEvtTree);
    
}
	
	//=====================================================//
	//-------------------------->Q-Cummulants--------------//
	//=====================================================//
Float_t  AliAnalysisTaskEbyEChargeSep::Cal2PartCor(Float_t cos2x,Float_t sin2y, Float_t nt)   // same charge (++, --)
{
  Float_t qc =0.;
  qc =cos2x*cos2x+sin2y*sin2y;
  //  if((nt*(nt-1.0)) > 0) qc = (qc-nt)/(nt*(nt-1.0);
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
  Float_t q4c =0.;
  float  qc  = cos2x*cos2x+sin2y*sin2y;
  float  q2c = cos4x*cos4x+sin4y*sin4y;
  float  q3c = cos4x*(cos2x*cos2x-sin2y*sin2y)+2.0*sin4y*cos2x*sin2y;
  //   float q3c = cos4x*(cos2x**2.-sin2y**2.);   
  float wt1 = (nt*(nt-1.)*(nt-2.)*(nt-3.));
  if(wt1>0) q4c = (qc*qc+q2c-2.*q3c-4.*(nt-2.)*qc)/float(wt1);
  //  if(((nt-1.)*(nt-2.)) >0) {              
  float wt = (nt-1.)*(nt-2.);
  if(wt>0) q4c +=2./float(wt);
  return q4c;
}

//--------------------------> Calculated Eventt Phi via Dumbell Method
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
    return;
}

//_____________________________________________________________
void AliAnalysisTaskEbyEChargeSep::testmax(Int_t win, Float_t *TempTotChrgPos,Float_t *TempTotChrgNeg)
{
  Float_t temp[361],temn[361];
  for(Int_t k=1;k<361;k++)
    {
      temp[k]=0.;
      temn[k]=0.;
    }

  for(Int_t k=1;k<=360;k++)
    {
      for(Int_t l=k;l<k+win;l++)
	{
	  temp[k] += TempTotChrgPos[l];
	  temn[k] += TempTotChrgNeg[l];
	}
    }

  Float_t rp,rpos=0.,rneg=0.,rpmax=0.0,asy=0.0;
  Int_t maxp=0,maxn=0,minp=0,minn=0; //asy=0.0; 

  for(Int_t k=1;k<=360;k++)
    {
      //if(k >181-win/2 && k <= 361-win/2)continue;
      rpos=0.0;
      rneg=0.0;
      if((temp[k]+temn[k])>0)  rpos=temp[k]/(temp[k]+temn[k]);
      if(k<181 && (temp[k+180]+temn[k+180])>0)  rneg=temn[k+180]/(temp[k+180]+temn[k+180]);
      if(k>180 && (temp[k-180]+temn[k-180])>0)  rneg=temn[k-180]/(temp[k-180]+temn[k-180]);

      if(rpos <0.5 || rneg <0.5)continue;
      if(rpos >=0  && rneg >=0)      rp=rpos+rneg;

      Float_t  maxpn =temp[k]-temn[k];
      Float_t  minnp;
      if(k<181) minnp =temn[k+180]-temp[k+180];
      if(k>180) minnp =temn[k-180]-temp[k-180];
      if((maxpn+minnp) ==0)asy=99;
      if((maxpn+minnp)>0)  asy = (maxpn-minnp)/(maxpn+minnp);
      if(asy <0)asy =-asy;
      if(asy >0.25)continue;

      if(rp >=rpmax)
	{
	  maxp=temp[k];
	  maxn=temn[k];
	  if(k<181)
            {
	      minp=temp[k+180];
	      minn=temn[k+180];
            }
	  if(k>180)
            {
	      minp=temp[k-180];
	      minn=temn[k-180];
	    }
	  rpmax=rp;
	  dbOld=rp;
	  maxpos=k;
	  //  asyOld=asy;
	  //  fmaxpos=k;
	  //maxpos +=win/2-1;
	  //if(maxpos >=360)maxpos=maxpos-360;
	}
    }
  /*  asy = (float(maxp-maxn)-float(minn-minp))/(float(maxp-maxn)+float(minn-minp));
  if(maxp == 0 && maxn ==0 && minp == 0 && minn == 0)asyOld=99;
  asyOld=asy; */

  return;
}


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
    
    fVtrxX = VtrxX;
    fVtrxY = VtrxY;
    fVtrxZ = VtrxZ;

    //VZERO-AMPLITUDE
    AliAODVZERO *fAODV0 = fAOD->GetVZEROData();
    
    //Centrality Selection
    /*    AliCentrality *EvtCent = fAOD->GetCentrality();
    fEvtCentrality  = EvtCent->GetCentralityClass10("V0M");
    if(fEvtCentrality < 0 || fEvtCentrality > 8)return kFALSE;
    cout<< fEvtCentrality << endl;
    fEvtCentOrMult = Float_t(EvtCent->GetCentralityPercentile("V0M"));
    */
    float mult = fAODV0->GetMTotV0A();
    //QA Plots after above cuts are applied
    ((TH1F*)fOutputList->FindObject("fhMultV0ADet"))->Fill(fAODV0->GetMTotV0A());
    ((TH1F*)fOutputList->FindObject("fhMultV0CDet"))->Fill(fAODV0->GetMTotV0C());
    ((TH1F*)fOutputList->FindObject("fhMultV0ACDet"))->Fill(fAODV0->GetMTotV0A() + fAODV0->GetMTotV0C());
    ((TH1F*)fOutputList->FindObject("fhEvtZvtx"))->Fill(VtrxZ);
    //((TH1F*)fOutputList->FindObject("fhEvtCentrality"))->Fill(fEvtCentrality+0.5);
    //((TH1F*)fOutputList->FindObject("fhEvtCentralityPer"))->Fill(Float_t(EvtCent->GetCentralityPercentile("V0M")));
    
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

