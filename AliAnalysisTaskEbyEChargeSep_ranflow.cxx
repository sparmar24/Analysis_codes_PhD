
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
TRandom rr;

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
rphi(0),
maxp(0),
maxn(0),
minp(0),
minn(0),
dbOld(0),
dbpos(0),
dbneg(0),
asyold(0),
maxpos(0),
rpmax(0),
db90Asy025(0),
maxpos90(0),
mdb90Asy025(0),
mmaxpos90(0),
rdb90Asy025(0),
rmaxpos90(0),
cosC(0),
sinC(0),
mcosc(0),
msinc(0),
cos1(0),
cos2(0),
cos3(0),
cos4(0),
/*cos5(0),
cos6(0),
cos7(0),
cos8(0),*/
sin1(0),
sin2(0),
sin3(0),
sin4(0),
/*sin5(0),
sin6(0),
sin7(0),
sin8(0),*/
cos1p(0),
cos2p(0),
cos3p(0),
cos4p(0),
/*cos5p(0),
cos6p(0),
cos7p(0),
cos8p(0),*/
sin1p(0),
sin2p(0),
sin3p(0),
sin4p(0),
/*sin5p(0),
sin6p(0),
sin7p(0),
sin8p(0),*/
cos1n(0),
cos2n(0),
cos3n(0),
cos4n(0),
/*cos5n(0),
cos6n(0),
cos7n(0),
cos8n(0),*/
sin1n(0),
sin2n(0),
sin3n(0),
sin4n(0),
/*sin5n(0),
sin6n(0),
sin7n(0),
sin8n(0),*/
mcos1(0),
msin1(0),
mcos2(0),
msin2(0),
mcos3(0),
msin3(0),
mcos4(0),
msin4(0),
mcos1p(0),
mcos2p(0),
mcos3p(0),
mcos4p(0),
msin1p(0),
msin2p(0),
msin3p(0),
msin4p(0),
mcos1n(0),
mcos2n(0),
mcos3n(0),
mcos4n(0),
msin1n(0),
msin2n(0),
msin3n(0),
msin4n(0),
rcos1(0),
rsin1(0),
rcos2(0),
rsin2(0),
rcos3(0),
rsin3(0),
rcos4(0),
rsin4(0),
rcos1p(0),
rcos2p(0),
rcos3p(0),
rcos4p(0),
rcos1n(0),
rcos2n(0),
rcos3n(0),
rcos4n(0),
rsin1p(0),
rsin2p(0),
rsin3p(0),
rsin4p(0),
rsin1n(0),
rsin2n(0),
rsin3n(0),
rsin4n(0),
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
rphi(0),
maxp(0),
maxn(0),
minp(0),
minn(0),
dbOld(0),
dbpos(0),
dbneg(0),
asyold(0),
maxpos(0),
rpmax(0),
db90Asy025(0),
maxpos90(0),
mdb90Asy025(0),
mmaxpos90(0),
rdb90Asy025(0),
rmaxpos90(0),
cosC(0),
sinC(0),
mcosc(0),
msinc(0),
cos1(0),
cos2(0),
cos3(0),
cos4(0),
/*cos5(0),
cos6(0),
cos7(0),
cos8(0),*/
sin1(0),
sin2(0),
sin3(0),
sin4(0),
/*sin5(0),
sin6(0),
sin7(0),
sin8(0),*/
cos1p(0),
cos2p(0),
cos3p(0),
cos4p(0),
/*cos5p(0),
cos6p(0),
cos7p(0),
cos8p(0),*/
sin1p(0),
sin2p(0),
sin3p(0),
sin4p(0),
/*sin5p(0),
sin6p(0),
sin7p(0),
sin8p(0),*/
cos1n(0),
cos2n(0),
cos3n(0),
cos4n(0),
/*cos5n(0),
cos6n(0),
cos7n(0),
cos8n(0),*/
sin1n(0),
sin2n(0),
sin3n(0),
sin4n(0),
/*sin5n(0),
sin6n(0),
sin7n(0),
sin8n(0),*/
mcos1(0),
msin1(0),
mcos2(0),
msin2(0),
mcos3(0),
msin3(0),
mcos4(0),
msin4(0),
mcos1p(0),
mcos2p(0),
mcos3p(0),
mcos4p(0),
msin1p(0),
msin2p(0),
msin3p(0),
msin4p(0),
mcos1n(0),
mcos2n(0),
mcos3n(0),
mcos4n(0),
msin1n(0),
msin2n(0),
msin3n(0),
msin4n(0),
rcos1(0),
rsin1(0),
rcos2(0),
rsin2(0),
rcos3(0),
rsin3(0),
rcos4(0),
rsin4(0),
rcos1p(0),
rcos2p(0),
rcos3p(0),
rcos4p(0),
rcos1n(0),
rcos2n(0),
rcos3n(0),
rcos4n(0),
rsin1p(0),
rsin2p(0),
rsin3p(0),
rsin4p(0),
rsin1n(0),
rsin2n(0),
rsin3n(0),
rsin4n(0),
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
    fDBEvtTree->Branch("fCAODNumber", &fCAODNumber,  "fCAODNumber/I");
    fDBEvtTree->Branch("fTreeEventCounter", &fTreeEventCounter,"fTreeEventCounter/I");
    fDBEvtTree->Branch("fPeriodNo", &fPeriodNo,    "fPeriodNo/I");                
    fDBEvtTree->Branch("fOrbitNo",  &fOrbitNo,     "fOrbitNo/I");                    
    fDBEvtTree->Branch("fBunchCrossNo",&fBunchCrossNo, "fBunchCrossNo/I"); 
    fDBEvtTree->Branch("fVtrxZ",     &fVtrxZ,    "fVtrxZ/F");
    fDBEvtTree->Branch("fEvtCentrality",&fEvtCentrality,"fEvtCentrality/I"); 
    fDBEvtTree->Branch("fEvtCentOrMult",&fEvtCentOrMult,"fEvtCentOrMult/F");
    fDBEvtTree->Branch("fEvtCentCL1",&fEvtCentCL1,"fEvtCentCL1/F");
    fDBEvtTree->Branch("fEvtCentTRK",&fEvtCentTRK,"fEvtCentTRK/F");
    fDBEvtTree->Branch("fV0AMult",&fV0AMult,"fV0AMult/F");
    fDBEvtTree->Branch("fV0CMult",&fV0CMult,"fV0CMult/F");
    fDBEvtTree->Branch("fV0ACMult",&fV0ACMult,"fV0ACMult/F");
    /*    fDBEvtTree->Branch("V0AEvtAngleRad_raw",       &V0AEvtAngleRad_raw,       "V0AEvtAngleRad_raw/F");
	  fDBEvtTree->Branch("V0AEvtAngleRad_plain",     &V0AEvtAngleRad_plain,     "V0AEvtAngleRad_plain/F");*/
    fDBEvtTree->Branch("V0AEvtAngleRad_rec",       &V0AEvtAngleRad_rec,       "V0AEvtAngleRad_rec/F");
    /*    fDBEvtTree->Branch("V0AEvtAngleRad_align",     &V0AEvtAngleRad_align,     "V0AEvtAngleRad_align/F");
	  fDBEvtTree->Branch("V0AEvtAngleRad_twist",     &V0AEvtAngleRad_twist,     "V0AEvtAngleRad_twist/F");
	  fDBEvtTree->Branch("V0AEvtAngleRad_rescale",   &V0AEvtAngleRad_rescale,   "V0AEvtAngleRad_rescale/F");*/
    fDBEvtTree->Branch("V0AEvtAngleRad_latest",    &V0AEvtAngleRad_latest,    "V0AEvtAngleRad_latest/F");
    /*    fDBEvtTree->Branch("V0CEvtAngleRad_raw",       &V0CEvtAngleRad_raw,       "V0CEvtAngleRad_raw/F");
	  fDBEvtTree->Branch("V0CEvtAngleRad_plain",     &V0CEvtAngleRad_plain,     "V0CEvtAngleRad_plain/F");*/
    fDBEvtTree->Branch("V0CEvtAngleRad_rec",       &V0CEvtAngleRad_rec,       "V0CEvtAngleRad_rec/F");
    /*    fDBEvtTree->Branch("V0CEvtAngleRad_align",     &V0CEvtAngleRad_align,     "V0CEvtAngleRad_align/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_twist",     &V0CEvtAngleRad_twist,     "V0CEvtAngleRad_twist/F");
    fDBEvtTree->Branch("V0CEvtAngleRad_rescale",   &V0CEvtAngleRad_rescale,   "V0CEvtAngleRad_rescale/F");*/
    fDBEvtTree->Branch("V0CEvtAngleRad_latest",    &V0CEvtAngleRad_latest,    "V0CEvtAngleRad_latest/F");
    /*    fDBEvtTree->Branch("TPCEvtAngleRad_raw",       &TPCEvtAngleRad_raw,       "TPCEvtAngleRad_raw/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_plain",     &TPCEvtAngleRad_plain,     "TPCEvtAngleRad_plain/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_rec",       &TPCEvtAngleRad_rec,       "TPCEvtAngleRad_rec/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_align",     &TPCEvtAngleRad_align,     "TPCEvtAngleRad_align/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_twist",     &TPCEvtAngleRad_twist,     "TPCEvtAngleRad_twist/F");
    fDBEvtTree->Branch("TPCEvtAngleRad_rescale",   &TPCEvtAngleRad_rescale,   "TPCEvtAngleRad_rescale/F");*/
    fDBEvtTree->Branch("TPCEvtAngleRad_latest",    &TPCEvtAngleRad_latest,    "TPCEvtAngleRad_latest/F");
    fDBEvtTree->Branch("cosC",       &cosC,       "cosC/F");
    fDBEvtTree->Branch("sinC",       &sinC,       "sinC/F");
    fDBEvtTree->Branch("flMultiAllChrg", &flMultiAllChrg, "flMultiAllChrg/I");
    fDBEvtTree->Branch("flMultPosChrg",  &flMultPosChrg,  "flMultPosChrg/I");
    fDBEvtTree->Branch("flMultNegChrg",  &flMultNegChrg,  "flMultNegChrg/I");
    fDBEvtTree->Branch("db90Asy025",&db90Asy025,"db90Asy025/F");
    fDBEvtTree->Branch("maxpos90",&maxpos90,"maxpos90/I");
    fDBEvtTree->Branch("mdb90Asy025",&mdb90Asy025,"mdb90Asy025/F");
    fDBEvtTree->Branch("mmaxpos90",&mmaxpos90,"mmaxpos90/I");
    fDBEvtTree->Branch("rdb90Asy025",&rdb90Asy025,"rdb90Asy025/F");
    fDBEvtTree->Branch("rmaxpos90",&rmaxpos90,"rmaxpos90/I");
    fDBEvtTree->Branch("cos1",&cos1,"cos1/F");
    fDBEvtTree->Branch("cos2",&cos2,"cos2/F");
    fDBEvtTree->Branch("cos3",&cos3,"cos3/F");
    fDBEvtTree->Branch("cos4",&cos4,"cos4/F");
    /* fDBEvtTree->Branch("cos5",&cos5,"cos5/F");
    fDBEvtTree->Branch("cos6",&cos6,"cos6/F");
    fDBEvtTree->Branch("cos7",&cos7,"cos7/F");
    fDBEvtTree->Branch("cos8",&cos8,"cos8/F");*/
    fDBEvtTree->Branch("sin1",&sin1,"sin1/F");
    fDBEvtTree->Branch("sin2",&sin2,"sin2/F");
    fDBEvtTree->Branch("sin3",&sin3,"sin3/F");
    fDBEvtTree->Branch("sin4",&sin4,"sin4/F");
    /*fDBEvtTree->Branch("sin5",&sin5,"sin5/F");
    fDBEvtTree->Branch("sin6",&sin6,"sin6/F");
    fDBEvtTree->Branch("sin7",&sin7,"sin7/F");
    fDBEvtTree->Branch("sin8",&sin8,"sin8/F");*/
    fDBEvtTree->Branch("cos1p",&cos1p,"cos1p/F");
    fDBEvtTree->Branch("cos2p",&cos2p,"cos2p/F");
    fDBEvtTree->Branch("cos3p",&cos3p,"cos3p/F");
    fDBEvtTree->Branch("cos4p",&cos4p,"cos4p/F");
    /*fDBEvtTree->Branch("cos5p",&cos5p,"cos5p/F");
    fDBEvtTree->Branch("cos6p",&cos6p,"cos6p/F");
    fDBEvtTree->Branch("cos7p",&cos7p,"cos7p/F");
    fDBEvtTree->Branch("cos8p",&cos8p,"cos8p/F");*/
    fDBEvtTree->Branch("sin1p",&sin1p,"sin1p/F");
    fDBEvtTree->Branch("sin2p",&sin2p,"sin2p/F");
    fDBEvtTree->Branch("sin3p",&sin3p,"sin3p/F");
    fDBEvtTree->Branch("sin4p",&sin4p,"sin4p/F");
    /*fDBEvtTree->Branch("sin5p",&sin5p,"sin5p/F");
    fDBEvtTree->Branch("sin6p",&sin6p,"sin6p/F");
    fDBEvtTree->Branch("sin7p",&sin7p,"sin7p/F");
    fDBEvtTree->Branch("sin8p",&sin8p,"sin8p/F");*/
    fDBEvtTree->Branch("cos1n",&cos1n,"cos1n/F");
    fDBEvtTree->Branch("cos2n",&cos2n,"cos2n/F");
    fDBEvtTree->Branch("cos3n",&cos3n,"cos3n/F");
    fDBEvtTree->Branch("cos4n",&cos4n,"cos4n/F");
    /*fDBEvtTree->Branch("cos5n",&cos5n,"cos5n/F");
    fDBEvtTree->Branch("cos6n",&cos6n,"cos6n/F");
    fDBEvtTree->Branch("cos7n",&cos7n,"cos7n/F");
    fDBEvtTree->Branch("cos8n",&cos8n,"cos8n/F");*/
    fDBEvtTree->Branch("sin1n",&sin1n,"sin1n/F");
    fDBEvtTree->Branch("sin2n",&sin2n,"sin2n/F");
    fDBEvtTree->Branch("sin3n",&sin3n,"sin3n/F");
    fDBEvtTree->Branch("sin4n",&sin4n,"sin4n/F");
    /*fDBEvtTree->Branch("sin5n",&sin5n,"sin5n/F");
    fDBEvtTree->Branch("sin6n",&sin6n,"sin6n/F");
    fDBEvtTree->Branch("sin7n",&sin7n,"sin7n/F");
    fDBEvtTree->Branch("sin8n",&sin8n,"sin8n/F");*/
    //mix
    /* fDBEvtTree->Branch("mcos1", &mcos1, "mcos1/F");
    fDBEvtTree->Branch("mcos2", &mcos2, "mcos2/F");
    fDBEvtTree->Branch("mcos3", &mcos3, "mcos3/F");
    fDBEvtTree->Branch("mcos4", &mcos4, "mcos4/F");
    fDBEvtTree->Branch("msin1", &msin1, "msin1/F"); 
    fDBEvtTree->Branch("msin2", &msin2, "msin2/F");
    fDBEvtTree->Branch("msin3", &msin3, "msin3/F");
    fDBEvtTree->Branch("msin4", &msin4, "msin4/F");
    fDBEvtTree->Branch("mcos1p",&mcos1p,"mcos1p/F");
    fDBEvtTree->Branch("mcos2p",&mcos2p,"mcos2p/F");
    fDBEvtTree->Branch("mcos3p",&mcos3p,"mcos3p/F");
    fDBEvtTree->Branch("mcos4p",&mcos4p,"mcos4p/F");
    fDBEvtTree->Branch("msin1p",&msin1p,"msin1p/F");
    fDBEvtTree->Branch("msin2p",&msin2p,"msin2p/F");
    fDBEvtTree->Branch("msin3p",&msin3p,"msin3p/F");
    fDBEvtTree->Branch("msin4p",&msin4p,"msin4p/F");
    fDBEvtTree->Branch("mcos1n",&mcos1n,"mcos1n/F");
    fDBEvtTree->Branch("mcos2n",&mcos2n,"mcos2n/F");
    fDBEvtTree->Branch("mcos3n",&mcos3n,"mcos3n/F");
    fDBEvtTree->Branch("mcos4n",&mcos4n,"mcos4n/F");
    fDBEvtTree->Branch("msin1n",&msin1n,"msin1n/F");
    fDBEvtTree->Branch("msin2n",&msin2n,"msin2n/F");
    fDBEvtTree->Branch("msin3n",&msin3n,"msin3n/F");
    fDBEvtTree->Branch("msin4n",&msin4n,"msin4n/F");
*/
    // Random
    fDBEvtTree->Branch("rcos1", &rcos1, "rcos1/F");
    fDBEvtTree->Branch("rcos2", &rcos2, "rcos2/F");
    fDBEvtTree->Branch("rcos3", &rcos3, "rcos3/F");
    fDBEvtTree->Branch("rcos4", &rcos4, "rcos4/F");
    fDBEvtTree->Branch("rsin1", &rsin1, "rsin1/F");
    fDBEvtTree->Branch("rsin2", &rsin2, "rsin2/F");
    fDBEvtTree->Branch("rsin3", &rsin3, "rsin3/F");
    fDBEvtTree->Branch("rsin4", &rsin4, "rsin4/F");
    fDBEvtTree->Branch("rcos1p",&rcos1p,"rcos1p/F");
    fDBEvtTree->Branch("rcos2p",&rcos2p,"rcos2p/F");
    fDBEvtTree->Branch("rcos3p",&rcos3p,"rcos3p/F");
    fDBEvtTree->Branch("rcos4p",&rcos4p,"rcos4p/F");
    fDBEvtTree->Branch("rsin1p",&rsin1p,"rsin1p/F");
    fDBEvtTree->Branch("rsin2p",&rsin2p,"rsin2p/F");
    fDBEvtTree->Branch("rsin3p",&rsin3p,"rsin3p/F");
    fDBEvtTree->Branch("rsin4p",&rsin4p,"rsin4p/F");
    fDBEvtTree->Branch("rcos1n",&rcos1n,"rcos1n/F");
    fDBEvtTree->Branch("rcos2n",&rcos2n,"rcos2n/F");
    fDBEvtTree->Branch("rcos3n",&rcos3n,"rcos3n/F");
    fDBEvtTree->Branch("rcos4n",&rcos4n,"rcos4n/F");
    fDBEvtTree->Branch("rsin1n",&rsin1n,"rsin1n/F");
    fDBEvtTree->Branch("rsin2n",&rsin2n,"rsin2n/F");
    fDBEvtTree->Branch("rsin3n",&rsin3n,"rsin3n/F");
    fDBEvtTree->Branch("rsin4n",&rsin4n,"rsin4n/F");

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

    //cout<<fEvtCentrality<< endl;
    //float v2[8] = {0.029,0.04,0.055,0.07,0.077,0.08,0.076,0.07};
    float v2[8] = {0.0235,0.039,0.0605,0.079,0.0905,0.1135,0.116,0.115};
    float gg[8] = {0.0116,0.016,0.022,0.028,0.0308,0.032,0.0304,0.028};
    //TRandom rr;
   
    float  v2flow = rr.Gaus(v2[fEvtCentrality],gg[fEvtCentrality]);
    //cout<<" v2flow   "<<v2flow <<"    "<<fEvtCentrality<<endl;



 
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
    //cout << " fCAODNumber = " << fCAODNumber << "     "<<fCompareFileNo<<"     "<< fTreeEventCounter<< endl;
    
    if(IsEvtSctPlotReq){
      if(fRunNumPlotReq != fAOD->GetRunNumber())return;
      if(fFileNumPlotReq  != GetFileNumberOfCurrentAOD())return;
      if(fEvtNumPlotReq != fTreeEventCounter)return;
      //cout << fCAODNumber<<"    "<<fRunNumPlotReq <<"     "<< fFileNumPlotReq <<"     "<< fEvtNumPlotReq<<"    "<< fTreeEventCounter<< endl;
      // Going Back from here                                                               
    }
    
    
    //----------------|| Get corrected Event plane V0A, V0C and TPC -----||
    int myHarmonic = 2;

    //------------- V0A
    /*    const AliQnCorrectionsQnVector *myQnVectorV0_raw;
    double myEventPlaneV0_raw = 0.0;
    myQnVectorV0_raw = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","raw","raw");
    if (myQnVectorV0_raw != NULL) {
      myEventPlaneV0_raw = myQnVectorV0_raw->EventPlane(myHarmonic);
      V0AEvtAngleRad_raw = myEventPlaneV0_raw;
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_plain;
    double myEventPlaneV0_plain = 0.0;
    myQnVectorV0_plain = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","plain","plain");
    if (myQnVectorV0_plain != NULL) {
      myEventPlaneV0_plain = myQnVectorV0_plain->EventPlane(myHarmonic);
      V0AEvtAngleRad_plain = myEventPlaneV0_plain;
      }*/
    const AliQnCorrectionsQnVector *myQnVectorV0_rec;
    double myEventPlaneV0_rec = 0.0;
    myQnVectorV0_rec = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rec","rec");
    if (myQnVectorV0_rec != NULL) {
      myEventPlaneV0_rec = myQnVectorV0_rec->EventPlane(myHarmonic);
      V0AEvtAngleRad_rec = myEventPlaneV0_rec;
    }
    /*    const AliQnCorrectionsQnVector *myQnVectorV0_align;
    double myEventPlaneV0_align = 0.0;
    myQnVectorV0_align = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","align","align");
    if (myQnVectorV0_align != NULL) {
      myEventPlaneV0_align = myQnVectorV0_align->EventPlane(myHarmonic);
      V0AEvtAngleRad_align = myEventPlaneV0_align;
    }
    const AliQnCorrectionsQnVector *myQnVectorV0_twist;
    double myEventPlaneV0_twist = 0.0;
    myQnVectorV0_twist = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","twist","twist");
    if (myQnVectorV0_twist != NULL) {
      myEventPlaneV0_twist = myQnVectorV0_twist->EventPlane(myHarmonic);
      V0AEvtAngleRad_twist = myEventPlaneV0_twist;
      }
    const AliQnCorrectionsQnVector *myQnVectorV0_rescale;
    double myEventPlaneV0_rescale = 0.0;
    myQnVectorV0_rescale = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rescale","rescale");
    if (myQnVectorV0_rescale != NULL) {
      myEventPlaneV0_rescale = myQnVectorV0_rescale->EventPlane(myHarmonic);
      V0AEvtAngleRad_rescale = myEventPlaneV0_rescale;
    }*/
    const AliQnCorrectionsQnVector *myQnVectorV0_latest;
    double myEventPlaneV0_latest = 0.0;
    myQnVectorV0_latest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","latest","latest");
    if (myQnVectorV0_latest != NULL) {
      myEventPlaneV0_latest = myQnVectorV0_latest->EventPlane(myHarmonic);
      V0AEvtAngleRad_latest = myEventPlaneV0_latest;
    }
    //----------------V0C
    /* const AliQnCorrectionsQnVector *myQnVectorV0C_raw;
    double myEventPlaneV0C_raw = 0.0;
    myQnVectorV0C_raw = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","raw","raw");
    if (myQnVectorV0C_raw != NULL) {
      myEventPlaneV0C_raw = myQnVectorV0C_raw->EventPlane(myHarmonic);
      V0CEvtAngleRad_raw = myEventPlaneV0C_raw;
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_plain;
    double myEventPlaneV0C_plain = 0.0;
    myQnVectorV0C_plain = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","plain","plain");
    if (myQnVectorV0C_plain != NULL) {
      myEventPlaneV0C_plain = myQnVectorV0C_plain->EventPlane(myHarmonic);
      V0CEvtAngleRad_plain = myEventPlaneV0C_plain;
      }*/
    const AliQnCorrectionsQnVector *myQnVectorV0C_rec;
    double myEventPlaneV0C_rec = 0.0;
    myQnVectorV0C_rec = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rec","rec");
    if (myQnVectorV0C_rec != NULL) {
      myEventPlaneV0C_rec = myQnVectorV0C_rec->EventPlane(myHarmonic);
      V0CEvtAngleRad_rec = myEventPlaneV0C_rec;
    }
    /* const AliQnCorrectionsQnVector *myQnVectorV0C_align;
    double myEventPlaneV0C_align = 0.0;
    myQnVectorV0C_align = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","align","align");
    if (myQnVectorV0C_align != NULL) {
      myEventPlaneV0C_align = myQnVectorV0C_align->EventPlane(myHarmonic);
      V0CEvtAngleRad_align = myEventPlaneV0C_align;
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_twist;
    double myEventPlaneV0C_twist = 0.0;
    myQnVectorV0C_twist = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","twist","twist");
    if (myQnVectorV0C_twist != NULL) {
      myEventPlaneV0C_twist = myQnVectorV0C_twist->EventPlane(myHarmonic);
      V0CEvtAngleRad_twist = myEventPlaneV0C_twist;
    }
    const AliQnCorrectionsQnVector *myQnVectorV0C_rescale;
    double myEventPlaneV0C_rescale = 0.0;
    myQnVectorV0C_rescale = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rescale","rescale");
    if (myQnVectorV0C_rescale != NULL) {
      myEventPlaneV0C_rescale = myQnVectorV0C_rescale->EventPlane(myHarmonic);
      V0CEvtAngleRad_rescale = myEventPlaneV0C_rescale;
      }*/
    const AliQnCorrectionsQnVector *myQnVectorV0C_latest;
    double myEventPlaneV0C_latest = 0.0;
    myQnVectorV0C_latest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","latest","latest");
    if (myQnVectorV0C_latest != NULL) {
      myEventPlaneV0C_latest = myQnVectorV0C_latest->EventPlane(myHarmonic);
      V0CEvtAngleRad_latest = myEventPlaneV0C_latest;
    }

    //---------->>>TPC Histograms                                                                                       
    /*    const AliQnCorrectionsQnVector *myQnVectorTPC_raw;
    double myEventPlaneTPC_raw = 0.0;
    myQnVectorTPC_raw = fFlowQnVectorMgr->GetDetectorQnVector("TPC","raw","raw"); // Null                               
    if (myQnVectorTPC_raw != NULL) {
      myEventPlaneTPC_raw = myQnVectorTPC_raw->EventPlane(myHarmonic);
      TPCEvtAngleRad_raw = myEventPlaneTPC_raw;
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_plain;
    double myEventPlaneTPC_plain = 0.0;
    myQnVectorTPC_plain = fFlowQnVectorMgr->GetDetectorQnVector("TPC","plain","plain");
    if (myQnVectorTPC_plain != NULL) {
      myEventPlaneTPC_plain = myQnVectorTPC_plain->EventPlane(myHarmonic);
      TPCEvtAngleRad_plain = myEventPlaneTPC_plain;
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_rec;
    double myEventPlaneTPC_rec = 0.0;
    myQnVectorTPC_rec = fFlowQnVectorMgr->GetDetectorQnVector("TPC","rec","rec");
    if (myQnVectorTPC_rec != NULL) {
      myEventPlaneTPC_rec = myQnVectorTPC_rec->EventPlane(myHarmonic);
      TPCEvtAngleRad_rec = myEventPlaneTPC_rec;
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_align;
    double myEventPlaneTPC_align = 0.0;
    myQnVectorTPC_align = fFlowQnVectorMgr->GetDetectorQnVector("TPC","align","align"); // Null                         
    if (myQnVectorTPC_align != NULL) {
      myEventPlaneTPC_align = myQnVectorTPC_align->EventPlane(myHarmonic);
      TPCEvtAngleRad_align = myEventPlaneTPC_align;
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_twist;
    double myEventPlaneTPC_twist = 0.0;
    myQnVectorTPC_twist = fFlowQnVectorMgr->GetDetectorQnVector("TPC","twist","twist");
    if (myQnVectorTPC_twist != NULL) {
      myEventPlaneTPC_twist = myQnVectorTPC_twist->EventPlane(myHarmonic);
      TPCEvtAngleRad_twist = myEventPlaneTPC_twist;
    }
    const AliQnCorrectionsQnVector *myQnVectorTPC_rescale;
    double myEventPlaneTPC_rescale = 0.0;
    myQnVectorTPC_rescale = fFlowQnVectorMgr->GetDetectorQnVector("TPC","rescale","rescale"); // Null                   
    if (myQnVectorTPC_rescale != NULL) {
      myEventPlaneTPC_rescale = myQnVectorTPC_rescale->EventPlane(myHarmonic);
      TPCEvtAngleRad_rescale = myEventPlaneTPC_rescale;
      }*/
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
    Float_t qx,qy,q2x,q2y,qxp,qyp,q2xp,q2yp,q4x,q4y,q4xp,q4yp,q3xp,q5xp,q6xp,q7xp,q8xp,q3yp,q5yp,q6yp,q7yp,q8yp;
    Float_t qxn,qyn,q2xn,q2yn,q4yn,q4xn,fddb,q3xn,q5xn,q6xn,q7xn,q8xn,q3yn,q5yn,q6yn,q7yn,q8yn;
    Float_t q3x,q5x,q6x,q7x,q8x,q3y,q5y,q6y,q7y,q8y;
    Float_t cos2phi,sin2phi,cos4phi,sin4phi;
    Float_t cosch, sinch;
    Float_t phitr[5000],phitrr[2000];

    //mix charge
    Int_t chrm[5000], chr[5000], phibin =0, win;
    Float_t mqx =0.,mqy=0.,mq2x=0.0,mq2y=0.,mq4x=0.,mq4y=0.,mq3x=0.,mq3y=0.;
    Float_t mqxp =0.,mqyp=0.,mq2xp=0.0,mq2yp=0.0,mq4xp=0.,mq4yp=0.,mq3xp=0.,mq3yp=0.;
    Float_t mqxn =0.,mqyn=0.,mq2xn=0.0,mq2yn=0.0,mq4xn=0.,mq4yn=0.,mq3xn=0.,mq3yn=0.;
    Float_t mcosch=0.,msinch=0.;
    Float_t phi;
    //Random
    Float_t rqx =0.,rqy=0.,rq2x=0.0,rq2y=0.0,rq3x=0.,rq3y=0.,rq4x=0.,rq4y=0.;
    Float_t rqxp =0.,rqyp=0.,rq2xp=0.0,rq2yp=0.0,rq4xp=0.,rq4yp=0.,rq3xp=0.,rq3yp=0.;
    Float_t rqxn =0.,rqyn=0.,rq2xn=0.0,rq2yn=0.0,rq4xn=0.,rq4yn=0.,rq3xn=0.,rq3yn=0.;

	for(Int_t iPhiDegree=0; iPhiDegree<=360; iPhiDegree++){
	  NTrkAryPosChrgPhiRadBin[iPhiDegree]=0.;
	  NTrkAryNegChrgPhiRadBin[iPhiDegree]=0.;
        }
	win=0;
	qx =0.,qy=0.,q2x=0.0,q2y=0.0, q4x=0.0,q4y=0.0,q3x=0.0,q5x=0.0,q6x=0.0,q7x=0.0,q8x=0.0;
	q3y=0.0,q5y=0.0,q6y=0.0,q7y=0.0,q8y=0.0;
	qxp =0.,qyp=0.,q2xp=0.0,q2yp=0.0,q4xn=0.0,q4yn=0.0,q3xp=0.0,q5xp=0.0,q6xp=0.0,q7xp=0.0,q8xp=0.0;
	qxn =0.,qyn=0.,q2xn=0.0,q2yn=0.0,q4xp=0.0,q4yp=0.0,q3xn=0.0,q5xn=0.0,q6xn=0.0,q7xn=0.0,q8xn=0.0;
	q3yp=0.0,q5yp=0.0,q6yp=0.0,q7yp=0.0,q8yp=0.0,q3yn=0.0,q5yn=0.0,q6yn=0.0,q7yn=0.0,q8yn=0.0;
	cos2phi=0.0, sin2phi =0.0, cos4phi=0.0,sin4phi=0.0;
	cosch=0.,sinch=0.;
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
		}
		if(TorPChrg==0)continue;
		//cout<< " Charge >> " << TorPChrg <<endl;

                lMultiAllChrg++;
		phitr[lMultiAllChrg]=TorPPhi;
		chrm[lMultiAllChrg]=TorPChrg;
		chr[lMultiAllChrg]=TorPChrg;
		phitrr[lMultiAllChrg]=gRandom->Uniform(0,360);

		TrkAryAllChrgPhi[lMultiAllChrg]=TorPPhi;
		cosch += TorPChrg*cos(TorPPhi);
		sinch += TorPChrg*sin(TorPPhi);
		qx +=cos(TorPPhi);
		q2x +=cos(2.*TorPPhi);
		q3x +=cos(3.*TorPPhi);
		q4x +=cos(4.*TorPPhi);
		/*q5x +=cos(5.*TorPPhi);
		q6x +=cos(6.*TorPPhi);
		q7x +=cos(7.*TorPPhi);
		q8x +=cos(8.*TorPPhi);*/
		qy +=sin(TorPPhi);
		q2y +=sin(2.*TorPPhi);
		q3y +=sin(3.*TorPPhi);
		q4y +=sin(4.*TorPPhi);
		/*q5y +=sin(5.*TorPPhi);
		q6y +=sin(6.*TorPPhi);
		q7y +=sin(7.*TorPPhi);
		q8y +=sin(8.*TorPPhi);*/

                if(TorPChrg > 0){
		  qxp +=cos(TorPPhi);
		  q2xp +=cos(2.*TorPPhi);
		  q3xp +=cos(3.*TorPPhi);
		  q4xp +=cos(4.*TorPPhi);
		  /*q5xp +=cos(5.*TorPPhi);
                  q6xp +=cos(6.*TorPPhi);
		  q7xp +=cos(7.*TorPPhi);
		  q8xp +=cos(8.*TorPPhi);*/
		  qyp +=sin(TorPPhi);
		  q2yp +=sin(2.*TorPPhi);
		  q3yp +=sin(3.*TorPPhi);
		  q4yp +=sin(4.*TorPPhi);
                  /*q5yp +=sin(5.*TorPPhi);
                  q6yp +=sin(6.*TorPPhi);
		  q7yp +=sin(7.*TorPPhi);
                  q8yp +=sin(8.*TorPPhi);*/
		  lMultPosChrg++;
		  TrkAryPosChrgPhi[lMultPosChrg]=0.0;
		  TrkAryPosChrgPhi[lMultPosChrg]=TorPPhi;
                }else if(TorPChrg < 0){
		  qxn +=cos(TorPPhi);
		  q2xn +=cos(2.*TorPPhi);
		  q3xn +=cos(3.*TorPPhi);
		  q4xn +=cos(4.*TorPPhi);
		  /*q5xn +=cos(5.*TorPPhi);
		  q6xn +=cos(6.*TorPPhi);
		  q7xn +=cos(7.*TorPPhi);
		  q8xn +=cos(8.*TorPPhi);*/
		  qyn +=sin(TorPPhi);
		  q2yn +=sin(2.*TorPPhi);
		  q3yn +=sin(3.*TorPPhi);
		  q4yn +=sin(4.*TorPPhi);
		  /*q5yn +=sin(5.*TorPPhi);
		  q6yn +=sin(6.*TorPPhi);
		  q7yn +=sin(7.*TorPPhi);
		  q8yn +=sin(8.*TorPPhi);*/
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

	cosC=cosch;
	sinC=sinch;
	cos1 = qx;
	cos2 = q2x;
	cos3 = q3x;
	cos4 = q4x;
	/*cos5 = q5x;
	cos6 = q6x;
	cos7 = q7x;
	cos8 = q8x;*/
	sin1 =qy;
	sin2 =q2y;
	sin3 =q3y;
	sin4 =q4y;
	/*sin5 =q5y;
	sin6 =q6y;
	sin7 =q7y;
	sin8 =q8y;*/
	cos1p = qxp;
	cos2p = q2xp;
	cos3p = q3xp;
	cos4p = q4xp;
	/*cos5p = q5xp;
	cos6p = q6xp;
	cos7p = q7xp;
	cos8p = q8xp;*/
	sin1p =qyp;
	sin2p =q2yp;
	sin3p =q3yp;
	sin4p =q4yp;
	/*sin5p =q5yp;
	sin6p =q6yp;
	sin7p =q7yp;
	sin8p =q8yp;*/
	cos1n = qxn;
	cos2n = q2xn;
	cos3n = q3xn;
	cos4n = q4xn;
	/*cos5n = q5xn;
	cos6n = q6xn;
	cos7n = q7xn;
	cos8n = q8xn;*/
	sin1n =qyn;
	sin2n =q2yn;
	sin3n =q3yn;
        sin4n =q4yn;
        /*sin5n =q5yn;
        sin6n =q6yn;
        sin7n =q7yn;
        sin8n =q8yn;*/

	if(IsEvtSctPlotReq)return; // Going back if second run requested 

	//1) DB Extraction with phi window: +90Max             
	//        if(IsDBMethodOn){
	  Float_t TempTotPosChargArray[451],TempTotNegChargArray[451];
	  for(Int_t iPhi =1; iPhi<=360; iPhi++){
	    TempTotPosChargArray[iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
	    TempTotNegChargArray[iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
	    if(iPhi <=90){
	      TempTotPosChargArray[360+iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
	      TempTotNegChargArray[360+iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
	    }
	  }
	  //}
	  
	  //DBMethod random 
	  testmax(90,0.25,TempTotPosChargArray,TempTotNegChargArray);
          maxpos90=maxpos;
          db90Asy025=dbOld;
	  //cout<< "data 0.25  "<<db90Asy025<<"   "<<maxpos90 << endl;
	  
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
	    
	    mqx +=cos(phimix);
	    mq2x +=cos(2.*phimix);
	    mq3x +=cos(3.*phimix);
	    mq4x +=cos(4.*phimix);
	    mqy +=sin(phimix);
	    mq2y +=sin(2.*phimix);
	    mq3y +=sin(3.*phimix);
	    mq4y +=sin(4.*phimix);

	    if(chrm[m] >  0)
	      {
		mqxp +=cos(phimix);
		mq2xp +=cos(2.*phimix);
		mq3xp +=cos(3.*phimix);
		mq4xp +=cos(4.*phimix);
		mqyp +=sin(phimix);
		mq2yp +=sin(2.*phimix);
		mq3yp +=sin(3.*phimix);
		mq4yp +=sin(4.*phimix);
		NTrkAryPosChrgPhiRadBin[phibin]++;
	      }
	    if(chrm[m]<0 )
	      {
		mqxn +=cos(phimix);
		mq2xn +=cos(2.*phimix);
		mq3xn +=cos(3.*phimix);
		mq4xn +=cos(4.*phimix);
		mqyn +=sin(phimix);
		mq2yn +=sin(2.*phimix);
		mq3xn +=cos(3.*phimix);
		mq4yn +=sin(4.*phimix);
		NTrkAryNegChrgPhiRadBin[phibin]++;
	      }
	  } // trk loop

	mcos1 = mqx;
	mcos2 = mq2x;
	mcos3 = mq3x;
	mcos4 = mq4x;
	msin1 = mqy;
	msin2 = mq2y;
	msin3 = mq3y;
	msin4 = mq4y;
	mcos1p = mqxp;
	mcos2p = mq2xp;
	mcos3p = mq3xp;
	mcos4p = mq4xp;
	msin1p =mqyp;
	msin2p =mq2yp;
	msin3p =mq3yp;
	msin4p =mq4yp;
	mcos1n = mqxn;
	mcos2n = mq2xn;
	mcos3n = mq3xn;
	mcos4n = mq4xn;
	msin1n =mqyn;
	msin2n =mq2yn;
	msin3n =mq3yn;
	msin4n =mq4yn;
	
	for(Int_t iPhi =1; iPhi<=360; iPhi++){
	  TempTotPosChargArray[iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
	  TempTotNegChargArray[iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
	  if(iPhi <=90){
	    TempTotPosChargArray[360+iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
	    TempTotNegChargArray[360+iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
	  }
	}
	testmax(90,0.25,TempTotPosChargArray,TempTotNegChargArray);
	mmaxpos90=maxpos; 
        mdb90Asy025=dbOld;
	//cout<< "mix 0.25  "<< mdb90Asy025<<"   "<<mmaxpos90<<endl;

	
	//--------- Random Phi ------------//
	for(Int_t iPhiDegree=0; iPhiDegree<=360; iPhiDegree++){
          NTrkAryPosChrgPhiRadBin[iPhiDegree]=0.;
          NTrkAryNegChrgPhiRadBin[iPhiDegree]=0.;
        }

	Float_t phimixr;

	for (Int_t m=1;m<=lMultiAllChrg;m++)
          {
	    Float_t  phi1  =phitrr[m];
	    phimix = (phi1*TMath::Pi()/180.)-v2flow*sin(2*(phi1*TMath::Pi()/180.));
	    //phimixr = phimix-v2flow*sin(2*phimix); 
	    //cout<<"heyy  "<<fEvtCentrality<<"   "<<phimix<<"    "<<v2flow<<"   "<<phimixr<< endl;
	 
	    phibin = phimix*180./TMath::Pi()+1.0;
	    //phibin = phibin+1.0;
	    if(phibin > 360)phibin=360;
	    //cout<<"heyy  "<<fEvtCentrality<<"   "<<phimix<<"    "<<v2flow<<"   "<<phibin<< endl;

	    rqx +=cos(phimix);
	    rq2x +=cos(2.*phimix);
	    rq3x +=cos(3.*phimix);
	    rq4x +=cos(4.*phimix);
	    rqy +=sin(phimix);
	    rq2y +=sin(2.*phimix);
	    rq3y +=sin(3.*phimix);
	    rq4y +=sin(4.*phimix);

	    if(chr[m] >  0)
	      {
		NTrkAryPosChrgPhiRadBin[phibin]++;
		rqxp +=cos(phimix);
		rq2xp +=cos(2.*phimix);
		rq3xp +=cos(3.*phimix);
		rq4xp +=cos(4.*phimix);
		rqyp +=sin(phimix);
		rq2yp +=sin(2.*phimix);
		rq3yp +=sin(3.*phimix);
		rq4yp +=sin(4.*phimix);
	      } if(chr[m]<0 )
		  {
		    NTrkAryNegChrgPhiRadBin[phibin]++;
		    rqxn +=cos(phimix);
		    rq2xn +=cos(2.*phimix);
		    rq3xn +=cos(3.*phimix);
		    rq4xn +=cos(4.*phimix);
		    rqyn +=sin(phimix);
		    rq2yn +=sin(2.*phimix);
		    rq3yn +=sin(3.*phimix);
		    rq4yn +=sin(4.*phimix);
		  }
	  }

	rcos1 = rqx;
	rcos2 = rq2x;
	rcos3 = rq3x;
	rcos4 = rq4x;
	rsin1 = rqy;
	rsin2 = rq2y;
	rsin3 = rq3y;
	rsin4 = rq4y;
	rcos1p = rqxp;
	rcos2p = rq2xp;
	rcos3p = rq3xp;
	rcos4p = rq4xp;
	rsin1p =rqyp;
	rsin2p =rq2yp;
	rsin3p = rq3yp;
	rsin4p =rq4yp;
	rcos1n = rqxn;
	rcos2n = rq2xn;
	rcos3n = rq3xn;
	rcos4n = rq4xn;
	rsin1n =rqyn;
	rsin2n =rq2yn;
	rsin3n =rq3yn;
	rsin4n =rq4yn;

	for(Int_t iPhi =1; iPhi<=360; iPhi++){
          TempTotPosChargArray[iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
          TempTotNegChargArray[iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
          if(iPhi <=90){
            TempTotPosChargArray[360+iPhi]=NTrkAryPosChrgPhiRadBin[iPhi];
            TempTotNegChargArray[360+iPhi]=NTrkAryNegChrgPhiRadBin[iPhi];
          }
        }
        testmax(90,0.25,TempTotPosChargArray,TempTotNegChargArray);
        rmaxpos90=maxpos;
        rdb90Asy025=dbOld;
	//cout<< "random 0.25  "<< rdb90Asy025<<"   "<<rmaxpos90 << endl;

	 fDBEvtTree->Fill();
	 PostData(1, fOutputList);
	 PostData(2, fDBEvtTree);
	
}

//_____________________________________________________________                  
void AliAnalysisTaskEbyEChargeSep::testmax(Int_t win, Float_t asyCut, Float_t *TempTotChrgPos,Float_t *TempTotChrgNeg)
{
  Float_t temp[361],temn[361];
  Int_t  maxp,maxn,minp,minn;
  Float_t rp,rpos,rneg,rpmax,asy;

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

  rp=0.,rpos=0.,rneg=0.,rpmax=0.0,asy=0.0,dbOld=0.,dbpos=0.,dbneg=0.;  //,asyold=0.;           
  maxp=0,maxn=0,minp=0,minn=0; 

  for(Int_t k=1;k<=360;k++)
    {
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
      if(asy > asyCut)continue;
      //if(asy >0.25)continue;

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
	  dbpos=rpos;
	  dbneg=rneg;
          maxpos=k;
	  //asyold=asy;
	}
    }
  /*  asy = (float(maxp-maxn)-float(minn-minp))/(float(maxp-maxn)+float(minn-minp));
  if(maxp == 0 && maxn ==0 && minp == 0 && minn == 0)asyOld=99;
  asyOld=asy;*/

  return;
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
    if(TMath::Abs(VtrxZ) > fEvtabsZvtx)return kFALSE;
    
    fVtrxZ = VtrxZ;
    fEventZvtx = Vtx->GetZ(); //for event Pool
    
    if(Vtx->GetType()==AliAODVertex::kPileupSPD)return kFALSE;
    
    //VZERO-AMPLITUDE
    AliAODVZERO *fAODV0 = fAOD->GetVZEROData();
    
    //Centrality Selection
    AliCentrality *EvtCent = fAOD->GetCentrality();
    fEvtCentrality = EvtCent->GetCentralityClass10("V0M");
    if(fEvtCentrality < 0 || fEvtCentrality > 7)return kFALSE;
    //if(fEvtCentrality < 0)return kFALSE;

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

