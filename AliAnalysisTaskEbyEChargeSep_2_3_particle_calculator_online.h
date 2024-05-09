
#ifndef AliAnalysisTaskEbyEChargeSep_H
#define AliAnalysisTaskEbyEChargeSep_H

class TH1F;
class TH1D;
class TH2F;
class TTree;
class NTuple;
class TString;
class TParticle ;
class TObjArray;
class TClonesArray ;
class TDatabasePDG;
class TParticlePDG;
class AliInputEventHandler;
class AliAnalysisTaskSE;
class AliEventPoolManager;
class AliAODEvent;
class AliAODHandler;
class AliAODInputHandler;
class AliMCEvent;
class AliVParticle;
class AliMCParticle;
class AliStack;

class AliMixedEvtTracks;

#include <TClonesArray.h>
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODVertex.h"
#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskEbyEChargeSep : public AliAnalysisTaskSE{
    
public:
    
    AliAnalysisTaskEbyEChargeSep();
    AliAnalysisTaskEbyEChargeSep(const char *name);
    virtual ~AliAnalysisTaskEbyEChargeSep();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    Int_t GetFileNumberOfCurrentAOD();
    Bool_t IsTrackSelectionPassed(AliAODTrack *fTrack);
    Bool_t IsGeneratedPartPassed(AliAODMCParticle* mcPart);
    Bool_t IsEventSelectionPassed();
    
    //  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    //    fESDtrackCuts = trackCuts;}
    
    void SetCentralityEstimator(const char* CentEst){fCentEstimator = CentEst;}
    void SetAODProdcNumber(const char* AODnum){fAODProdNumber = AODnum;}
    void SetFilterBitNumber(Int_t Bit){fBit = Bit;}
    void SetSwitchToDBCalc(Bool_t Sel, Int_t wnd){IsDBMethodOn = Sel; fDBWindow = wnd;}
    void SetSwitchToCorrltrCalc(Bool_t Sel1,Bool_t Sel2, Bool_t Sel3){IsCorrltrMethodOn=Sel1; fCorrltrLS=Sel2; fCorrltrULS=Sel3;}
    
    void SetEvtZvtxValue(Double_t Sel){fEvtabsZvtx = Sel;}
    void SetTrkPtRange(Double_t Sel1, Double_t Sel2){fTrkPtMin = Sel1; fTrkPtMax = Sel2;}
    void SetTrkEtaRange(Double_t Sel1, Double_t Sel2){fTrkEtaMin = Sel1; fTrkEtaMax = Sel2;}
    void SetTrkabsDCAXY(Double_t Sel){fTrkabsDCAxy = Sel;}
    void SetTrkabsDCAZ(Double_t Sel){fTrkabsDCAz = Sel;}
    void SetTrkMinClusterTPC(Int_t Sel){fTrkTPCclus = Sel;}
    void SetMCAnalysis(Bool_t Sel1, Bool_t Sel2){fReadMC = Sel1, fMCKine=Sel2;}
    void SetMEAnalysis(Bool_t Sel1){fMixedEvent = Sel1;}
    void SetMEEvtTrkCriteria(Int_t Val1,Int_t Val2){fMaxEvtMixed = Val1; fMaxTrkDepth = Val2;}
    void SetMEnTrkPerEvt(Int_t Val1){fMEnTrkPerEvt = Val1;}
    
    void SetMEEvtTrkTREE(Bool_t Sel){fFillMETree = Sel;}
    
    
    void SetRunFile4ScatPlots(Bool_t YN, Int_t Sel0, Int_t Sel1, Int_t Sel2){
        IsEvtSctPlotReq=YN,
        fRunNumPlotReq  = Sel0,
        fFileNumPlotReq = Sel1,
        fEvtNumPlotReq  = Sel2;
    }
    
    void SetEfficiencyCorrection(Bool_t Sel1){fIsEffCorr = Sel1;}
    void SetInputEfficiencyFile(TH1D *hMap){if(fhEffValues)delete fhEffValues;fhEffValues=(TH1D*)hMap->Clone();}
    TObjArray* StoreMEEventsandTracks(TObjArray*  fSEMEEvtTracks);

    //Float_t CalEvtPhiDumbellMethod(Int_t fPhiWindow, Float_t *TempTotChargePos,Float_t *TempTotChargeNeg, Bool_t isMethod);
    void CalEvtPhiDumbellMethod(Int_t fPhiWindow, Float_t *TempTotChargePos,Float_t *TempTotChargeNeg);
    Float_t ExtractEvtPhiInProperRange(Float_t X, Float_t Y);
    Float_t Cal3PartCor(float x,float y,float x1,float y1, Int_t num1);
    Float_t Cal2PartCor(float x2,float y2, Float_t num2);
    Float_t Cal3PartCorUn(float xp,float yp,float x2p,float y2p, float xn,float yn,float x2n,float y2n, Int_t num1, Int_t num2);
    Float_t Cal4PartCor(float cos2x,float sin2y,float cos4x,float sin4y, Int_t num1);

    /*    Float_t CalCorlVia2PCLikeSignwEvtP(Int_t num,Float_t *array);
    
    void CalCorlVia2PCLikeSign(Int_t nTrk, Float_t EvtPhiRad, Float_t EvtPhiWeightRad, Float_t V0Plane,Float_t V0APlane,Float_t V0CPlane, Float_t *TrkArray, Bool_t isCorrNorPLS);
    void CalCorlVia2PCUnlikeSign(Int_t nTrkP, Int_t nTrkN, Float_t EvtPhiRad,Float_t EvtPhiWeightRad, Float_t V0Plane, Float_t V0APlane, Float_t V0CPlane, Float_t *TrkPArray, Float_t *TrkNArray, Bool_t isCorrNorPULS);
    */
    void OfflineMETree(AliAODEvent* aod, TObjArray* fTrackClnArray);
    
    // Getters
    Double_t GetTrackWeight(Int_t PhiBin);
    
    //Booliens
    Bool_t IsTrackEffMapLoaded(){if(fhEffValues) return kTRUE; else return kFALSE;}
    
private:
    
    AliAnalysisTaskEbyEChargeSep(const AliAnalysisTaskEbyEChargeSep&);
    AliAnalysisTaskEbyEChargeSep& operator=(const AliAnalysisTaskEbyEChargeSep&);
    TObjArray* CloneAcceptAndReduceTracks(TObjArray* fMEtracks);

    
    AliAODEvent *fAOD;  //! AOD Event object
    TList *fOutputList; //! Output list
    TClonesArray *farrayMC;//!
    Float_t fEvtCentrality; //
    Float_t fEvtCentOrMult; //
    TString fCentEstimator;
    Int_t fEventCounter; //
    TH1I* fhTotEvent;//!
    TH1I* fhTotTrack;//!
    TH1I* fhTotPart;//!
    TClonesArray* fmcArray;//!
    Bool_t IsEvtSctPlotReq;
    Int_t fRunNumPlotReq;
    Int_t fFileNumPlotReq;
    Int_t fEvtNumPlotReq;
    TString fAODProdNumber;
    Int_t fCAODNumber;
    Int_t fRunNumber;
    TTree *fDBEvtTree; //!
    TTree *fTreeTr;//!
    Int_t fBit;
    Bool_t IsDBMethodOn;
    Bool_t IsCorrltrMethodOn;
    Bool_t fCorrltrLS;
    Bool_t fCorrltrULS;
    Bool_t fReadMC;
    Bool_t fMCKine;
    Int_t fDBWindow;
    Float_t fEvtDBValue;
    Float_t fEvtRndmValue;
    Float_t DBPMValueMethod;
    Float_t DBPMValueRndm;
    Int_t   iPhiMax;
    Int_t   fMaxPhi;
    Float_t fASymtry;
    Float_t ASymtryPar;
    Float_t fEvtq2Cor;
    Float_t fEvtq2CorPos;
    Float_t fEvtq2CorNeg;
    Float_t fEvtq3Cor;
    Float_t fEvtq3CorPos;
    Float_t fEvtq3CorNeg;
    Float_t fEvtq3CorUnNPP;
    Float_t fEvtq3CorUnPNN;
    Float_t fEvtq4Cor;
    Float_t fEvtq4CorPos;
    Float_t fEvtq4CorNeg;
    Int_t fCompareFileNo;
    Int_t fTreeEventCounter;
    Double_t fEvtabsZvtx;
    Double_t fEventZvtx;
    Double_t fTrkPtMin;
    Double_t fTrkPtMax;
    Double_t fTrkEtaMin;
    Double_t fTrkEtaMax;
    Double_t fTrkabsDCAxy;
    Double_t fTrkabsDCAz;
    UInt_t fTrkTPCclus;
    Float_t fCorrCos3LSPos;
    Float_t fCorrCos3LSNeg;
    Float_t fCorrCos3ULSPos;
    Float_t fCorrCos3ULSNeg;
    Int_t fPhiRecoDgr;
    Int_t fPhiKineDgr;
    TH1D  *fhEffValues; 
    Bool_t fIsEffCorr;
    AliEventPoolManager *fPoolMngr;//!
    AliEventPool *fEvtPool;//!
    Bool_t fMixedEvent;
    Int_t fMaxTrkDepth;
    Int_t fMaxEvtMixed;
    Int_t fMEnTrkPerEvt;
    Bool_t fFillMETree;
    //Tree
    Float_t ftPt_Tr;
    Float_t ftPhi_Tr;
    Float_t ftEta_Tr;
    Float_t ftCharge_Tr;
    Float_t ftMltOrCnt_Tr;
    Float_t ftZvtx_Tr;
    UInt_t ftPeriod_Tr;
    UInt_t ftOrbit_Tr;
    UShort_t ftBC_Tr;
    ClassDef(AliAnalysisTaskEbyEChargeSep, 1); // #of revison
};



//2. Below sub-class is just to Store track for ME
class AliMixedEvtTracks : public AliVParticle{
    
public:
    AliMixedEvtTracks(Float_t eta, Float_t phi, Float_t pt, Float_t p, Short_t charge):
    fEta(eta),
    fPhi(phi),
    fPt(pt),
    fP(p),
    fCharge(charge)
    {
     // default constructor
    }
    
    ~AliMixedEvtTracks(){}
    
    virtual Double_t Eta()   const {return fEta; }
    virtual Double_t Phi()   const {return fPhi; }
    virtual Double_t Pt()    const {return fPt; }
    virtual Double_t P()     const {return fP; }
    virtual Short_t Charge() const {return fCharge; }
    virtual Bool_t IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }

    //Unimplement proproties
    // kinematics
    virtual Double_t Px() const {AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const {AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const {AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Xv() const {AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const {AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const {AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
    virtual Double_t OneOverPt()  const {AliFatal("Not implemented"); return 0; }
    virtual Double_t Theta()      const {AliFatal("Not implemented"); return 0; }
    virtual Double_t E()          const {AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const {AliFatal("Not implemented"); return 0; }
    virtual Double_t Y()          const {AliFatal("Not implemented"); return 0; }
    virtual Int_t   GetLabel()    const {AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const {AliFatal("Not implemented"); return 0; }
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }

    
private:
    
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fPt;       // Pt
    Float_t fP;       // P
    Short_t fCharge;   // charge
    ClassDef(AliMixedEvtTracks, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};



//3. Sub-Class for Offline TREE but skipped (Not Used so far)
class AliMEOfflineTreeParam : public TObject{
public:
    
    AliMEOfflineTreeParam():
    phi_Tr(0),
    eta_Tr(0),
    pT_Tr(0),
    mult_Tr(0),
    charge_Tr(0),
    zVtx_Tr(0),
    period_Tr(0),
    orbit_Tr(0),
    BC_Tr(0)
    {
     //default constructor
    }
    
    ~AliMEOfflineTreeParam(){}

private:

    Float_t  phi_Tr;
    Float_t  eta_Tr;
    Float_t  pT_Tr;
    Short_t  charge_Tr;
    Float_t  mult_Tr;
    Float_t  zVtx_Tr;
    UInt_t   period_Tr;
    UInt_t   orbit_Tr;
    UShort_t BC_Tr;
    ClassDef(AliMEOfflineTreeParam, 1);
};

#endif


