
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
    
    void SetCentralityEstimator(const char* CentEst){fCentEstimator = CentEst;}
    void SetAODProdcNumber(const char* AODnum){fAODProdNumber = AODnum;}
    void SetFilterBitNumber(Int_t Bit){fBit = Bit;}
    void SetSwitchToDBCalc(Bool_t Sel, Int_t wnd){IsDBMethodOn = Sel; fDBWindow = wnd;}
    
    void SetEvtZvtxValue(Double_t Sel){fEvtabsZvtx = Sel;}
    void SetTrkPtRange(Double_t Sel1, Double_t Sel2){fTrkPtMin = Sel1; fTrkPtMax = Sel2;}
    void SetTrkEtaRange(Double_t Sel1, Double_t Sel2){fTrkEtaMin = Sel1; fTrkEtaMax = Sel2;}
    void SetTrkabsDCAXY(Double_t Sel){fTrkabsDCAxy = Sel;}
    void SetTrkabsDCAZ(Double_t Sel){fTrkabsDCAz = Sel;}
    void SetTrkMinClusterTPC(Int_t Sel){fTrkTPCclus = Sel;}
    void SetMCAnalysis(Bool_t Sel1, Bool_t Sel2){fReadMC = Sel1, fMCKine=Sel2;}
    void SetMEAnalysis(Bool_t Sel1){fMixedEvent = Sel1;}

    void SetEfficiencyCorrection(Bool_t Sel1){fIsEffCorr = Sel1;}
    void SetInputEfficiencyFile(TH1D *hMap){if(fhEffValues)delete fhEffValues;fhEffValues=(TH1D*)hMap->Clone();}

    void SetExpectedCorrectionPass(const char *pass) { fExpectedCorrectionPass = pass; }
    void SetAlternativeCorrectionPass(const char *pass) { fAlternativeCorrectionPass = pass; }

    /*    Float_t Cal3PartCor(float x,float y,float x1,float y1, Int_t num1);
    Float_t Cal2PartCor(float x2,float y2, Float_t num2);
    Float_t Cal3PartCorUn(float xp,float yp,float x2p,float y2p, float xn,float yn,float x2n,float y2n, Int_t num1, Int_t num2);
    Float_t Cal4PartCor(float cos2x,float sin2y,float cos4x,float sin4y, Int_t num1);
    void CalEvtPhiDumbellMethod(Int_t fPhiWindow, Float_t *TempTotChargePos,Float_t *TempTotChargeNeg);*/
    //    void testmax(Int_t win, Float_t asyCut, Float_t *TempTotChrgPos,Float_t *TempTotChrgNeg);
    //void CalCorrelatorVia2PCLikeSign(Int_t nTrk, Float_t *TrkArray, Float_t V0APlane,Float_t V0CPlane,Float_t TPCPlane);

    // Getters                                                                                   
    Double_t GetTrackWeight(Int_t PtBin);

    //Booliens                                                                                    
    Bool_t IsTrackEffMapLoaded(){if(fhEffValues) return kTRUE; else return kFALSE;}
            
private:
    
    AliAnalysisTaskEbyEChargeSep(const AliAnalysisTaskEbyEChargeSep&);
    AliAnalysisTaskEbyEChargeSep& operator=(const AliAnalysisTaskEbyEChargeSep&);
    TObjArray* CloneAcceptAndReduceTracks(TObjArray* fMEtracks);

    
    AliAODEvent *fAOD;  //! AOD Event object
    TList *fOutputList; //! Output list

    AliQnCorrectionsManager *fFlowQnVectorMgr;
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask;
    TString fExpectedCorrectionPass;
    TString fAlternativeCorrectionPass;


    TClonesArray *farrayMC;//!
    Float_t fEvtCentrality; //
    Float_t fEvtCentOrMult; //
    Int_t fRunNumber; //
    TString fCentEstimator;
    Int_t fEventCounter; //
    TH1I* fhTotEvent;//!
    TH1I* fhTotTrack;//!
    TH1I* fhTotPart;//!
    TClonesArray* fmcArray;//!
    TTree *fDBEvtTree; //!  
    Bool_t IsEvtSctPlotReq;
    Float_t fEvtCentCL1;
    Float_t fEvtCentTRK;
    Float_t fV0AMult; 
    Float_t fV0CMult;
    Float_t fV0ACMult;
    Int_t fBunchCrossNo;                                                           
    Int_t fOrbitNo;                                   
    Int_t fPeriodNo;           
    TString fAODProdNumber;
    Int_t fCAODNumber;
    Int_t fRunNumPlotReq;
    Int_t fFileNumPlotReq;
    Int_t fEvtNumPlotReq;
    Int_t fCompareFileNo;
    Int_t fTreeEventCounter;
    Int_t fBit;
    Bool_t IsDBMethodOn;
    Bool_t fReadMC;
    Bool_t fMCKine;
    Int_t fDBWindow;
    Double_t fEvtabsZvtx;
    Double_t fEventZvtx;
    Double_t fTrkPtMin;
    Double_t fTrkPtMax;
    Double_t fTrkEtaMin;
    Double_t fTrkEtaMax;
    Double_t fTrkabsDCAxy;
    Double_t fTrkabsDCAz;
    UInt_t fTrkTPCclus;
    Int_t flMultiAllChrg;
    Int_t flMultPosChrg;
    Int_t flMultNegChrg;
    Int_t fPhiRecoDgr;
    Int_t fPhiKineDgr;
    TH1D  *fhEffValues;//
    Bool_t fIsEffCorr;
    Bool_t fMixedEvent;
    Float_t fVtrxZ; // 
    Float_t V0AEvtAngleRad_latest;
    Float_t V0CEvtAngleRad_latest;
    Float_t TPCEvtAngleRad_latest;


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


#endif


