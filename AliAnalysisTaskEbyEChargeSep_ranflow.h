
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
    void testmax(Int_t win, Float_t asyCut, Float_t *TempTotChrgPos,Float_t *TempTotChrgNeg);

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
    Int_t fEvtCentrality; //
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
    Float_t dbOld;
    Float_t dbpos;
    Float_t dbneg;
    Float_t asyold;
    Int_t   maxp;
    Int_t   maxn;
    Int_t   minp;
    Int_t   minn;
    Int_t maxpos; 
    Float_t db90Asy025;
    Int_t maxpos90;
    Float_t mdb90Asy025;
    Int_t mmaxpos90;
    Int_t rmaxpos90;
    Float_t rdb90Asy025;
    Float_t cosC;
    Float_t sinC;
    Float_t mcosc;
    Float_t msinc;
    Float_t   mcos1, mcos2, mcos3, mcos4, msin1, msin2, msin3, msin4;
    Float_t rcos1,rcos2,rcos3,rcos4,rsin1,rsin2,rsin3,rsin4;
    Float_t         mcos1p;
    Float_t         mcos2p,mcos3p;
    Float_t         mcos4p;
    Float_t         msin1p;
    Float_t         msin2p,msin3p;
    Float_t         msin4p;
    Float_t         mcos1n;
    Float_t         mcos2n,mcos3n;
    Float_t         mcos4n;
    Float_t         msin1n;
    Float_t         msin2n,msin3n;
    Float_t         msin4n;
    Float_t         cos1,cos2,cos3,cos4; //,cos5,cos6,cos7,cos8;
    Float_t         sin1,sin2,sin3,sin4; //,sin5,sin6,sin7,sin8;
    Float_t         cos1p,cos2p,cos3p,cos4p; //,cos5p,cos6p,cos7p,cos8p;
    Float_t         sin1p,sin2p,sin3p,sin4p; //,sin5p,sin6p,sin7p,sin8p;
    Float_t         cos1n,cos2n,cos3n,cos4n; //,cos5n,cos6n,cos7n,cos8n;
    Float_t         sin1n,sin2n,sin3n,sin4n; //,sin5n,sin6n,sin7n,sin8n;
    Float_t         rcos1p,rcos2p,rcos4p,rcos1n,rcos2n,rcos4n,rcos3p,rcos3n;
    Float_t         rsin1p,rsin2p,rsin4p,rsin1n,rsin2n,rsin4n,rsin3p,rsin3n;
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
    Int_t rphi;
    Float_t rpmax;
    float v2[8],gg[8];
    Float_t V0AEvtAngleRad_raw;
    Float_t V0AEvtAngleRad_plain;
    Float_t V0AEvtAngleRad_rec;
    Float_t V0AEvtAngleRad_align;
    Float_t V0AEvtAngleRad_twist;
    Float_t V0AEvtAngleRad_rescale;
    Float_t V0AEvtAngleRad_latest;
    Float_t V0CEvtAngleRad_raw;
    Float_t V0CEvtAngleRad_plain;
    Float_t V0CEvtAngleRad_rec;
    Float_t V0CEvtAngleRad_align;
    Float_t V0CEvtAngleRad_twist;
    Float_t V0CEvtAngleRad_rescale;
    Float_t V0CEvtAngleRad_latest;
    Float_t TPCEvtAngleRad_raw;
    Float_t TPCEvtAngleRad_plain;
    Float_t TPCEvtAngleRad_rec;
    Float_t TPCEvtAngleRad_align;
    Float_t TPCEvtAngleRad_twist;
    Float_t TPCEvtAngleRad_rescale;
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


