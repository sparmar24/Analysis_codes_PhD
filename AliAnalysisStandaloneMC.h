#ifndef AliAnalysisStandaloneMC_cxx
#define AliAnalysisStandaloneMC_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TH2F;
class TTree;
class AliESDEvent;
class AliESDPmdTrack;
class AliESDVertex;
class AliESDVZERO;
class AliESDFMD;

#include "AliAnalysisTaskSE.h"

Int_t out[8][10],outh[25],oute[20];
Float_t ma1, pmax1, cmax1;

class AliAnalysisStandaloneMC : public AliAnalysisTaskSE {
 public:
 AliAnalysisStandaloneMC() 
   : AliAnalysisTaskSE(), 
    fESD(0), 
    fOutputList(0),
    fTreeOut(),  
    fHistTotEvent(0),
    fCentrality(0),
    fCentrality1(0),
    fHistEta(0),
    fchTrk(0),
    fpmdTrk(0),
    chTrk(0),
    pmdTrk(0),
    fHistCent(0){
  }

  AliAnalysisStandaloneMC(const char *name);
  virtual ~AliAnalysisStandaloneMC() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  /*
    void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;
  }
  */
 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TTree       *fTreeOut;
  TH1F        *fHistTotEvent; //total event
  TH1F *fHistCent;
  TH1F *fHistEta;
  Int_t chTrk;
  Int_t pmdTrk;
  Int_t fchTrk;
  Int_t fpmdTrk;
  Float_t fCentrality1;
  Int_t fCentrality;

  TString fCentralityEstimator;
 
  AliAnalysisStandaloneMC(const AliAnalysisStandaloneMC&); // not implemented
  AliAnalysisStandaloneMC& operator=(const AliAnalysisStandaloneMC&); // not implemented
  
  ClassDef(AliAnalysisStandaloneMC, 1); // example of analysis
};

#endif
