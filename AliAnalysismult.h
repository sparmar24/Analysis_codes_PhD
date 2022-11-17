#ifndef AliAnalysismult_cxx
#define AliAnalysismult_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TH2F;
class TH1D;
class TH2D;
class TTree;
class NTuple;
class TProfile;
class TAxis;
class AliPhysicsSelection;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"
const Int_t nCentralityBins =5;
Int_t out[8][10],outh[25],oute[20],outf[3],outp[10][10];
/*
// Float_t pf[14],mf[20];
Float_t pf[0]=0.896530628,pf[1]=1.14103901,pf[2]=0.996145129,pf[3]=1.45946848;
Float_t pf[4]=1.04595244,pf[5]=0.871626973,pf[6]=0.922899187,pf[7]=0.896530628;
Float_t pf[8]=1.04595244,pf[9]=0.883903444,pf[10]=0.896530628,pf[11]=1.26526499;
Float_t pf[12]=1.02880561,pf[13]=0.950865805;
Float_t mf[0]=0.997369885,mf[1]=0.96375072,mf[2]=0.960154593,mf[3]=0.942569375;
Float_t mf[4]=0.949525595,mf[5]=0.951280713,mf[6]=0.978408515,mf[7]=0.993519068;
Float_t mf[8]=1.03341937,mf[9]=1.04178715,mf[10]=1.07441103,mf[11]=1.07665873;
Float_t mf[12]=1.05243945,mf[13]=1.01307654;


const Float_t fpmd[20];ffmd[20];
fpmd(0)=0.896530628;fpmd(1)=1.14103901;fpmd(2)=0.996145129;fpmd(3)=1.45946848;
fpmd(4)=1.04595244;fpmd(5)=0.871626973;fpmd(6)=0.922899187;fpmd(7)=0.896530628;
fpmd(8)=1.04595244;fpmd(9)=0.883903444;fpmd(10)=0.896530628;fpmd(11)=1.26526499;
fpmd(12)=1.02880561;fpmd(13)=0.950865805;
ffmd(0)=0.997369885;ffmd(1)=0.96375072;ffmd(2)=0.960154593;ffmd(3)=0.942569375;
ffmd(4)=0.949525595;ffmd(5)=0.951280713;ffmd(6)=0.978408515;ffmd(7)=0.993519068;
ffmd(8)=1.03341937;ffmd(9)=1.04178715;ffmd(10)=1.07441103;ffmd(11)=1.07665873;
ffmd(12)=1.01307654;ffmd(13)=1.05243945;
*/
class AliAnalysismult : public AliAnalysisTaskSE {
 public:
 AliAnalysismult() : AliAnalysisTaskSE(), fESD(0), fOutputList(0),fTreeOut() {}
  AliAnalysismult(const char *name);
  virtual ~AliAnalysismult() {}
  void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;}
  
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;    //! ESD object
  TTree   *fTreeOut;
  TList       *fOutputList; //! Output list
  
  TH1F        *fHistPhiBin[nCentralityBins];
  TH1F        *fHistPhiFull[nCentralityBins];
  TH1F        *fHistPhiFMD[nCentralityBins];
  TH1F        *fHistPhiFMD253[nCentralityBins];
  TH1F        *fHistFMDmult[nCentralityBins];
  TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD" 
  AliAnalysismult(const AliAnalysismult&); // not implemented
  AliAnalysismult& operator=(const AliAnalysismult&); // not implemented
  
  ClassDef(AliAnalysismult, 1); // example of analysis
};

#endif
