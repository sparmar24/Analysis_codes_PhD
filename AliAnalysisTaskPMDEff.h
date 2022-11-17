#ifndef AliAnalysisTaskPMDEff_cxx
#define AliAnalysisTaskPMDEff_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDPmdTrack;
class AliESDVertex;
class AliESDVZERO;
class AliESDFMD;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPMDEff : public AliAnalysisTaskSE {
 public:
 AliAnalysisTaskPMDEff() 
   : AliAnalysisTaskSE(), 
    fESD(0), 
    fOutputList(0),  
    fHistTotEvent(0),
    fHistXYPre(0), 
    fHistXYCpv(0),
    fHistCent(0),
    fHistGammaLikeMipCut1Ncell2(0),
    fHistGammaTrueMipCut1Ncell2(0),
    fHistGammaLikeMipCut2Ncell2(0),
    fHistGammaTrueMipCut2Ncell2(0),
    fHistGammaLikeMipCut1Ncell0(0),
    fHistGammaTrueMipCut1Ncell0(0),
    fHistGammaLikeMipCut2Ncell0(0),
    fHistGammaTrueMipCut2Ncell0(0),
    fHistMultPhotonInc(0),
    fCentralityEstimator("V0M"){
    for(Int_t i=0; i<8; i++){
      fHistGammaLikeEtaBinMipCut1Ncell2[i] = 0;
      fHistGammaTrueEtaBinMipCut1Ncell2[i] = 0;
      fHistGammaLikeEtaBinMipCut2Ncell2[i] = 0;
      fHistGammaTrueEtaBinMipCut2Ncell2[i] = 0;
      fHistGammaLikeEtaBinMipCut1Ncell0[i] = 0;
      fHistGammaTrueEtaBinMipCut1Ncell0[i] = 0;
      fHistGammaLikeEtaBinMipCut2Ncell0[i] = 0;
      fHistGammaTrueEtaBinMipCut2Ncell0[i] = 0;
      fHistMultPhotonIncEtaBin[i] = 0;
      fHistXYPreEtaBin[i]=0; 
      fHistXYCpvEtaBin[i]=0;
    }
  }

  AliAnalysisTaskPMDEff(const char *name);
  virtual ~AliAnalysisTaskPMDEff() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

 
    void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;
  }
 
 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TH1F        *fHistTotEvent; //total event
  TH2F *fHistXYPre;
  TH2F *fHistXYCpv;
  TH1F *fHistCent;
  TH1F *fHistGammaLikeEtaBinMipCut1Ncell2[8];
  TH1F *fHistGammaTrueEtaBinMipCut1Ncell2[8];
  TH1F *fHistGammaLikeEtaBinMipCut2Ncell2[8];
  TH1F *fHistGammaTrueEtaBinMipCut2Ncell2[8];
  TH1F *fHistGammaLikeEtaBinMipCut1Ncell0[8];
  TH1F *fHistGammaTrueEtaBinMipCut1Ncell0[8];
  TH1F *fHistGammaLikeEtaBinMipCut2Ncell0[8];
  TH1F *fHistGammaTrueEtaBinMipCut2Ncell0[8];
  TH1F *fHistMultPhotonIncEtaBin[8];
  TH2F *fHistXYPreEtaBin[8]; 
  TH2F *fHistXYCpvEtaBin[8];
  TH1F *fHistGammaLikeMipCut1Ncell2;
  TH1F *fHistGammaTrueMipCut1Ncell2;
  TH1F *fHistGammaLikeMipCut2Ncell2;
  TH1F *fHistGammaTrueMipCut2Ncell2;
  TH1F *fHistGammaLikeMipCut1Ncell0;
  TH1F *fHistGammaTrueMipCut1Ncell0;
  TH1F *fHistGammaLikeMipCut2Ncell0;
  TH1F *fHistGammaTrueMipCut2Ncell0;
  //__________MC_______
  TH1F *fHistMultPhotonInc;

  TString fCentralityEstimator;
 
  AliAnalysisTaskPMDEff(const AliAnalysisTaskPMDEff&); // not implemented
  AliAnalysisTaskPMDEff& operator=(const AliAnalysisTaskPMDEff&); // not implemented
  
  ClassDef(AliAnalysisTaskPMDEff, 1); // example of analysis
};

#endif
