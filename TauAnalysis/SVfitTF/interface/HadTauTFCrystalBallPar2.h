#ifndef TauAnalysis_SVfitTF_HadTauTFCrystalBallPar2_h
#define TauAnalysis_SVfitTF_HadTauTFCrystalBallPar2_h

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <string>
#include <vector>

class HadTauTFCrystalBallPar2
{
  public:
    HadTauTFCrystalBallPar2(int decayMode, int parNumber, const std::string& jetParFileName);
    HadTauTFCrystalBallPar2(const HadTauTFCrystalBallPar2& cbPar);
    virtual ~HadTauTFCrystalBallPar2();
    double operator()(double genPt, double genEta) const; // return CB parameters vs (pt, eta)

  private:
    int parNumber_;
    int decayMode_;
    std::string jetParFileName_;
    JetCorrectorParameters* jetPar_;
    FactorizedJetCorrector* jetCorrector_;
    std::vector<JetCorrectorParameters> vPar_;
};

#endif

