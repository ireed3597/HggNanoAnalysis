#ifndef TauAnalysis_SVfitTF_HadTauTFCrystalBall2_h
#define TauAnalysis_SVfitTF_HadTauTFCrystalBall2_h

// This class returns the CrystalBall probability
// value for parameters produced with the JetMET packadge.
// Author: Betty Calpas, Christian Veelken
// Email: betty.calpas@cern.ch, christian.veelken@cern.ch

#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBallPar2.h"
#include <string>
#include <map>
#include <vector>

class HadTauTFCrystalBall2 : public HadTauTFBase
{
 public:

  HadTauTFCrystalBall2(const std::string& = "TauAnalysis/SVfitTF/data/L2L3Corr/");

  virtual ~HadTauTFCrystalBall2();

  double operator()(double recPt, double genPt, double genEta) const; // return CB probability

  virtual double integral(double recPt_min, double recPt_max, double genPt, double genEta) const;

  void setDecayMode(int decayMode) const;

  virtual HadTauTFCrystalBall2* Clone(const std::string& label) const;

  const HadTauTFCrystalBallPar2* getPar(int par);

 private:

  bool calibrated_; // calibrated or uncalibrated

  int cbParSize_;

  int nDecayMode_;

  enum { kAll, kOneProng0Pi0, kOneProng1Pi0, kTwoProng0Pi0, kTwoProng1Pi0, kThreeProng0Pi0, kThreeProng1Pi0 };

  mutable int decayMode_;

  mutable double *xx_; // response

  mutable double *pp_; // response correction parameters array

  typedef std::vector<HadTauTFCrystalBallPar2*> vHadTauTFCrystalBallParPtr;

  mutable vHadTauTFCrystalBallParPtr thePar_;

  std::map<int, vHadTauTFCrystalBallParPtr> mapPar_; // <decayMode, parametersVec>

  mutable double genPt_cache_;

  mutable double genEta_cache_;

};

#endif
