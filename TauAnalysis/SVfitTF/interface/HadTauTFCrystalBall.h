#ifndef TauAnalysis_SVfitTF_HadTauTFCrystalBall_h
#define TauAnalysis_SVfitTF_HadTauTFCrystalBall_h

// This class returns the CrystalBall probability
// value for parameters produced with the JetMET packadge.
// Author: Betty Calpas, Christian Veelken
// Email: betty.calpas@cern.ch, christian.veelken@cern.ch

#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBallPar.h"

#include <string>
#include <map>
#include <vector>

class HadTauTFCrystalBall : public HadTauTFBase
{
 public:
  HadTauTFCrystalBall(const std::string& = "TauAnalysis/SVfitTF/data/L2L3Corr/"); // default constructor
  virtual ~HadTauTFCrystalBall(); // default destructor
  double operator()(double recPt, double genPt, double genEta) const; // call the function "fnc_dscb" which return the CrystalBall probalility value

  double integral(double recPt_low, double recPt_up, double genPt, double genEta) const;

  void setDecayMode(int decayMode) const;

  const HadTauTFCrystalBallPar* getPar(int par);

  virtual HadTauTFCrystalBall* Clone(const std::string& label) const;

 private:
  int cbParSize_;  // nb of CrystalBall parameter
  int nDecayMode_; // nb of different supported tau decay modes
  enum { kAll, kOneProng0Pi0, kOneProng1Pi0, kOneProng2Pi0, kThreeProng0Pi0 };
  mutable int decayMode_;  // reconstructed decay mode of given tau (to be set for each event by calling "setDecayMode" function)
  mutable double* xx_;
  mutable double* pp_; // point to the 7 parameter function value
  typedef std::vector<HadTauTFCrystalBallPar*> vHadTauTFCrystalBallParPtr;
  mutable vHadTauTFCrystalBallParPtr thePar_;
  std::map<int, vHadTauTFCrystalBallParPtr> mapPar_; // <decayMode, parametersVec>
  mutable double genPt_cache_;
  mutable double genEta_cache_;
};

#endif
