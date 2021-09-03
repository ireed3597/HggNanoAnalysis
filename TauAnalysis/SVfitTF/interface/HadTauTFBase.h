#ifndef TauAnalysis_SVfitTF_HadTauTFBase_h
#define TauAnalysis_SVfitTF_HadTauTFBase_h

/** \class HadTauTFBase
 *
 * Abstract base class defining interface for hadronic tau pT transfer functions.
 *
 * \author Christian Veelken, Tallinn
 *
 */

#include <string>

class HadTauTFBase
{
 public:
  HadTauTFBase() {}
  virtual ~HadTauTFBase() {}
  virtual double operator()(double recPt, double genPt, double genEta) const = 0;

  virtual double integral(double recPt_min, double recPt_max, double genPt, double genEta) const = 0;

  virtual void setDecayMode(int decayMode) const {}

  virtual HadTauTFBase* Clone(const std::string& label) const = 0;
};

#endif
