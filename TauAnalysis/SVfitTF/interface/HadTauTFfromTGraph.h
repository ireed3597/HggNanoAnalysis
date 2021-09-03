#ifndef TauAnalysis_SVfitTF_HadTauTFfromTGraph_h
#define TauAnalysis_SVfitTF_HadTauTFfromTGraph_h

/** \class HadTauTFfromTGraph
 *
 * Read transfer functions for hadronic tau pT as TGraph objects from ROOT file.
 *
 * \author Christian Veelken, Tallinn
 *
 */

#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"

#include <TGraph.h>

#include <map>
#include <assert.h>

class HadTauTFfromTGraph : public HadTauTFBase
{
 public:
  HadTauTFfromTGraph();
  HadTauTFfromTGraph(const TGraph* resolution);
  HadTauTFfromTGraph(const std::map<int, const TGraph*>& resolutionMap);
  ~HadTauTFfromTGraph();
  double operator()(double recPt, double genPt, double genEta) const;

  double integral(double recPt_low, double recPt_up, double genPt, double genEta) const;

  void setDecayMode(int decayMode) const;

  HadTauTFfromTGraph* Clone(const std::string& label) const;

 private:
  struct ResolutionMapEntry
  {
    ResolutionMapEntry(const TGraph* resolution)
      : resolution_(0),
        cdf_(0)
    {
      assert(resolution);
      resolution_ = (TGraph*)resolution->Clone(Form("%s_cloned", resolution->GetName()));
      int numPoints = resolution_->GetN();
      cdf_ = new TGraph(numPoints);
      double integral = 0.;
      for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
        double x, y;
        resolution_->GetPoint(idxPoint, x, y);
        integral += y;
        cdf_->SetPoint(idxPoint, x, integral);
      }
    }
    ~ResolutionMapEntry()
    {
      delete resolution_;
      delete cdf_;
    }
    const TGraph* resolution_;
    TGraph* cdf_;
  };

  std::map<int, const ResolutionMapEntry*> resolutionMap_; // key = decayMode
  bool hasResolutionMap_;
  mutable int theDecayMode_;
  mutable const ResolutionMapEntry* theResolution_;
  mutable double xMin_; // range of validity of TGraph object (x-coordinate of first and last point)
  mutable double xMax_;
};

#endif
