#include "TauAnalysis/SVfitTF/interface/HadTauTFfromTGraph.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <iostream>
#include <string>
#include <assert.h>
#include <math.h>

namespace
{
  void initializeRange(const TGraph* graph, double& xMin, double& xMax)
  {
    int numPoints = graph->GetN();
    assert(numPoints >= 1);
    for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
      double x, y;
      graph->GetPoint(idxPoint, x, y);
      if ( idxPoint == 0 || x < xMin ) xMin = x;
      if ( idxPoint == 0 || x > xMax ) xMax = x;
    }
  }
}

HadTauTFfromTGraph::HadTauTFfromTGraph()
  : hasResolutionMap_(false),
    theDecayMode_(-1),
    theResolution_(0),
    xMin_(-1.),
    xMax_(-1.)
{}

HadTauTFfromTGraph::HadTauTFfromTGraph(const TGraph* resolution)
  : hasResolutionMap_(false),
    theDecayMode_(-1),
    theResolution_(new ResolutionMapEntry(resolution))
{
  resolutionMap_[-1] = theResolution_;
  initializeRange(theResolution_->resolution_, xMin_, xMax_);
}

HadTauTFfromTGraph::HadTauTFfromTGraph(const std::map<int, const TGraph*>& resolutionMap)
  : hasResolutionMap_(true),
    theDecayMode_(-1),
    theResolution_(0)
{
  for ( std::map<int, const TGraph*>::const_iterator resolutionEntry = resolutionMap.begin();
  resolutionEntry != resolutionMap.end(); ++resolutionEntry ) {
    int decayMode = resolutionEntry->first;
    const TGraph* resolution = resolutionEntry->second;
    resolutionMap_[decayMode] = new ResolutionMapEntry(resolution);
  }
}

HadTauTFfromTGraph::~HadTauTFfromTGraph()
{
  if ( hasResolutionMap_ ) {
    for ( std::map<int, const ResolutionMapEntry*>::iterator it = resolutionMap_.begin();
          it != resolutionMap_.end(); ++it ) {
      delete it->second;
    }
  } else {
    delete theResolution_;
  }
}

void HadTauTFfromTGraph::setDecayMode(int decayMode) const
{
  if ( hasResolutionMap_ ) {
    if ( resolutionMap_.find(decayMode) == resolutionMap_.end() ) {
      std::cerr << "No tau pT transfer functions defined for decayMode = " << decayMode << " !!" << std::endl;
      assert(0);
    }
    theResolution_ = resolutionMap_.find(decayMode)->second;
    initializeRange(theResolution_->resolution_, xMin_, xMax_);
  }
}

double HadTauTFfromTGraph::operator()(double recPt, double genPt, double genEta) const
{
  if ( !theResolution_ ) {
    std::cerr << "No tau pT transfer functions defined, call 'setDecayMode' function first !!" << std::endl;
    assert(0);
  }
  if ( genPt > 0. ) {
    double x = recPt/genPt;
    if ( x < 0. ) return 0.;
    if ( x < xMin_ ) x = xMin_;
    if ( x > xMax_ ) x = xMax_;
    return theResolution_->resolution_->Eval(x);
  } else {
    return 0.;
  }
}

double HadTauTFfromTGraph::integral(double recPt_low, double recPt_up, double genPt, double genEta) const
{
  if ( !(recPt_low < recPt_up) ) return 0.;

  if ( !theResolution_ ) {
    std::cerr << "No tau pT transfer functions defined, call 'setDecayMode' function first !!" << std::endl;
    assert(0);
  }
  if ( genPt > 0. ) {
    double x_low = recPt_low/genPt;
    if ( x_low < xMin_ ) x_low = xMin_;
    if ( x_low > xMax_ ) x_low = xMax_;
    double x_up = recPt_up/genPt;
    if ( x_up < xMin_ ) x_up = xMin_;
    if ( x_up > xMax_ ) x_up = xMax_;
    return (theResolution_->cdf_->Eval(x_up) - theResolution_->cdf_->Eval(x_low));
  } else {
    return 0.;
  }
}

HadTauTFfromTGraph* HadTauTFfromTGraph::Clone(const std::string& label) const
{
  HadTauTFfromTGraph* clone = new HadTauTFfromTGraph();
  for ( std::map<int, const ResolutionMapEntry*>::const_iterator resolutionEntry = resolutionMap_.begin();
        resolutionEntry != resolutionMap_.end(); ++resolutionEntry ) {
    int decayMode = resolutionEntry->first;
    const TGraph* resolution = resolutionEntry->second->resolution_;
    clone->resolutionMap_[decayMode] = new ResolutionMapEntry(resolution);
  }
  clone->hasResolutionMap_ = hasResolutionMap_;
  clone->theDecayMode_ = theDecayMode_;
  if ( clone->resolutionMap_.find(theDecayMode_) != clone->resolutionMap_.end() ) {
    clone->theResolution_ = clone->resolutionMap_.find(theDecayMode_)->second;
    initializeRange(clone->theResolution_->resolution_, clone->xMin_, clone->xMax_);
  } else {
    clone->theResolution_ = 0;
  }
  return clone;
}
