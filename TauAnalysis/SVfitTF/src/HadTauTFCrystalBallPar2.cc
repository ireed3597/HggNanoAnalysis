#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBallPar2.h"

#include <iostream>

using namespace std;

HadTauTFCrystalBallPar2::HadTauTFCrystalBallPar2(int decayMode, int parNumber, const std::string& jetParFileName):
    parNumber_{parNumber},
    decayMode_{decayMode},
    jetParFileName_{jetParFileName},
    jetPar_{0},
    jetCorrector_{0}
{
  //std::cout << "<HadTauTFCrystalBallPar2::HadTauTFCrystalBallPar2>:" << std::endl;

  // Create the JetCorrectorParameter objects, the order does not matter.
  jetPar_ = new JetCorrectorParameters(jetParFileName_.data());

  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
  vPar_.push_back(*jetPar_);

  // Construct a FactorizedJetCorrector object with the vector
  jetCorrector_ = new FactorizedJetCorrector(vPar_);
}


HadTauTFCrystalBallPar2::HadTauTFCrystalBallPar2(const HadTauTFCrystalBallPar2& cbPar):
    parNumber_{cbPar.parNumber_},
    decayMode_{cbPar.decayMode_},
    jetParFileName_{cbPar.jetParFileName_},
    jetPar_{0},
    jetCorrector_{0}
{
  //std::cout << "<HadTauTFCrystalBallPar2::HadTauTFCrystalBallPar2>:" << std::endl;

  // Create the JetCorrectorParameter objects, the order does not matter.
  //std::cout << " initializing L2L3 corrections..." << std::endl;
  jetPar_ = new JetCorrectorParameters(jetParFileName_.data());
  //std::cout << "done." << std::endl;

  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
  vPar_.push_back(*jetPar_);

  // Construct a FactorizedJetCorrector object with the vector
  //std::cout << " initializing FactorizedJetCorrector..." << std::endl;
  jetCorrector_ = new FactorizedJetCorrector(vPar_);
  //std::cout << "done." << std::endl;
}

HadTauTFCrystalBallPar2::~HadTauTFCrystalBallPar2()
{
  //std::cout << "<HadTauTFCrystalBallPar2::~HadTauTFCrystalBallPar2>:" << std::endl;

  delete jetPar_;
  delete jetCorrector_;
}

double HadTauTFCrystalBallPar2::operator()(double genPt, double genEta) const
{
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite

  // Set the necessary values, for every jet, inside the jet loop
  jetCorrector_->setJetPt(genPt);
  jetCorrector_->setJetEta(genEta);

  // Get a vector of the individual correction factors
  double par = jetCorrector_->getCorrection();

  return par;
}
