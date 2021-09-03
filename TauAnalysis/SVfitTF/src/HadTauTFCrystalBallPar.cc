#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBallPar.h"

#include <iostream>

using namespace std;

HadTauTFCrystalBallPar::HadTauTFCrystalBallPar(int parNumber, int decayMode, const std::string& jetParFileName)
  : parNumber_(parNumber),
    decayMode_(decayMode),
    jetParFileName_(jetParFileName),
    jetPar_(0),
    l2JetPar_(0),
    l3JetPar_(0),
    jetCorrector_(0)
{
  //std::cout << "<HadTauTFCrystalBallPar::HadTauTFCrystalBallPar>:" << std::endl;

  // Create the JetCorrectorParameter objects, the order does not matter.
  jetPar_ = new JetCorrectorParameters(jetParFileName_.data());

  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
  vPar_.push_back(*jetPar_);

  // Construct a FactorizedJetCorrector object with the vector
  jetCorrector_ = new FactorizedJetCorrector(vPar_);
}

HadTauTFCrystalBallPar::HadTauTFCrystalBallPar(int parNumber, int decayMode, const std::string& l2JetParFileName, const std::string& l3JetParFileName)
  : parNumber_(parNumber),
    decayMode_(decayMode),
    l2JetParFileName_(l2JetParFileName),
    l3JetParFileName_(l3JetParFileName),
    jetPar_(0),
    l2JetPar_(0),
    l3JetPar_(0),
    jetCorrector_(0)
{
  //std::cout << "<HadTauTFCrystalBallPar::HadTauTFCrystalBallPar>:" << std::endl;

  // Create the JetCorrectorParameter objects, the order does not matter.
  //std::cout << " initializing L2 corrections..." << std::endl;
  l2JetPar_ = new JetCorrectorParameters(l2JetParFileName_.data());
  //std::cout << " initializing L3 corrections..." << std::endl;
  l3JetPar_ = new JetCorrectorParameters(l3JetParFileName_.data());
  //std::cout << "done." << std::endl;

  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
  vPar_.push_back(*l2JetPar_);
  vPar_.push_back(*l3JetPar_);

  // Construct a FactorizedJetCorrector object with the vector
  //std::cout << " initializing FactorizedJetCorrector..." << std::endl;
  jetCorrector_ = new FactorizedJetCorrector(vPar_);
  //std::cout << "done." << std::endl;
}

HadTauTFCrystalBallPar::HadTauTFCrystalBallPar(const HadTauTFCrystalBallPar& cbPar)
  : parNumber_(cbPar.parNumber_),
    decayMode_(cbPar.decayMode_),
    l2JetParFileName_(cbPar.l2JetParFileName_),
    l3JetParFileName_(cbPar.l3JetParFileName_),
    jetPar_(0),
    l2JetPar_(0),
    l3JetPar_(0),
    jetCorrector_(0)
{
  //std::cout << "<HadTauTFCrystalBallPar::HadTauTFCrystalBallPar>:" << std::endl;

  // Create the JetCorrectorParameter objects, the order does not matter.
  //std::cout << " initializing L2 corrections..." << std::endl;
  l2JetPar_ = new JetCorrectorParameters(l2JetParFileName_.data());
  //std::cout << " initializing L3 corrections..." << std::endl;
  l3JetPar_ = new JetCorrectorParameters(l3JetParFileName_.data());
  //std::cout << "done." << std::endl;

  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
  vPar_.push_back(*l2JetPar_);
  vPar_.push_back(*l3JetPar_);

  // Construct a FactorizedJetCorrector object with the vector
  //std::cout << " initializing FactorizedJetCorrector..." << std::endl;
  jetCorrector_ = new FactorizedJetCorrector(vPar_);
  //std::cout << "done." << std::endl;
}

HadTauTFCrystalBallPar::~HadTauTFCrystalBallPar()
{
  //std::cout << "<HadTauTFCrystalBallPar::~HadTauTFCrystalBallPar>:" << std::endl;

  delete jetPar_;
  delete l2JetPar_;
  delete l3JetPar_;
  delete jetCorrector_;
}

double HadTauTFCrystalBallPar::operator()(double genPt, double genEta) const
{
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite

  // Set the necessary values, for every jet, inside the jet loop
  jetCorrector_->setJetPt(genPt);
  jetCorrector_->setJetEta(genEta);

  // Get a vector of the individual correction factors
  double correction = jetCorrector_->getCorrection();

  return correction;
}
