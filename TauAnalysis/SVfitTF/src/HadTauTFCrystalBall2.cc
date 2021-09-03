#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"
#include "TauAnalysis/SVfitTF/interface/crystalBall.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include <iostream>
#include <string>
#include <assert.h>
#include <math.h>
#include "TMath.h"

namespace
{
  std::string findFile(const std::string& fileName)
  {
    edm::FileInPath inputFile(fileName);
    if ( inputFile.fullPath() == "" ) {
      std::cerr << "Error: Cannot find file = " << fileName << " !!" << std::endl;
      assert(0);
    }
    return inputFile.fullPath().data();
  }
}

HadTauTFCrystalBall2::HadTauTFCrystalBall2(const std::string& inputFilePath)
  : cbParSize_{19},
    nDecayMode_{7},
    decayMode_{kAll},
    xx_{0},
    pp_{0},
    genPt_cache_{-1.},
    genEta_cache_{0.}
{
  //std::cout << "<HadTauTFCrystalBall2::HadTauTFCrystalBall2>:" << std::endl;

  xx_ = new double[1];
  pp_ = new double[cbParSize_];

  std::map<int, std::vector<std::string>> parDir; // par number, directory

  std::vector<std::string> vdecay{"ak5tauHPSlooseCombDBcorrAll", // uncalibrated
            "ak5tauHPSlooseCombDBcorrOneProng0Pi0",
            "ak5tauHPSlooseCombDBcorrOneProng1Pi0",
            "ak5tauHPSlooseCombDBcorrTwoProng0Pi0",
            "ak5tauHPSlooseCombDBcorrTwoProng1Pi0",
            "ak5tauHPSlooseCombDBcorrThreeProng0Pi0",
            "ak5tauHPSlooseCombDBcorrThreeProng1Pi0"
          //"ak5tauHPSlooseCombDBcorrAlll2", // calibrated
            //"ak5tauHPSlooseCombDBcorrOneProng0Pi0l2",
            //"ak5tauHPSlooseCombDBcorrOneProng1Pi0l2",
            //"ak5tauHPSlooseCombDBcorrTwoProng0Pi0l2",
            //"ak5tauHPSlooseCombDBcorrTwoProng1Pi0l2",
            //"ak5tauHPSlooseCombDBcorrThreeProng0Pi0l2",
            //"ak5tauHPSlooseCombDBcorrThreeProng1Pi0l2",
          };

  for ( int i = 0; i < cbParSize_; i++ ){
    for ( auto &decay : vdecay ){
      parDir[i].push_back(inputFilePath + decay + "/" + std::to_string(i) + "/");
    }
  }

  std::map<int, std::string> fileCorr;
  fileCorr[kAll]            = "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrAll.txt";
  fileCorr[kOneProng0Pi0]   = "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrOneProng0Pi0.txt";
  fileCorr[kOneProng1Pi0]   = "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrOneProng1Pi0.txt";
  fileCorr[kTwoProng0Pi0]   = "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrTwoProng0Pi0.txt";
  fileCorr[kTwoProng1Pi0]   = "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrTwoProng1Pi0.txt";
  fileCorr[kThreeProng0Pi0] = "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrThreeProng0Pi0.txt";
  fileCorr[kThreeProng1Pi0] = "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrThreeProng1Pi0.txt";

  for ( int decayMode = 0; decayMode < nDecayMode_; ++decayMode ){
    for ( int par = 0; par < cbParSize_; ++par ){
      std::string jetParFileName = findFile(parDir[par][decayMode] + fileCorr[decayMode]); // each parameter has it own directory
      //std::cout << "reading tauES correction parameters for decayMode = " << decayMode << ", par #" << par << " from files:" << std::endl;
      //std::cout << " L2 = " << jetParFileName << std::endl;
      mapPar_[decayMode].push_back(new HadTauTFCrystalBallPar2(decayMode, par, jetParFileName));
    }
  }

  setDecayMode(decayMode_);
}

HadTauTFCrystalBall2::~HadTauTFCrystalBall2()
{
  //std::cout << "<HadTauTFCrystalBall2::~HadTauTFCrystalBall2>:" << std::endl;

  delete [] xx_;
  delete [] pp_;

  for ( int decayMode = 0; decayMode < nDecayMode_; ++decayMode ){
    for ( int par = 0; par < cbParSize_; ++par ){
      delete mapPar_[decayMode][par];
    }
  }
}

void HadTauTFCrystalBall2::setDecayMode(int decayMode) const
{
  //std::cout << "<HadTauTFCrystalBall2::setDecayMode>:" << std::endl;
  //std::cout << " decayMode = " << decayMode << std::endl;

  decayMode_ = decayMode;
  int idxDecayMode = -1;

  if      ( decayMode_ == reco::PFTau::kOneProng0PiZero   ) idxDecayMode = kOneProng0Pi0; // 0 <=> 1
  else if ( decayMode_ == reco::PFTau::kOneProng1PiZero   ) idxDecayMode = kOneProng1Pi0; // 1 <=> 2
  else if ( decayMode_ == reco::PFTau::kTwoProng0PiZero   ) idxDecayMode = kTwoProng0Pi0; // 5 <=> 3
  else if ( decayMode_ == reco::PFTau::kTwoProng1PiZero   ) idxDecayMode = kTwoProng1Pi0; // 6 <=> 4
  else if ( decayMode_ == reco::PFTau::kThreeProng0PiZero ) idxDecayMode = kThreeProng0Pi0; // 10 <=> 5
  else if ( decayMode_ == reco::PFTau::kThreeProng1PiZero ) idxDecayMode = kThreeProng1Pi0; // 11 <=> 6
  else {
    //std::cerr << "Warning: No transfer function defined for decay mode = " << decayMode_ << " !!";
    idxDecayMode = kAll; // 0
  }
  assert(mapPar_.find(idxDecayMode) != mapPar_.end());
  //std::cout << "set decay mode to: "<<idxDecayMode<<std::endl;

  thePar_.resize(cbParSize_);
  for ( int par = 0; par < cbParSize_; ++par ) {
    thePar_[par] = mapPar_.find(idxDecayMode)->second[par];
  }
  //std::cout << "setting decay mode done!! "<<std::endl;
}

double HadTauTFCrystalBall2::operator()(double recPt, double genPt, double genEta) const
{
  xx_[0] = recPt/genPt;

  if ( genPt != genPt_cache_ || genEta != genEta_cache_ ) {
    for ( int iPar = 0; iPar < cbParSize_; ++iPar ) {
      pp_[iPar] = (*thePar_[iPar])(genPt, genEta);
    }
    genPt_cache_ = genPt;
    genEta_cache_ = genEta;
  }

  double retVal = crystalBall(xx_, pp_);
  double normVal = normalizedCrystalBall(pp_);
  retVal /= normVal;

  if ( isinf(retVal) || isnan(retVal) ) {
    std::cout << "<HadTauTFCrystalBall2::operator()>:" << std::endl;
    std::cout << " pT: rec = " << recPt << ", gen = " << genPt << std::endl;
    std::cout << " eta = " << genEta << std::endl;
    std::cout << "xx = " << xx_[0] << std::endl;
    for ( int iPar = 0; iPar < cbParSize_; ++iPar ) {
      pp_[iPar] = (*thePar_[iPar])(genPt, genEta);
      std::cout << "pp(" << iPar << ") = " << pp_[iPar] << std::endl;
    }
  }

  return retVal;
}

double HadTauTFCrystalBall2::integral(double recPt_low, double recPt_up, double genPt, double genEta) const
{
  if ( !(recPt_low < recPt_up) ) return 0.;

  double x_low = recPt_low/genPt;
  double x_up = recPt_up/genPt;

  if ( genPt != genPt_cache_ || genEta != genEta_cache_ ) {
    for ( int iPar = 0; iPar < cbParSize_; ++iPar ) {
      pp_[iPar] = (*thePar_[iPar])(genPt, genEta);
    }
    genPt_cache_ = genPt;
    genEta_cache_ = genEta;
  }

  double xx[1];
  xx[0] = (x_up + x_low)/(2.*genPt);
  double integral = crystalBall(xx, pp_);

  if ( isinf(integral) || isnan(integral) ) {
    std::cout << "<HadTauTFCrystalBall2::integral()>:" << std::endl;
    std::cout << " pT: rec = " << recPt_low << ".." << recPt_up << ", gen = " << genPt << std::endl;
    std::cout << " eta = " << genEta << std::endl;
    std::cout << "x = " << x_low << ".." << x_up << std::endl;
    for ( int iPar = 0; iPar < cbParSize_; ++iPar ) {
      pp_[iPar] = (*thePar_[iPar])(genPt, genEta);
      std::cout << "pp(" << iPar << ") = " << pp_[iPar] << std::endl;
    }
  }

  return integral;
}

HadTauTFCrystalBall2* HadTauTFCrystalBall2::Clone(const std::string& label) const
{
  HadTauTFCrystalBall2* clone = new HadTauTFCrystalBall2();
  clone->cbParSize_ = cbParSize_;
  clone->nDecayMode_ = nDecayMode_;
  clone->decayMode_ = decayMode_;
  clone->xx_ = new double[1];
  clone->xx_[0] = xx_[0];
  clone->pp_ = new double[cbParSize_];
  for ( int iPar = 1; iPar < cbParSize_; ++iPar ) {
    clone->pp_[iPar] = pp_[iPar];
  }
  for ( std::map<int, vHadTauTFCrystalBallParPtr>::const_iterator entryPar = mapPar_.begin();
  entryPar != mapPar_.end(); ++entryPar ) {
    const vHadTauTFCrystalBallParPtr& cbPars = entryPar->second;
    vHadTauTFCrystalBallParPtr cbPars_cloned;
    for ( int iPar = 0; iPar < cbParSize_; ++iPar ) {
      const HadTauTFCrystalBallPar2* cbPar = cbPars[iPar];
      cbPars_cloned.push_back(new HadTauTFCrystalBallPar2(*cbPar));
    }
    clone->mapPar_[decayMode_] = cbPars_cloned;
  }
  clone->setDecayMode(decayMode_);
  clone->genPt_cache_ = genPt_cache_;
  clone->genEta_cache_ = genEta_cache_;
  return clone;
}
