#ifndef FastMTT_FastMTT_H
#define FastMTT_FastMTT_H

#include <string>
#include <vector>
#include <tuple>
#include <bitset>
#include "TMatrixD.h"
#include "TBenchmark.h"
#include "Math/LorentzVector.h"

namespace classic_svFit{
  class MeasuredTauLepton;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

}

namespace ROOT{
  namespace Math{
    class Minimizer;
    class Functor;
  }
}
class TVector2;



namespace fastMTT {
  double likelihoodFunc(double *x, double *par);

  enum likelihoodComponent{MET, MASS, PX, PY, ENERGY, IP};
}

class Likelihood{

public:

  Likelihood();

  ~Likelihood();

  double value(const double *x) const;

  void setLeptonInputs(const classic_svFit::LorentzVector & aLeg1P4,
                       const classic_svFit::LorentzVector & aLeg2P4,
                       int aLeg1DecayType, int aLeg2DecayType,
		       int aLeg1DecayMode, int aLeg2DecayMode);

  void setMETInputs(const classic_svFit::LorentzVector & aMET,
                    const TMatrixD& aCovMET);

  void setParameters(const std::vector<double> & parameters);

  void enableComponent(fastMTT::likelihoodComponent aCompIndex);

  void disableComponent(fastMTT::likelihoodComponent aCompIndex);

  double massLikelihood(const double & m) const;

  double ptLikelihood(const double & pTTauTau, int type) const;

  double metTF(const classic_svFit::LorentzVector & metP4,
               const classic_svFit::LorentzVector & nuP4,
               const TMatrixD& covMET) const;

private:

  classic_svFit::LorentzVector leg1P4, leg2P4;
  classic_svFit::LorentzVector recoMET;
  
  TMatrixD covMET;
 
  double mVis, mVisLeg1, mVisLeg2;
  
  int leg1DecayType, leg2DecayType;
  int leg1DecayMode, leg2DecayMode;

  std::vector<double> parameters;

  ///Bit word coding enabled likelihood components
  std::bitset<128> compnentsBitWord;

  //precomputed values used to reduce the number of redundant calculations
  double mVis1OverTauSquare;
  double mVis2OverTauSquare;
  using PowTable = std::array<double, 5u>; //first powers of a number
  //array with dimensions 2x3x5 used to store the powers of the coefficients of two vectors
  std::array<std::array<PowTable,3u>,2u> allpTpows;
  
};
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
class FastMTT {

 public:

  FastMTT();

  ~FastMTT();

  void initialize();

  ///Run fastMTT algorithm for given input
  void run(const std::vector<classic_svFit::MeasuredTauLepton>&,
	   const double &, const double &, const TMatrixD&);

  ///Set likelihood shape parameters. Two parameters are expected:
  ///power of 1/mVis, and scaling factor of mTest
  ///in case of less than two parameters the default values are used:
  ///{6, 1.0/1.15};
  void setLikelihoodParams(const std::vector<double> & aPars);

  ///Set the bit word for switching on particular likelihood components
  void enableComponent(fastMTT::likelihoodComponent aCompIndex);

  ///Set the bit word for switching off particular likelihood components
  void disableComponent(fastMTT::likelihoodComponent aCompIndex);

  ///Retrieve the four momentum corresponding to the likelihood maximum
  const classic_svFit::LorentzVector & getBestP4() const { return bestP4; }

  ///Retrieve the tau1 four momentum corresponding to the likelihood maximum
  const classic_svFit::LorentzVector & getTau1P4() const { return tau1P4; }

  ///Retrieve the tau1 four momentum corresponding to the likelihood maximum
  const classic_svFit::LorentzVector & getTau2P4() const { return tau2P4; }

  ///Retrieve x values corresponding to the likelihood maximim
  std::tuple<double, double> getBestX() const;

  ///Retrieve the best likelihood value
  double getBestLikelihood() const;

  ///Retrieve the likelihood value for given x values
  double getLikelihoodForX(double *x) const;
  
  ///Retrieve the CPU timing for given methods
  ///Possible values:
  /// scan
  /// minimize
  double getCpuTime(const std::string & method);

  ///Retrieve the CPU timing for given methods
  ///Possible values:
  /// scan
  /// minimize
  double getRealTime(const std::string & method);

 private:

  static bool compareLeptons(const classic_svFit::MeasuredTauLepton& measuredTauLepton1,
			     const classic_svFit::MeasuredTauLepton& measuredTauLepton2);
  

  ///Run the minimalization procedure for given inputs.
  ///Results are stored in internal variables accesed by
  ///relevant get methods.
  void minimize();

  ///Run a scan over x1 and x2 [0,1] rectangle for given inputs.
  ///Results are stored in internal variables accesed by
  ///relevant get methods.
  void scan();

   // Minimizer types and algorithms.
   // minimizerName               minimizerAlgorithm
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic
   std::string minimizerName;
   std::string  minimizerAlgorithm;

   ROOT::Math::Minimizer* minimizer;

   ///Minimum location
   std::vector<double> minimumPosition;

   ///Mimimum value
   double minimumValue;
   
   ///Dimension of minimalization space
   unsigned int nVariables;

   ///Names of variables to be minimized
   std::vector<std::string> varNames;
   
  ///Values of variables to be minimized

   std::vector<double> variables;
   
   ///Step sizes for each minimized variable
   std::vector<double> stepSizes;
   
   ROOT::Math::Functor *likelihoodFunctor;
   
   Likelihood myLikelihood;
   
   classic_svFit::LorentzVector tau1P4, tau2P4, bestP4;
   
   TBenchmark clock;
   
   int verbosity;
   
};

#endif
