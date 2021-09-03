#include "../TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "../TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "../TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "../TauAnalysis/ClassicSVfit/interface/FastMTT.h"

double SVfit_mass(float measuredMETx, float measuredMETy, float covMET_XX, float covMET_XY, float covMET_YY,int tauDecay_mode1, int tauDecay_mode2, int category1, int category2, float tau1_pt,float  tau1_eta, float tau1_phi, float tau1_mass, float tau2_pt, float tau2_eta, float tau2_phi, float tau2_mass);

std::vector<double>  SVfit_ditau_p4(float measuredMETx, float measuredMETy, float covMET_XX, float covMET_XY, float covMET_YY,int tauDecay_mode1, int tauDecay_mode2, int category1, int category2, float tau1_pt,float  tau1_eta, float tau1_phi, float tau1_mass, float tau2_pt, float tau2_eta, float tau2_phi, float tau2_mass);

std::vector<std::vector<double>> SVfit_all_p4(float measuredMETx, float measuredMETy, float covMET_XX, float covMET_XY, float covMET_YY,int tauDecay_mode1, int tauDecay_mode2, int category1, int category2, float tau1_pt,float  tau1_eta, float tau1_phi, float tau1_mass, float tau2_pt, float tau2_eta, float tau2_phi, float tau2_mass);
