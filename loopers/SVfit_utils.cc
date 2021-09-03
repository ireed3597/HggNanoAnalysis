//#include "SVfit_utils.h"
#include "../TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "../TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "../TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "../TauAnalysis/ClassicSVfit/interface/FastMTT.h"

//using namespace classic_svFit;
using namespace std;

double SVfit_mass(float measuredMETx, float measuredMETy, float covMET_XX, float covMET_XY, float covMET_YY,int tauDecay_mode1, int tauDecay_mode2, int category1, int category2, float tau1_pt,float  tau1_eta, float tau1_phi, float tau1_mass, float tau2_pt, float tau2_eta, float tau2_phi, float tau2_mass)
{
	double diTauSVfit_mass;

	TMatrixD covMET(2, 2);
	covMET[0][0] = covMET_XX;
	covMET[0][1] = covMET_XY;
	covMET[1][0] = covMET_XY;
	covMET[1][1] = covMET_YY;

	// define lepton four vectors
	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

	if ( category1 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau1_pt, tau1_eta, tau1_phi, 0.10566));
	}
	else if ( category1 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau1_pt, tau1_eta, tau1_phi, 0.51100e-3));
	}
	else if ( category1 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau1_pt, tau1_eta, tau1_phi, tau1_mass, tauDecay_mode1));
	}

	if ( category2 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau2_pt, tau2_eta, tau2_phi, 0.10566));
	}
	else if ( category2 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau2_pt, tau2_eta, tau2_phi, 0.51100e-3));
	}
	else if ( category2 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau2_pt, tau2_eta, tau2_phi, tau2_mass, tauDecay_mode2));
	}
	//measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  tau2_pt, tau2_eta, tau2_phi,  tau2_mass, tauDecay_mode2));

	//Run FastMTT
	FastMTT aFastMTTAlgo;
	aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
	classic_svFit::LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
	diTauSVfit_mass = ttP4.M();

	return diTauSVfit_mass;
}


classic_svFit::LorentzVector  SVfit_ditau_p4(float measuredMETx, float measuredMETy, float covMET_XX, float covMET_XY, float covMET_YY,int tauDecay_mode1, int tauDecay_mode2, int category1, int category2, float tau1_pt,float  tau1_eta, float tau1_phi, float tau1_mass, float tau2_pt, float tau2_eta, float tau2_phi, float tau2_mass)
{
	double diTauSVfit_mass;

	TMatrixD covMET(2, 2);
	covMET[0][0] = covMET_XX;
	covMET[0][1] = covMET_XY;
	covMET[1][0] = covMET_XY;
	covMET[1][1] = covMET_YY;

	// define lepton four vectors
	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

	if ( category1 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau1_pt, tau1_eta, tau1_phi, 0.10566));
	}
	else if ( category1 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau1_pt, tau1_eta, tau1_phi, 0.51100e-3));
	}
	else if ( category1 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau1_pt, tau1_eta, tau1_phi, tau1_mass, tauDecay_mode1));
	}

	if ( category2 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau2_pt, tau2_eta, tau2_phi, 0.10566));
	}
	else if ( category2 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau2_pt, tau2_eta, tau2_phi, 0.51100e-3));
	}
	else if ( category2 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau2_pt, tau2_eta, tau2_phi, tau2_mass, tauDecay_mode2));
	}
	//measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  tau2_pt, tau2_eta, tau2_phi,  tau2_mass, tauDecay_mode2));

	//Run FastMTT
	FastMTT aFastMTTAlgo;
	aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
	classic_svFit::LorentzVector ttP4 = aFastMTTAlgo.getBestP4();

	return ttP4;
}

vector<double>  SVfit_ditau_p4_v2(float measuredMETx, float measuredMETy, float covMET_XX, float covMET_XY, float covMET_YY,int tauDecay_mode1, int tauDecay_mode2, int category1, int category2, float tau1_pt,float  tau1_eta, float tau1_phi, float tau1_mass, float tau2_pt, float tau2_eta, float tau2_phi, float tau2_mass)
{
	double diTauSVfit_mass;

	TMatrixD covMET(2, 2);
	covMET[0][0] = covMET_XX;
	covMET[0][1] = covMET_XY;
	covMET[1][0] = covMET_XY;
	covMET[1][1] = covMET_YY;

	// define lepton four vectors
	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

	if ( category1 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau1_pt, tau1_eta, tau1_phi, 0.10566));
	}
	else if ( category1 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau1_pt, tau1_eta, tau1_phi, 0.51100e-3));
	}
	else if ( category1 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau1_pt, tau1_eta, tau1_phi, tau1_mass, tauDecay_mode1));
	}

	if ( category2 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau2_pt, tau2_eta, tau2_phi, 0.10566));
	}
	else if ( category2 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau2_pt, tau2_eta, tau2_phi, 0.51100e-3));
	}
	else if ( category2 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau2_pt, tau2_eta, tau2_phi, tau2_mass, tauDecay_mode2));
	}
	//measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  tau2_pt, tau2_eta, tau2_phi,  tau2_mass, tauDecay_mode2));

	//Run FastMTT
	FastMTT aFastMTTAlgo;
	aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
	classic_svFit::LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
	//diTauSVfit_mass = ttP4.M();
	
	vector<double> tt_p4(4,0.);
	tt_p4[0]	= ttP4.Px();
	tt_p4[1]	= ttP4.Py();
	tt_p4[2]	= ttP4.Pz();
	tt_p4[3]	= ttP4.E();
	return tt_p4;
}


vector<classic_svFit::LorentzVector>  SVfit_all_p4(float measuredMETx, float measuredMETy, float covMET_XX, float covMET_XY, float covMET_YY,int tauDecay_mode1, int tauDecay_mode2, int category1, int category2, float tau1_pt,float  tau1_eta, float tau1_phi, float tau1_mass, float tau2_pt, float tau2_eta, float tau2_phi, float tau2_mass)
{
	double diTauSVfit_mass;

	TMatrixD covMET(2, 2);
	covMET[0][0] = covMET_XX;
	covMET[0][1] = covMET_XY;
	covMET[1][0] = covMET_XY;
	covMET[1][1] = covMET_YY;

	// define lepton four vectors
	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

	if ( category1 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau1_pt, tau1_eta, tau1_phi, 0.10566));
	}
	else if ( category1 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau1_pt, tau1_eta, tau1_phi, 0.51100e-3));
	}
	else if ( category1 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau1_pt, tau1_eta, tau1_phi, tau1_mass, tauDecay_mode1));
	}

	if ( category2 == 1 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, tau2_pt, tau2_eta, tau2_phi, 0.10566));
	}
	else if ( category2 == 2 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, tau2_pt, tau2_eta, tau2_phi, 0.51100e-3));
	}
	else if ( category2 == 3 ){
		measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, tau2_pt, tau2_eta, tau2_phi, tau2_mass, tauDecay_mode2));
	}
	//measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  tau2_pt, tau2_eta, tau2_phi,  tau2_mass, tauDecay_mode2));

	//Run FastMTT
	FastMTT aFastMTTAlgo;
	aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
	classic_svFit::LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
	//diTauSVfit_mass = ttP4.M();
	
	//classic_svFit::LorentzVector tt_p4;
	//tt_p4[0]	= ttP4.Px();
	//tt_p4[1]	= ttP4.Py();
	//tt_p4[2]	= ttP4.Pz();
	//tt_p4[3]	= ttP4.E();

	classic_svFit::LorentzVector tau1P4 = aFastMTTAlgo.getTau1P4();	
	//vector<double> tau1_p4(4,0.);
	//tau1_p4[0]	= tau1P4.Px();
	//tau1_p4[1]	= tau1P4.Py();
	//tau1_p4[2]	= tau1P4.Pz();
	//tau1_p4[3]	= tau1P4.E();

	classic_svFit::LorentzVector tau2P4 = aFastMTTAlgo.getTau2P4();	
	//vector<double> tau2_p4(4,0.);
	//tau2_p4[0]	= tau2P4.Px();
	//tau2_p4[1]	= tau2P4.Py();
	//tau2_p4[2]	= tau2P4.Pz();
	//tau2_p4[3]	= tau2P4.E();
	
	return  { ttP4 , tau1P4 , tau2P4 } ;
}
