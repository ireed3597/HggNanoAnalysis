#include "utils.h"
#include "parameters.h"

using namespace tas;

double deltaR( LorentzVector v1 , LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaR( v1 , v2);
}

double deltaR_v1( LorentzVector * v1 , LorentzVector v2){
	LorentzVector v1_1;
	v1_1.SetPt( v1->pt() );
	v1_1.SetEta( v1->eta());
	v1_1.SetPhi( v1->phi());
	v1_1.SetM( v1->mass() );
	return ROOT::Math::VectorUtil::DeltaR( v1_1 , v2);
}

double deltaR_v2( LorentzVector * v1 , LorentzVector * v2){
	LorentzVector v1_1;
	v1_1.SetPt( v1->pt() );
	v1_1.SetEta( v1->eta());
	v1_1.SetPhi( v1->phi());
	v1_1.SetM( v1->mass() );
	LorentzVector v2_1;
	v2_1.SetPt( v2->pt() );
	v2_1.SetEta( v2->eta());
	v2_1.SetPhi( v2->phi());
	v2_1.SetM( v2->mass() );
	return ROOT::Math::VectorUtil::DeltaR( v1_1 , v2_1 );
}

void clear_branches(){

	t_run					= -9;
	t_lumiBlock				= -9;
	t_event					= -9;
	t_MET_pt				= -9;
	t_MET_phi				= -9;
	t_weight				= -9;

	mgg						= -9;
	n_electrons				= -9;
	n_muons					= -9;
	n_taus 					= -9;
	n_jets 					= -9;
	n_isoTrks				= -9;

	g1_ptmgg				= -9;
	g1_pt					= -9;
	g1_eta					= -9;
	g1_phi					= -9;
	g1_idmva				= -9;
	g1_pixVeto				= -9;

	g2_ptmgg				= -9;
	g2_pt					= -9;
	g2_eta					= -9;
	g2_phi					= -9;
	g2_idmva				= -9;
	g2_pixVeto				= -9;

	gg_pt					= -9;
	gg_eta					= -9;
	gg_phi					= -9;
	gg_dR					= -9;

	lep1_pt					= -9;
	lep1_eta				= -9;
	lep1_phi				= -9;
	lep1_charge				= -9;
	lep1_pdgID				= -9;
	lep1_tightID			= -9;
	lep1_id_vs_e			= -9;
	lep1_id_vs_m			= -9;
	lep1_id_vs_jet			= -9;

	lep2_pt        			= -9;
	lep2_eta       			= -9;
	lep2_phi       			= -9;
	lep2_charge    			= -9;
	lep2_pdgID     			= -9;
	lep2_tightID   			= -9;
	lep2_id_vs_e   			= -9;
	lep2_id_vs_m   			= -9;
	lep2_id_vs_jet 			= -9;

	jet1_pt					= -9;
	jet1_eta				= -9;
	jet1_bTag				= -9;
	jet1_id					= -9;

	jet2_pt					= -9;
	jet2_eta				= -9;
	jet2_bTag				= -9;
	jet2_id					= -9;

	pt_tautauSVFitLoose		= -9;
	eta_tautauSVFitLoose	= -9;
	phi_tautauSVFitLoose	= -9;
	m_tautauSVFitLoose		= -9;
	dR_tautauSVFitLoose		= -9;
	dR_ggtautauSVFitLoose	= -9;
	dPhi_MET_tau1			= -9;

	m_tautau_vis			= -9;
	pt_tautau_vis			= -9;
	eta_tautau_vis			= -9;
	phi_tautau_vis			= -9;
}
