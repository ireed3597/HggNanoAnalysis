#include <iostream>

using namespace std;

const float mHiggs					= 125;
const float trans_eta_low 			= 1.4442;
const float trans_eta_high 			= 1.566;
const float mZ_veto_low				= 70;
const float mZ_veto_up				= 110;

//di-photon selection
const float mgg_lower 				= 100;
const float mgg_upper 				= 180;
const float mgg_sideband_lower 		= 120;
const float mgg_sideband_upper 		= 130;
const float lead_pt_mgg_cut			= 0.33;
const float sublead_pt_mgg_cut		= 0.25;
const float pho_idmva_cut			= -0.7;
const float pho_eveto_cut			= 0.5;
const float pho_eta_cut				= 2.5;
const float pho_pt_cut				= 25;

//electron selection
const float ele_pt					= 7;
const float ele_eta					= 2.5;
const float ele_dxy					= 0.045;
const float ele_dz					= 0.2;
const float ele_pfRelIso			= 0.3;
const float ele_dR_pho				= 0.2;

//muon selection
const float muon_pt					= 5;
const float muon_eta				= 2.4;
const float muon_dxy				= 0.045;
const float muon_dz					= 0.2;
const float muon_pfRelIso			= 0.3;
const float muon_dR_pho				= 0.2;

//tau selection
const float tau_pt					= 18;
const float tau_eta					= 2.3;
const float tau_dz					= 0.2;
const int 	tau_deepID_e			= 1;
const int 	tau_deepID_m			= 0;
const int 	tau_deepID_j			= 7;
const float tau_dR_pho				= 0.2;
const float tau_dR_lep				= 0.2;

//IsoTracks
const float isoTrk_dR				= 0.2;

//Tree branches
int 			t_run;
int				t_lumiBlock;
int				t_event;
float 			t_MET_pt;
float 			t_MET_phi;
float 			t_weight;

float 			mgg;
int				n_electrons;
int				n_muons;
int				n_taus;
int				n_jets;
int 			n_isoTrks;

bool			cat1;
bool			cat2;
bool			cat3;
bool			cat4;
bool			cat5;
bool			cat6;
bool			cat7;
bool			cat8;

float 			g1_ptmgg;
float 			g1_pt;
float 			g1_eta;
float 			g1_eta_bdt;
float 			g1_phi;
float 			g1_idmva;
float 			g1_pixVeto;

float 			g2_ptmgg;
float 			g2_pt;
float 			g2_eta;
float 			g2_eta_bdt;
float 			g2_phi;
float 			g2_idmva;
float 			g2_pixVeto;

float 			gg_pt;
float 			gg_eta;
float 			gg_eta_bdt;
float 			gg_phi;
float 			gg_dR;
float 			gg_dPhi;

float 			lep1_pt				;
float 			lep1_eta			;
float 			lep1_eta_bdt		;
float 			lep1_phi			;
float 			lep1_charge			;
float 			lep1_pdgID			;
float 			lep1_tightID		;
float 			lep1_id_vs_e		;
float 			lep1_id_vs_m		;
float 			lep1_id_vs_jet		;

float 			lep2_pt;
float 			lep2_eta;
float 			lep2_eta_bdt		;
float 			lep2_phi;
float 			lep2_charge;
float 			lep2_pdgID;
float 			lep2_tightID;
float 			lep2_id_vs_e;
float 			lep2_id_vs_m;
float 			lep2_id_vs_jet;

float 			jet1_pt			;
float 			jet1_eta		;
float 			jet1_eta_bdt	;
float 			jet1_bTag		;
float 			jet1_id			;

float 			jet2_pt			;
float 			jet2_eta		;
float 			jet2_eta_bdt	;
float 			jet2_bTag		;
float 			jet2_id			;

float 			pt_tautauSVFitLoose		;
float 			eta_tautauSVFitLoose	;
float 			eta_tautauSVFitLoose_bdt;
float 			phi_tautauSVFitLoose	;
float 			m_tautauSVFitLoose		;
float 			dR_tautauSVFitLoose		;
float 			dR_ggtautauSVFitLoose	;
float 			dPhi_MET_tau1			;

float			m_tautau_vis			;
float			pt_tautau_vis			;
float			eta_tautau_vis			;
float			eta_tautau_vis_bdt		;
float			phi_tautau_vis			;
