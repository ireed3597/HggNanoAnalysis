#include <iostream>

using namespace std;

const float mHiggs					= 125;
const float mZ							= 91.19;
const float trans_eta_low 			= 1.4442;
const float trans_eta_high 			= 1.566;
const float mZ_veto_low				= 80;
const float mZ_veto_up				= 100;
const float mllg_window				= 25;

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

//jet selection
const float jet_pt					= 25;
const float jet_eta					= 2.4;
const float jet_neEmEF				= 0.99;
const float jet_neHEF				= 0.99;
const float jet_chHEF				= 0;
const float jet_chEmEF				= 0.99;
const int 	jet_nConstituents		= 1;
const float jet_dR_pho				= 0.4;
const float jet_dR_lep				= 0.4;
const float jet_dR_tau				= 0.4;

//bJet MediumWorkingPoint
const float bTag_medium_WP[3] = { 0.3093 , 0.3033, 0.2770 };

//IsoTracks
const float isoTrk_dR				= 0.2;

//Tree branches
int 			t_run;
int				t_lumiBlock;
int				t_event;
float 			t_MET_pt;
float 			t_MET_phi;
float 			t_weight;
int 			process_id;
int 			category;

float 			mgg;
int				n_electrons;
int				n_muons;
int				n_taus;
int				n_jets;
int				n_bjets;
int 			n_isoTrks;

float			lep12_dphi;
float			lep12_deta;
float			lep12_deta_bdt;
float			lep12_dr;

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
bool 			g1_pixVeto;

float 			g2_ptmgg;
float 			g2_pt;
float 			g2_eta;
float 			g2_eta_bdt;
float 			g2_phi;
float 			g2_idmva;
bool 			g2_pixVeto;

float 			gg_pt;
float 			gg_ptmgg;
float 			gg_eta;
float 			gg_eta_bdt;
float 			gg_phi;
float 			gg_dR;
float 			gg_dPhi;
float 			gg_hel;
float 			gg_hel_phys;
float 			gg_tt_CS;
float 			gg_tt_hel;
float 			gg_tt_hel_phys;

float 			lep1_pt				;
float 			lep1_eta			;
float 			lep1_eta_bdt		;
float 			lep1_phi			;
float 			lep1_charge			;
float 			lep1_pdgID			;
float 			lep1_tightID		;
Int_t 			lep1_id_vs_e		;
Int_t 			lep1_id_vs_m		;
Int_t 			lep1_id_vs_jet		;

float 			lep2_pt;
float 			lep2_eta;
float 			lep2_eta_bdt		;
float 			lep2_phi;
float 			lep2_charge;
float 			lep2_pdgID;
float 			lep2_tightID;
Int_t 			lep2_id_vs_e;
Int_t 			lep2_id_vs_m;
Int_t 			lep2_id_vs_jet;

float 			jet1_pt			;
float 			jet1_eta		;
float 			jet1_eta_bdt	;
float 			jet1_phi		;
float 			jet1_bTag		;
int 			jet1_id			;

float 			jet2_pt			;
float 			jet2_eta		;
float 			jet2_eta_bdt	;
float 			jet2_phi		;
float 			jet2_bTag		;
int 			jet2_id			;

float			max_bTag		;
float			tt_hel					;
float			tt_hel_phys				;

float			m_tautau_vis			;
float			pt_tautau_vis			;
float			eta_tautau_vis			;
float			eta_tautau_vis_bdt		;
float			phi_tautau_vis			;

float			MET_gg_dPhi				;
float			MET_ll_dPhi				;
float			dPhi_MET_l				;
float			ll_dPhi					;
float			ll_dEta					;
float			ll_dR					;
float			m_Z						;

float			dZ						;
float			g1_energyErr	;
float			g2_energyErr	;
float			max_g_ptmgg		;
float			min_g_ptmgg		;
float			max_g_idmva		;
float			min_g_idmva		;
float			lep2_pfRelIso03_all	;
float			lep2_pfRelIso03_chg	;
float			max_lep_pt		;
float			min_lep_pt		;

float 		pt_tautauSVFitLoose		;
float 		eta_tautauSVFitLoose	;
float 		eta_tautauSVFitLoose_bdt;
float 		phi_tautauSVFitLoose	;
float 		m_tautauSVFitLoose		;
float 		dR_tautauSVFitLoose		;
float 		dR_ggtautauSVFitLoose	;
float 		dPhi_tautauSVFitLoose	;
float 		dPhi_ggtautauSVFitLoose	;

float			tau1_pt_SVFit		;
float			tau1_eta_SVFit	;
float			tau1_phi_SVFit	;
float			tau1_m_SVFit		;
float			tau2_pt_SVFit		;
float			tau2_eta_SVFit	;
float			tau2_phi_SVFit	;
float			tau2_m_SVFit		;

float 		pt_tautau_sntMtt		;
float 		eta_tautau_sntMtt	;
float 		eta_tautau_sntMtt_bdt;
float 		phi_tautau_sntMtt	;
float 		m_tautau_sntMtt		;
float 		dR_tautau_sntMtt		;
float 		dR_ggtautau_sntMtt	;
float 		dPhi_tautau_sntMtt	;
float 		dPhi_ggtautau_sntMtt	;

float			tau1_pt_sntMtt		;
float			tau1_eta_sntMtt	;
float			tau1_phi_sntMtt	;
float			tau1_m_sntMtt		;
float			tau2_pt_sntMtt		;
float			tau2_eta_sntMtt	;
float			tau2_phi_sntMtt	;
float			tau2_m_sntMtt		;

float			mX		;
float 		m_llg_lead;
float 		m_llg_subl;
