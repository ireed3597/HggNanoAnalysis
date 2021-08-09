#include <iostream>

using namespace std;

const float trans_eta_low 			= 1.4442;
const float trans_eta_high 			= 1.566;

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

