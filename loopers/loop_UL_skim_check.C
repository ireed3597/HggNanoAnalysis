#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TObjString.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"

#include "../NanoCORE/Nano.cc"
#include "../NanoCORE/tqdm.h"
#include "../NanoCORE/utils.cc"
#include "SVfit_utils.cc"

#include <string>
#include <iostream>
#include <iomanip>

#define SUM(vec) std::accumulate((vec).begin(), (vec).end(), 0);
#define SUM_GT(vec,num) std::accumulate((vec).begin(), (vec).end(), 0, [](float x,float y){return ((y > (num)) ? x+y : x); });
#define COUNT_GT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x > (num); });
#define COUNT_LT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x < (num); });

#define H1(name,nbins,low,high) TH1F *h_##name = new TH1F(#name,#name,nbins,low,high);

// #define DEBUG

struct debugger { template<typename T> debugger& operator , (const T& v) { cerr<<v<<" "; return *this; } } dbg;
#ifdef DEBUG
    #define debug(args...) do {cerr << #args << ": "; dbg,args; cerr << endl;} while(0)
#else
    #define debug(args...)
#endif

using namespace std;
using namespace tas;

int ScanChain( TChain *ch, string proc, string str_year, float scale_factor = 1, bool resonant = false ) {

	int year;
	if ( str_year == "2016_APV") year = 2016;
	else { year = stoi(str_year); }

	TString file_name = proc + "_skimCheck_" +  str_year;

	TFile* f1 = new TFile("outputs_UL/" + file_name + ".root", "RECREATE");
	H1(mgg, 60, 100 , 180 );
	H1(mgg_1t0l, 60, 100 , 180 );
	H1(mgg_1t1l, 60, 100 , 180 );
	H1(mgg_2t0l, 60, 100 , 180 );
	H1(mgg_0t2l, 60, 100 , 180 );
	H1(mgg_1t0l_iso, 60, 100 , 180 );

	H1(iso_pt, 60, 0 , 50 );
	H1(iso_dxy, 60, -0.3 , 0.3 );
	H1(iso_dz, 60, -0.3 , 0.3 );
	H1(iso_eta, 60, -3 , 3 );
	H1(iso_fromPV, 3, -1 , 1 );
	H1(iso_miniPFRelIso_all, 60, -1.5 , 1.5 );
	H1(iso_miniPFRelIso_chg, 60, -1.5 , 1.5 );
	H1(iso_pdgID, 425, -212 , 212 );
	H1(iso_pfRelIso03_all, 60, -1.5 , 1.5 );
	H1(iso_pfRelIso03_chg, 60, -1.5 , 1.5 );
	H1(iso_phi, 60, -3.5 , 3.5 );
	H1(iso_nIso, 10, 0 , 11);

	TTree *out_tree	=	new TTree("Events","output tree");
	out_tree->Branch("year"					,	&year					,	"year/I"			);
	out_tree->Branch("run"					,	&t_run					,	"run/I"				);
	out_tree->Branch("lumiBlock"			,	&t_lumiBlock			,	"lumiBlock/I"		);
	out_tree->Branch("event"				,	&t_event				,	"event/I"			);
	out_tree->Branch("process_id"			,	&process_id				,	"process_id/I"		);
	out_tree->Branch("MET_pt"				,	&t_MET_pt				,	"MET_pt/F"			);
	out_tree->Branch("MET_phi"				,	&t_MET_phi				,	"MET_phi/F"			);
	out_tree->Branch("weight"				,	&t_weight				,	"weight/F"			);
	out_tree->Branch("Category"				,	&category				,	"Category/I"		);

	out_tree->Branch("cat1"					,	&cat1					,	"cat1/B"			);
	out_tree->Branch("cat2"					,	&cat2					,	"cat2/B"			);
	out_tree->Branch("cat3"					,	&cat3					,	"cat3/B"			);
	out_tree->Branch("cat4"					,	&cat4					,	"cat4/B"			);
	out_tree->Branch("cat5"					,	&cat5					,	"cat5/B"			);
	out_tree->Branch("cat6"					,	&cat6					,	"cat6/B"			);
	out_tree->Branch("cat7"					,	&cat7					,	"cat7/B"			);
	out_tree->Branch("cat8"					,	&cat8					,	"cat8/B"			);

	out_tree->Branch("n_electrons"			,	&n_electrons			,	"n_electrons/I"		);
	out_tree->Branch("n_muons"				,	&n_muons				,	"n_muons/I"			);
	out_tree->Branch("n_taus"				,	&n_taus					,	"n_taus/I"			);
	out_tree->Branch("n_isoTrks"			,	&n_isoTrks				,	"n_isoTrks/I"		);
	out_tree->Branch("n_jets"				,	&n_jets					,	"n_jets/I"			);
	out_tree->Branch("n_bjets"				,	&n_bjets				,	"n_bjets/I"			);

	out_tree->Branch("m_Z"						,	&m_Z					, 	"m_Z/F"				);	  
	out_tree->Branch("mX"							,	&mX						, 	"mX/F"				);	  
	out_tree->Branch("m_llg_lead"			,	&m_llg_lead		, 	"m_llg_lead/F"				);	  
	out_tree->Branch("m_llg_subl"			,	&m_llg_subl		, 	"m_llg_subl/F"				);	  

	out_tree->Branch("MET_gg_dPhi"			,	&MET_gg_dPhi			, 	"MET_gg_dPhi/F"		);	  
	out_tree->Branch("MET_ll_dPhi"			,	&MET_ll_dPhi			, 	"MET_ll_dPhi/F"		);	  
	out_tree->Branch("dPhi_MET_l"			,	&dPhi_MET_l				, 	"dPhi_MET_l/F"		);	  

	out_tree->Branch("lep12_dphi"			,	&lep12_dphi				, 	"lep12_dphi/F"		);	  
	out_tree->Branch("lep12_deta"			,	&lep12_deta				, 	"lep12_deta/F"		);	  
	out_tree->Branch("lep12_deta_bdt"		,	&lep12_deta_bdt			, 	"lep12_deta_bdt/F"	);	  
	out_tree->Branch("lep12_dr"				,	&lep12_dr				, 	"lep12_dr/F"		);	  

	out_tree->Branch("g1_ptmgg"				,	&g1_ptmgg				,	"g1_ptmgg/F"		);
	out_tree->Branch("g1_pt"				,	&g1_pt					, 	"g1_pt/F"			);
	out_tree->Branch("g1_eta"				,	&g1_eta					, 	"g1_eta/F"			);	  
	out_tree->Branch("g1_eta_bdt"			,	&g1_eta_bdt				, 	"g1_eta_bdt/F"		);	  
	out_tree->Branch("g1_phi"				,	&g1_phi					,  	"g1_phi/F"			);
	out_tree->Branch("g1_idmva"				,	&g1_idmva				,	"g1_idmva/F"		);
	out_tree->Branch("g1_pixVeto"			,	&g1_pixVeto				,	"g1_pixVeto/B"		);
	out_tree->Branch("g1_energyErr"			,	&g1_energyErr			,  	"g1_energyErr/F"	);
	out_tree->Branch("g2_ptmgg"				,	&g2_ptmgg				,	"g2_ptmgg/F"		);
	out_tree->Branch("g2_pt"				,	&g2_pt					, 	"g2_pt/F"			);
	out_tree->Branch("g2_eta"				,	&g2_eta					, 	"g2_eta/F"			);	  
	out_tree->Branch("g2_eta_bdt"			,	&g2_eta_bdt				, 	"g2_eta_bdt/F"		);	  
	out_tree->Branch("g2_phi"				,	&g2_phi					,  	"g2_phi/F"			);
	out_tree->Branch("g2_idmva"				,	&g2_idmva				,	"g2_idmva/F"		);
	out_tree->Branch("g2_pixVeto"			,	&g2_pixVeto				,	"g2_pixVeto/B"		);
	out_tree->Branch("g2_energyErr"			,	&g2_energyErr			,  	"g2_energyErr/F"	);

	out_tree->Branch("max_g_ptmgg"			,	&max_g_ptmgg			,	"max_g_ptmgg/F"		);
	out_tree->Branch("min_g_ptmgg"			,	&min_g_ptmgg			,	"min_g_ptmgg/F"		);
	out_tree->Branch("max_g_idmva"			,	&max_g_idmva			,	"max_g_idmva/F"		);
	out_tree->Branch("min_g_idmva"			,	&min_g_idmva			,	"min_g_idmva/F"		);

	out_tree->Branch("gg_pt"				,	&gg_pt					, 	"gg_pt/F"			);
	out_tree->Branch("gg_ptmgg"				,	&gg_ptmgg				, 	"gg_ptmgg/F"		);
	out_tree->Branch("gg_eta"				,	&gg_eta					, 	"gg_eta/F"			);	  
	out_tree->Branch("gg_eta_bdt"			,	&gg_eta_bdt				, 	"gg_eta_bdt/F"		);	  
	out_tree->Branch("gg_phi"				,	&gg_phi					,  	"gg_phi/F"			);
	out_tree->Branch("gg_dR"				,	&gg_dR					,	"gg_dR/F"			);
	out_tree->Branch("gg_dPhi"				,	&gg_dPhi				,	"gg_dPhi/F"			);
	out_tree->Branch("gg_hel"				,	&gg_hel					,	"gg_hel/F"			);
	out_tree->Branch("gg_hel_phys"			,	&gg_hel_phys			,	"gg_hel_phys/F"		);
	out_tree->Branch("gg_tt_CS"				,	&gg_tt_CS				,	"gg_tt_CS/F"		);
	out_tree->Branch("gg_tt_hel"			,	&gg_tt_hel				,	"gg_tt_hel/F"		);
	out_tree->Branch("gg_tt_hel_phys"		,	&gg_tt_hel_phys			,	"gg_tt_hel_phys/F"	);
	out_tree->Branch("mgg"					,	&mgg					,	"mgg/F"				);
	out_tree->Branch("CMS_hgg_mass"			,	&mgg					,	"CMS_hgg_mass/F"	);

	out_tree->Branch("lep1_pt"				,	&lep1_pt				, 	"lep1_pt/F"					);	  
	out_tree->Branch("lep1_eta"				,	&lep1_eta				,  	"lep1_eta/F"				);
	out_tree->Branch("lep1_eta_bdt"			,	&lep1_eta_bdt			,  	"lep1_eta_bdt/F"			);
	out_tree->Branch("lep1_phi"				,	&lep1_phi				,	"lep1_phi/F"				);
	out_tree->Branch("lep1_charge"			,	&lep1_charge			, 	"lep1_charge/I"				);	  
	out_tree->Branch("lep1_pdgID"			,	&lep1_pdgID				,  	"lep1_pdgID/F"				);
	out_tree->Branch("lep1_tightID"			,	&lep1_tightID			,	"lep1_tightID/F"			);
	out_tree->Branch("lep1_id_vs_e"			,	&lep1_id_vs_e			, 	"lep1_id_vs_e/I");	  
	out_tree->Branch("lep1_id_vs_m"			,	&lep1_id_vs_m			,  	"lep1_id_vs_m/I");
	out_tree->Branch("lep1_id_vs_jet"		,	&lep1_id_vs_jet			,	"lep1_id_vs_jet/I");
	out_tree->Branch("lep2_pt"				,	&lep2_pt				, 	"lep2_pt/F"					);	  
	out_tree->Branch("lep2_eta"				,	&lep2_eta				,  	"lep2_eta/F"				);
	out_tree->Branch("lep2_eta_bdt"			,	&lep2_eta_bdt			,  	"lep2_eta_bdt/F"			);
	out_tree->Branch("lep2_phi"				,	&lep2_phi				,	"lep2_phi/F"				);
	out_tree->Branch("lep2_charge"			,	&lep2_charge			, 	"lep2_charge/I"				);	  
	out_tree->Branch("lep2_pdgID"			,	&lep2_pdgID				,  	"lep2_pdgID/F"				);
	out_tree->Branch("lep2_tightID"			,	&lep2_tightID			,	"lep2_tightID/F"			);
	out_tree->Branch("lep2_id_vs_e"			,	&lep2_id_vs_e			, 	"lep2_id_vs_e/I"			);	  
	out_tree->Branch("lep2_id_vs_m"			,	&lep2_id_vs_m			,  	"lep2_id_vs_m/I"			);
	out_tree->Branch("lep2_id_vs_jet"		,	&lep2_id_vs_jet			,	"lep2_id_vs_jet/I"			);
	out_tree->Branch("lep2_pfRelIso03_all"	,	&lep2_pfRelIso03_all	,	"lep2_pfRelIso03_all/F"		);
	out_tree->Branch("lep2_pfRelIso03_chg"	,	&lep2_pfRelIso03_chg	,	"lep2_pfRelIso03_chg/F"		);
	out_tree->Branch("max_lep_pt"			,	&max_lep_pt				, 	"max_lep_pt/F"				);	  
	out_tree->Branch("min_lep_pt"			,	&min_lep_pt				, 	"min_lep_pt/F"				);	  

	out_tree->Branch("jet1_pt"				,	&jet1_pt				, 	"jet1_pt/F"					);	  
	out_tree->Branch("jet1_eta"				,	&jet1_eta				, 	"jet1_eta/F"				);	  
	out_tree->Branch("jet1_eta_bdt"			,	&jet1_eta_bdt			,  	"jet1_eta_bdt/F"			);
	out_tree->Branch("jet1_phi"				,	&jet1_phi				, 	"jet1_phi/F"				);	  
	out_tree->Branch("jet1_bTag"			,	&jet1_bTag				, 	"jet1_bTag/F"				);	  
	out_tree->Branch("jet1_id"				,	&jet1_id				, 	"jet1_id/I"					);	  
	out_tree->Branch("jet2_pt"				,	&jet2_pt				, 	"jet2_pt/F"					);	  
	out_tree->Branch("jet2_eta"				,	&jet2_eta				, 	"jet2_eta/F"				);	  
	out_tree->Branch("jet2_eta_bdt"			,	&jet2_eta_bdt			,  	"jet2_eta_bdt/F"			);
	out_tree->Branch("jet2_phi"				,	&jet2_phi				, 	"jet2_phi/F"				);	  
	out_tree->Branch("jet2_bTag"			,	&jet2_bTag				, 	"jet2_bTag/F"				);	  
	out_tree->Branch("jet2_id"				,	&jet2_id				, 	"jet2_id/I"					);	  

	out_tree->Branch("max_bTag"				,	&max_bTag				, 	"max_bTag/F"				);	  

	out_tree->Branch("tau1_pt_SVFit"		,	&tau1_pt_SVFit			,	"tau1_pt_SVFit/F"			);	  
	out_tree->Branch("tau1_eta_SVFit"		,	&tau1_eta_SVFit			,	"tau1_eta_SVFit/F"			);	  
	out_tree->Branch("tau1_phi_SVFit"		,	&tau1_phi_SVFit			,	"tau1_phi_SVFit/F"			);	  
	out_tree->Branch("tau1_m_SVFit"			,	&tau1_m_SVFit			,	"tau1_m_SVFit/F"			);	  
	out_tree->Branch("tau2_pt_SVFit"		,	&tau2_pt_SVFit			,	"tau2_pt_SVFit/F"			);	  
	out_tree->Branch("tau2_eta_SVFit"		,	&tau2_eta_SVFit			,	"tau2_eta_SVFit/F"			);	  
	out_tree->Branch("tau2_phi_SVFit"		,	&tau2_phi_SVFit			,	"tau2_phi_SVFit/F"			);	  
	out_tree->Branch("tau2_m_SVFit"			,	&tau2_m_SVFit			,	"tau2_m_SVFit/F"			);	  

	out_tree->Branch("pt_tautau_SVFit"		,	&pt_tautauSVFitLoose		,	"pt_tautau_SVFit/F"			);	  
	out_tree->Branch("eta_tautau_SVFit"		,	&eta_tautauSVFitLoose		, 	"eta_tautau_SVFit/F"		);	  
	out_tree->Branch("eta_tautau_SVFit_bdt"	,	&eta_tautauSVFitLoose_bdt	, 	"eta_tautau_SVFit_bdt/F"		);	  
	out_tree->Branch("phi_tautau_SVFit"		,	&phi_tautauSVFitLoose		, 	"phi_tautau_SVFit/F"		);	  
	out_tree->Branch("m_tautau_SVFit"		,	&m_tautauSVFitLoose			, 	"m_tautau_SVFit/F"			);	  
	out_tree->Branch("dR_tautau_SVFit"		,	&dR_tautauSVFitLoose		, 	"dR_tautau_SVFit/F"			);	  
	out_tree->Branch("dR_ggtautau_SVFit"	,	&dR_ggtautauSVFitLoose		, 	"dR_ggtautau_SVFit/F"		);	  
	out_tree->Branch("dPhi_tautau_SVFit"	,	&dPhi_tautauSVFitLoose		, 	"dPhi_tautau_SVFit/F"		);
	out_tree->Branch("dPhi_ggtautau_SVFit"	,	&dPhi_ggtautauSVFitLoose	, 	"dPhi_ggtautau_SVFit/F"		);	  

	out_tree->Branch("tt_hel"				,	&tt_hel  					, 	"tt_hel/F"					);	  
	out_tree->Branch("tt_hel_phys"			,	&tt_hel_phys				, 	"tt_hel_phys/F"				);	  
	out_tree->Branch("m_tautau_vis"			,	&m_tautau_vis  				, 	"m_tautau_vis/F"			);	  
	out_tree->Branch("pt_tautau_vis"		,	&pt_tautau_vis 				, 	"pt_tautau_vis/F"			);	  
	out_tree->Branch("eta_tautau_vis"		,	&eta_tautau_vis				, 	"eta_tautau_vis/F"			);	  
	out_tree->Branch("eta_tautau_vis_bdt"	,	&eta_tautau_vis_bdt			, 	"eta_tautau_vis_bdt/F"		);	  
	out_tree->Branch("phi_tautau_vis"		,	&phi_tautau_vis				, 	"phi_tautau_vis/F"			);	  

	out_tree->Branch("dZ"					,	&dZ  						, 	"dZ/F"						);	  

	bool ggf_samples = false;
	//define process ids
	if (proc.find(std::string("HH_ggWW_semileptonic")) != std::string::npos)	process_id = -4;
	else if (proc.find(std::string("HH_ggWW_dileptonic")) != std::string::npos) 	process_id = -3;
	//else if (proc.find(std::string("HH_ggZZ")) != std::string::npos) 				process_id = -2;
	else if (proc.find(std::string("HH_ggZZ_2l2q")) != std::string::npos) 				process_id = -6;
	else if (proc.find(std::string("HH_ggZZ_4l")) != std::string::npos) 				process_id = -5;
	else if (proc.find(std::string("HH_ggTauTau")) != std::string::npos) 			process_id = -1;
	else if (proc.find(std::string("Data")) != std::string::npos) 				process_id = 0;
	else if (proc.find(std::string("ZGamma")) != std::string::npos) 				process_id = 2;
	else if (proc.find(std::string("DiPhoton")) != std::string::npos) 			process_id = 3;
	else if (proc.find(std::string("WGamma")) != std::string::npos) 				process_id = 4;
	else if (proc.find(std::string("TTbar")) != std::string::npos) 				process_id = 5;
	else if (proc.find(std::string("TTGamma")) != std::string::npos) 				process_id = 6;
	else if (proc.find(std::string("TTGG")) != std::string::npos) 				process_id = 7;
	else if (proc.find(std::string("GJets")) != std::string::npos)				process_id = 8;
	else if (proc.find(std::string("VH")) != std::string::npos) 					process_id = 9;
	else if (proc.find(std::string("ttH")) != std::string::npos) 					process_id = 10;
	else if (proc.find(std::string("ggH")) != std::string::npos) 					process_id = 11;
	else if (proc.find(std::string("VBFH")) != std::string::npos) 				process_id = 12;
	else if (proc.find(std::string("WW")) != std::string::npos) 				process_id = 14;
	else if (proc.find(std::string("WZ")) != std::string::npos) 				process_id = 15;
	else if (proc.find(std::string("ZZ")) != std::string::npos) 				process_id = 16;
	//non-SM couplings
	else if (proc.find(std::string("ggf_c0_HHggtautau")) != std::string::npos) 		{		process_id = -10;		ggf_samples = true;		}
	else if (proc.find(std::string("ggf_c1_HHggtautau")) != std::string::npos) 		{		process_id = -11;		ggf_samples = true;		}
	else if (proc.find(std::string("ggf_c2p45_HHggtautau")) != std::string::npos) {				process_id = -12;		ggf_samples = true;		}
	else if (proc.find(std::string("ggf_c5_HHggtautau")) != std::string::npos) 		{		process_id = -13;		ggf_samples = true;		}
	else if (proc.find(std::string("vbf_cv_0p5_c2v_1_c3_1_HHggtautau")) != std::string::npos) 				process_id = -14;
	else if (proc.find(std::string("vbf_cv_1p5_c2v_1_c3_1_HHggtautau")) != std::string::npos) 				process_id = -15;
	else if (proc.find(std::string("vbf_cv_1_c2v_0_c3_1_HHggtautau")) != std::string::npos) 				process_id = -16;
	else if (proc.find(std::string("vbf_cv_1_c2v_1_c3_0_HHggtautau")) != std::string::npos) 				process_id = -17;
	else if (proc.find(std::string("vbf_cv_1_c2v_1_c3_1_HHggtautau")) != std::string::npos) 				process_id = -18;
	else if (proc.find(std::string("vbf_cv_1_c2v_1_c3_2_HHggtautau")) != std::string::npos) 				process_id = -19;
	else if (proc.find(std::string("vbf_cv_1_c2v_2_c3_1_HHggtautau")) != std::string::npos) 				process_id = -20;

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();
    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    tqdm bar;

    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile *file = TFile::Open( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");

        tree->SetCacheSize(128*1024*1024);
        tree->SetCacheLearnEntries(100);

        //auto psRead = new TTreePerfStats("readPerf", tree);
				nt.SetYear(year);
        nt.Init(tree);

        for( unsigned int loop_event = 0; loop_event < tree->GetEntriesFast(); ++loop_event) {

         nt.GetEntry(loop_event);

         nEventsTotal++;
         bar.progress(nEventsTotal, nEventsChain);

				clear_branches();

				if ( proc != "Data" && ggf_samples & ( fabs(genWeight()) >= 1 ) )  continue;

				//////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////							photon selection
				vector<int> pho_cands;
				vector<float> pho_pt_cands;
				for (unsigned int i=0; i<nPhoton(); i++){
					if ( Photon_electronVeto().at(i)  >=0.5 &&
						 (	(	Photon_isScEtaEB().at(i)	 && Photon_r9().at(i) > 0.85 )			//pho_EB_highR9
						||	(	Photon_isScEtaEE().at(i)	 && Photon_r9().at(i) > 0.90 )			//pho_EE_highR9
						||	(	Photon_isScEtaEB().at(i)	 && Photon_r9().at(i) < 0.85 && Photon_r9().at(i) > 0.5 && Photon_sieie().at(i) < 0.015 && Photon_trkSumPtHollowConeDR03().at(i) < 6.0  && ( Photon_photonIso().at(i) - 0.16544*fixedGridRhoAll() ) < 4.0 )			//pho_EB_lowR9
						||	(	Photon_isScEtaEE().at(i)	 && Photon_r9().at(i) < 0.90 && Photon_r9().at(i) > 0.8 && Photon_sieie().at(i) < 0.035 && Photon_trkSumPtHollowConeDR03().at(i) < 6.0  && ( Photon_photonIso().at(i) - 0.13212*fixedGridRhoAll() ) < 4.0 )			/*pho_EE_lowR9 */ )
//						||	(	Photon_isScEtaEB().at(i)	 && Photon_r9().at(i) < 0.85 && Photon_r9().at(i) > 0.5 && Photon_sieie().at(i) < 0.015 && Photon_trkSumPtHollowConeDR03().at(i) < 6.0  && ( Photon_pfPhoIso03().at(i) - 0.16544*fixedGridRhoFastjetAll() ) < 4.0 )			//pho_EB_lowR9
//						||	(	Photon_isScEtaEE().at(i)	 && Photon_r9().at(i) < 0.90 && Photon_r9().at(i) > 0.8 && Photon_sieie().at(i) < 0.035 && Photon_trkSumPtHollowConeDR03().at(i) < 6.0  && ( Photon_pfPhoIso03().at(i) - 0.13212*fixedGridRhoFastjetAll() ) < 4.0 )			/*pho_EE_lowR9 */ )
					 	&& 	Photon_hoe().at(i) < 0.08
						&&	fabs(Photon_eta().at(i)) < 2.5
						&&	(	fabs(Photon_eta().at(i)) < 1.442 || fabs(Photon_eta().at(i)) > 1.566 )
						&&	( Photon_r9().at(i) > 0.8 || Photon_chargedHadronIso().at(i) < 20 || Photon_chargedHadronIso().at(i) / Photon_pt().at(i) < 0.3 )
						&&  ( Photon_isScEtaEB().at(i) || Photon_isScEtaEE().at(i) )
						&&    Photon_mvaID().at(i) > pho_idmva_cut
						){
							pho_cands.push_back(i);					
							pho_pt_cands.push_back( Photon_pt().at(i) );					
					}
				}
				if ( pho_cands.size() < 2 ) continue;

				sort(pho_pt_cands.begin(), pho_pt_cands.end(), greater<float>());
				int gHidx[2] = {-1,-1};
				for (unsigned int i=0; i<pho_cands.size();i++){
					if ( Photon_pt().at(pho_cands[i]) == pho_pt_cands[0] && Photon_pt().at(pho_cands[i]) > 35 ) gHidx[0]	= pho_cands[i];
					if ( Photon_pt().at(pho_cands[i]) == pho_pt_cands[1] && Photon_pt().at(pho_cands[i]) > 25 ) gHidx[1]	= pho_cands[i];
				}
				if ( gHidx[0] < 0 || gHidx[1] < 0 ) continue;

			//di-photon selection
			mgg = (float)(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]) ).M();
			if ( mgg < mgg_lower || mgg > mgg_upper ) continue;
			if ( (proc == "Data" || !resonant) && mgg > mgg_sideband_lower && mgg < mgg_sideband_upper ) continue;

			if ( Photon_pt().at(gHidx[0]) < pho_pt_cut 				|| Photon_pt().at(gHidx[1]) < pho_pt_cut ) continue;
			if ( Photon_pt().at(gHidx[0]) / mgg < lead_pt_mgg_cut 	|| Photon_pt().at(gHidx[1]) / mgg < sublead_pt_mgg_cut ) continue;
			if ( Photon_mvaID().at(gHidx[0]) < pho_idmva_cut 			|| Photon_mvaID().at(gHidx[1]) < pho_idmva_cut ) continue;
			if ( Photon_electronVeto().at(gHidx[0]) < pho_eveto_cut 	|| Photon_electronVeto().at(gHidx[1]) < pho_eveto_cut ) continue;
			if ( fabs(Photon_eta().at(gHidx[0])) > pho_eta_cut 		|| fabs(Photon_eta().at(gHidx[1])) > pho_eta_cut ) continue;
			if ( fabs(Photon_eta().at(gHidx[0])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx[0])) < trans_eta_high ) continue;
			if ( fabs(Photon_eta().at(gHidx[1])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx[1])) < trans_eta_high ) continue;

			//trigger requirements
			//if ( year == 2016 && !HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() ) continue;
			//if ( (year == 2017 || year == 2018 ) && !HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() ) continue;

			vector<int> sel_eles;
			for(unsigned int i=0; i<nElectron(); i++){
				if (Electron_pt().at(i) > ele_pt && fabs(Electron_eta().at(i)) < ele_eta && ( fabs(Electron_eta().at(i)) < trans_eta_low || fabs(Electron_eta().at(i)) > trans_eta_high ) && fabs(Electron_dxy().at(i)) < ele_dxy && fabs(Electron_dz().at(i)) < ele_dz 
				&& ( (Electron_pfRelIso03_all().at(i) < ele_pfRelIso && Electron_mvaFall17V2noIso_WP90().at(i)) || (Electron_mvaFall17V2Iso_WP90().at(i) ) ) 
				&& deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx[0]) ) > ele_dR_pho && deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx[1]) ) > ele_dR_pho 
				){
					sel_eles.push_back(i);
				}
			}
			
			vector<int> sel_muons;
			for(unsigned int i=0; i<nMuon(); i++){
				if (Muon_pt().at(i) > muon_pt && fabs(Muon_eta().at(i)) < muon_eta && fabs(Muon_dxy().at(i)) < muon_dxy && fabs(Muon_dz().at(i)) < muon_dz 
				&& Muon_pfRelIso03_all().at(i) < muon_pfRelIso 
				&& Muon_isGlobal().at(i) 
				&& deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx[0]) ) > muon_dR_pho && deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx[1]) ) > muon_dR_pho 
				){
					sel_muons.push_back(i);
				}
			}
			
			vector<int> sel_taus;
			for(unsigned int i=0; i<nTau(); i++){
				if (Tau_pt().at(i) > tau_pt && fabs(Tau_eta().at(i)) < tau_eta && /*Tau_idDecayModeNewDMs().at(i) &&*/ fabs(Tau_dz().at(i)) < tau_dz 
				&& Tau_idDeepTau2017v2p1VSe().at(i) >= tau_deepID_e && Tau_idDeepTau2017v2p1VSmu().at(i) >= tau_deepID_m && Tau_idDeepTau2017v2p1VSjet().at(i) >= tau_deepID_j 
				&& deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx[0]) ) > tau_dR_pho && deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx[1]) ) > tau_dR_pho 
				){

					bool overlap = false;
					for (unsigned int j=0; j<sel_eles.size(); j++){
						if ( deltaR( Tau_p4().at(i) , Electron_p4().at(sel_eles.at(j)) ) < tau_dR_lep ) overlap = true;
					}
					for (unsigned int j=0; j<sel_muons.size(); j++){
						if ( deltaR( Tau_p4().at(i) , Muon_p4().at(sel_muons.at(j)) ) < tau_dR_lep ) overlap = true;
					}

					if ( !overlap ) sel_taus.push_back(i);
				}
			}

			//Isolated Tracks selection
            vector<int> sel_isoTracks;
            for (unsigned int i=0; i<nIsoTrack(); i++){
                if ( IsoTrack_isPFcand().at(i) && IsoTrack_fromPV().at(i) ){
                    LorentzVector *iso_track = new LorentzVector;
                    iso_track->SetXYZT( IsoTrack_pt().at(i)* TMath::Cos(IsoTrack_phi().at(i)) , IsoTrack_pt().at(i)*TMath::Sin( IsoTrack_phi().at(i)), IsoTrack_pt().at(i)*TMath::SinH( IsoTrack_eta().at(i)),  IsoTrack_pt().at(i)*TMath::CosH( IsoTrack_eta().at(i) ) );
                    if ( deltaR( iso_track , Photon_p4().at(gHidx[0]) ) > isoTrk_dR  && deltaR( iso_track , Photon_p4().at(gHidx[1]) ) > isoTrk_dR  ){
						bool iso = true;
						for (unsigned int j=0; j<sel_eles.size(); j++){
                    		if ( deltaR( iso_track , Electron_p4().at(sel_eles.at(j)) ) < isoTrk_dR ){
								iso = false;
                    		}
                		}
						for (unsigned int j=0; j<sel_muons.size(); j++){
                    		if ( deltaR( iso_track , Muon_p4().at(sel_muons.at(j)) ) < isoTrk_dR ){
								iso = false;
                    		}
                		}
						for (unsigned int j=0; j<sel_taus.size(); j++){
                    		if ( deltaR( iso_track , Tau_p4().at(sel_taus.at(j)) ) < isoTrk_dR ){
								iso = false;
                    		}
                		}
						if ( iso ) sel_isoTracks.push_back( i );
					}
				}
			}
			//Jet Selection			
			vector<int> sel_jets;
//			for(unsigned int i=0; i<nJet(); i++){
//				if (Jet_pt().at(i) > jet_pt && fabs(Jet_eta().at(i)) < jet_eta && Jet_neEmEF().at(i) < jet_neEmEF && Jet_neHEF().at(i) < jet_neHEF && Jet_chHEF()[i] > jet_chHEF && Jet_chEmEF()[i] < jet_chEmEF && Jet_nConstituents()[i] > jet_nConstituents && deltaR( Jet_p4().at(i) , Photon_p4().at(gHidx[0]) ) > jet_dR_pho && deltaR( Jet_p4().at(i) , Photon_p4().at(gHidx[1]) ) > jet_dR_pho ){
//
//					bool overlap = false;
//					for (unsigned int j=0; j<sel_eles.size(); j++){
//						if ( !overlap && deltaR( Jet_p4().at(i) , Electron_p4().at(sel_eles.at(j)) ) < jet_dR_lep ){
//							overlap = true;
//							break;
//						}
//					}
//					for (unsigned int j=0; j<sel_muons.size(); j++){
//						if ( !overlap && deltaR( Jet_p4().at(i) , Muon_p4().at(sel_muons.at(j)) ) < jet_dR_lep ){
//							overlap = true;
//							break;
//						}
//					}
//					for (unsigned int j=0; j<sel_taus.size(); j++){
//						if ( !overlap && deltaR( Jet_p4().at(i) , Tau_p4().at(sel_taus.at(j)) ) < jet_dR_tau ){
//							overlap = true;
//							break;
//						}
//					}
//
//					if ( !overlap ){
//						sel_jets.push_back(i);
//						if ( Jet_btagDeepFlavB().at(i) > max_bTag ) max_bTag = Jet_btagDeepFlavB().at(i);
//					}
//				}
//			}
//
//			//bJet Selection			
			vector<int> sel_bJets;
			for(unsigned int i=0; i<sel_jets.size(); i++){
				if ( Jet_btagDeepFlavB().at(sel_jets[i]) > bTag_medium_WP[year - 2016] ) sel_bJets.push_back(sel_jets[i]);
			}


			//Z veto cut
			bool Z_cand = false;
			if ( sel_eles.size() >= 2 ){
				for (unsigned int i=0; i<sel_eles.size(); i++){
					for (unsigned int j=i+1; j<sel_eles.size(); j++){
						if ( (Electron_p4().at(sel_eles[i]) + Electron_p4().at(sel_eles[j])).M() > mZ_veto_low  && (Electron_p4().at(sel_eles[i]) + Electron_p4().at(sel_eles[j])).M() < mZ_veto_up 
							&& Electron_charge().at(sel_eles[i]) * Electron_charge().at(sel_eles[j]) < 0 ){
							Z_cand = true;
							break;
						}
					}
				}
			}
			if ( sel_muons.size() >= 2 ){
				for (unsigned int i=0; i<sel_muons.size(); i++){
					for (unsigned int j=i+1; j<sel_muons.size(); j++){
						if ( (Muon_p4().at(sel_muons[i]) + Muon_p4().at(sel_muons[j])).M() > mZ_veto_low  && (Muon_p4().at(sel_muons[i]) + Muon_p4().at(sel_muons[j])).M() < mZ_veto_up 
							&& Muon_charge().at(sel_muons[i]) * Muon_charge().at(sel_muons[j]) < 0 ){
							Z_cand = true;
							break;
						}
					}
				}
			}
			if ( Z_cand ) continue;

			n_electrons		= sel_eles.size();
			n_muons			= sel_muons.size();
			n_taus			= sel_taus.size();
			n_isoTrks		= sel_isoTracks.size();
			n_jets			= sel_jets.size();
			n_bjets			= sel_bJets.size();

			vector<int> h_cand1, h_cand2;
			vector<vector<int>> raw_results = categorise( sel_eles, sel_muons, sel_taus, sel_isoTracks );
			sel_eles		= raw_results[0];
			sel_muons		= raw_results[1];
			sel_taus		= raw_results[2];
			sel_isoTracks	= raw_results[3];
			h_cand1			= raw_results[4];
			h_cand2			= raw_results[5];

			if ( h_cand1[1] == -1 && h_cand2[1] == -1  ) continue;

			if ( h_cand1[1] == 2 && h_cand2[1] == 1  ) cat1 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 0  ) cat2 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 2  ) cat3 = true;
			if ( h_cand1[1] == 1 && h_cand2[1] == 1  ) cat4 = true;
			if ( h_cand1[1] == 0 && h_cand2[1] == 0  ) cat5 = true;
			if ( h_cand1[1] == 1 && h_cand2[1] == 0  ) cat6 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 3  ) cat7 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == -1 ) cat8 = true;

//			category = 99;
//			if ( cat1 ) category = 1;
//			if ( cat2 ) category = 2;
//			if ( cat3 ) category = 3;
//			if ( cat4 ) category = 4;
//			if ( cat5 ) category = 5;
//			if ( cat6 ) category = 6;
//			if ( cat7 ) category = 7;
//			if ( cat8 ) category = 8;
//
//			vector<classic_svFit::LorentzVector> svFit_res;
//			classic_svFit::LorentzVector diTau_p4, tau1_p4, tau2_p4;
//			float METx	= MET_pt() * TMath::Cos(MET_phi());
//			float METy	= MET_pt() * TMath::Sin(MET_phi());
//
//			if ( category == 1 ){
//				svFit_res = SVfit_all_p4( METx, METy, MET_covXX() , MET_covXY(), MET_covYY(), -1 , Tau_decayMode()[h_cand1[0]], 1 , 3, Muon_pt()[h_cand2[0]], Muon_eta()[h_cand2[0]], Muon_phi()[h_cand2[0]], -1 , Tau_pt()[h_cand1[0]], Tau_eta()[h_cand1[0]], Tau_phi()[h_cand1[0]], Tau_mass()[h_cand1[0]] );
//			}
//			if ( category == 2 ){
//				svFit_res = SVfit_all_p4( METx, METy, MET_covXX() , MET_covXY(), MET_covYY(), -1 , Tau_decayMode()[h_cand1[0]], 2 , 3, Electron_pt()[h_cand2[0]], Electron_eta()[h_cand2[0]], Electron_phi()[h_cand2[0]], -1 , Tau_pt()[h_cand1[0]], Tau_eta()[h_cand1[0]], Tau_phi()[h_cand1[0]], Tau_mass()[h_cand1[0]] );
//			}
//			if ( category == 3 ){
//				svFit_res = SVfit_all_p4( METx, METy, MET_covXX() , MET_covXY(), MET_covYY(), Tau_decayMode()[h_cand1[0]], Tau_decayMode()[h_cand2[0]], 3 , 3, Tau_pt()[h_cand1[0]], Tau_eta()[h_cand1[0]], Tau_phi()[h_cand1[0]], Tau_mass()[h_cand1[0]], Tau_pt()[h_cand2[0]], Tau_eta()[h_cand2[0]], Tau_phi()[h_cand2[0]], Tau_mass()[h_cand2[0]] );
//			}
//			if ( category == 4 ){
//				svFit_res = SVfit_all_p4( METx, METy, MET_covXX() , MET_covXY(), MET_covYY(), -1 , -1 , 1 , 1, Muon_pt()[h_cand1[0]], Muon_eta()[h_cand1[0]], Muon_phi()[h_cand1[0]], -1 ,Muon_pt()[h_cand2[0]], Muon_eta()[h_cand2[0]], Muon_phi()[h_cand2[0]], -1 );
//			}
//			if ( category == 5 ){
//				svFit_res = SVfit_all_p4( METx, METy, MET_covXX() , MET_covXY(), MET_covYY(), -1 , -1 , 2 , 2, Electron_pt()[h_cand1[0]], Electron_eta()[h_cand1[0]], Electron_phi()[h_cand1[0]], -1 , Electron_pt()[h_cand2[0]], Electron_eta()[h_cand2[0]], Electron_phi()[h_cand2[0]], -1 );
//			}
//			if ( category == 6 ){
//				svFit_res = SVfit_all_p4( METx, METy, MET_covXX() , MET_covXY(), MET_covYY(), -1 , -1 , 1 , 2, Muon_pt()[h_cand1[0]], Muon_eta()[h_cand1[0]], Muon_phi()[h_cand1[0]], -1,  Electron_pt()[h_cand2[0]], Electron_eta()[h_cand2[0]], Electron_phi()[h_cand2[0]], -1 );
//			}
//			if ( category == 7 ){
//				int isoTrk_svfit_code = -1;
//				if ( fabs(IsoTrack_pdgId()[h_cand2[1]]) == 11 ) isoTrk_svfit_code = 2;
//				if ( fabs(IsoTrack_pdgId()[h_cand2[1]]) == 13 ) isoTrk_svfit_code = 1;
//				if ( fabs(IsoTrack_pdgId()[h_cand2[1]]) != 11 && fabs(IsoTrack_pdgId()[h_cand2[1]]) != 13 ) isoTrk_svfit_code = 3;
//				//assuming, if hadronic IsoTrack, one-prong WITH neutral pions (larger Br) and massless
//				svFit_res = SVfit_all_p4( METx, METy, MET_covXX() , MET_covXY(), MET_covYY(), 1 , Tau_decayMode()[h_cand1[0]] , isoTrk_svfit_code , 3, IsoTrack_pt()[h_cand2[0]], IsoTrack_eta()[h_cand2[0]], IsoTrack_phi()[h_cand2[0]], 0.0, Tau_pt()[h_cand1[0]], Tau_eta()[h_cand1[0]], Tau_phi()[h_cand1[0]], Tau_mass()[h_cand1[0]] );
//			}
//
//
//			if ( category < 8 ){
//				diTau_p4	= svFit_res[0];
//				tau1_p4		= svFit_res[1];
//				tau2_p4		= svFit_res[2];
//
//				tau1_pt_SVFit	= tau1_p4.pt();
//				tau1_eta_SVFit	= tau1_p4.eta();
//				tau1_phi_SVFit	= tau1_p4.phi();
//				tau1_m_SVFit	= tau1_p4.M();
//				tau2_pt_SVFit	= tau2_p4.pt();
//				tau2_eta_SVFit	= tau2_p4.eta();
//				tau2_phi_SVFit	= tau2_p4.phi();
//				tau2_m_SVFit	= tau2_p4.M();
//			}
//
			float weight = 1.;
			if ( proc != "Data" ) weight = genWeight() * scale_factor;
//
//			t_run			= run();
//			t_lumiBlock		= luminosityBlock();
//			t_event			= event();
//			t_MET_pt		= MET_pt();
//			t_MET_phi		= MET_phi();
			t_weight		= weight;
//			MET_gg_dPhi		= deltaPhi( MET_phi() , (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).phi() );
//
//			gg_pt			=	(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).pt() ;
//			gg_ptmgg		=	gg_pt / mgg;
//			gg_eta			=	(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).eta() ;
//			//gg_eta_bdt		=	gg_eta * sgn( gg_eta ) ;
//			gg_phi			=	(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).phi() ;
//			gg_dR			=	deltaR(Photon_p4().at(gHidx[0]) , Photon_p4().at(gHidx[1])) ;
//			gg_dPhi			=	deltaPhi(Photon_p4().at(gHidx[0]) , Photon_p4().at(gHidx[1])) ;
//			gg_hel_phys		= 	fabs(helicityCosTheta_phys( Photon_p4().at(gHidx[0]), Photon_p4().at(gHidx[1]) ) ); 
//			bool roll		= 	rand() % 2 == 0;
//			if (roll ) 		gg_hel			= 	fabs( helicityCosTheta( Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]), Photon_p4().at(gHidx[0]) ) ) ;
//			else	{		gg_hel			= 	fabs( helicityCosTheta( Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]), Photon_p4().at(gHidx[1]) ) ) ; }
//
//			g1_ptmgg		=	Photon_pt().at(gHidx[0]) / mgg;
//			g1_pt			=	Photon_pt().at(gHidx[0]) ;
//			g1_eta			=	Photon_eta().at(gHidx[0]) ;
//			g1_eta_bdt		=	g1_eta * sgn( gg_eta ) ;
//			g1_phi			=	Photon_phi().at(gHidx[0]) ;
//			g1_idmva		=	Photon_mvaID().at(gHidx[0]) ;
//			g1_pixVeto		=   Photon_pixelSeed().at(gHidx[0]) ;
//			g1_energyErr	=   Photon_energyErr().at(gHidx[0]) ;
//
//			g2_ptmgg		=	Photon_pt().at(gHidx[1]) / mgg;
//			g2_pt			=	Photon_pt().at(gHidx[1]) ;
//			g2_eta			=	Photon_eta().at(gHidx[1]) ;
//			g2_eta_bdt		=	g2_eta * sgn( gg_eta ) ;
//			g2_phi			=	Photon_phi().at(gHidx[1]) ;
//			g2_idmva		=	Photon_mvaID().at(gHidx[1]) ;
//			g2_pixVeto		=   Photon_pixelSeed().at(gHidx[1]) ;
//			g2_energyErr	=   Photon_energyErr().at(gHidx[1]) ;
//
//			if ( g1_pt > g2_pt ){
//				max_g_ptmgg = g1_pt / mgg;
//				min_g_ptmgg = g2_pt / mgg;
//			}
//			else{
//				max_g_ptmgg = g2_pt / mgg;
//				min_g_ptmgg = g1_pt / mgg;
//			}
//			if ( g1_idmva > g2_idmva ){
//				max_g_idmva = g1_idmva;
//				min_g_idmva = g2_idmva;
//			}
//			else{
//				max_g_idmva = g2_idmva;
//				min_g_idmva = g1_idmva;
//			}
//
//			LorentzVector lep1_p4, lep2_p4;
//
//			if ( cat1 ){
//				lep1_p4			=	Muon_p4()[h_cand2[0]];
//				lep1_pt			=	Muon_pt()[h_cand2[0]];
//				lep1_eta		=	Muon_eta()[h_cand2[0]];
//				lep1_eta_bdt	=	Muon_eta()[h_cand2[0]] * sgn( gg_eta );
//				lep1_phi		=	Muon_phi()[h_cand2[0]];	
//				lep1_charge		=	Muon_charge()[h_cand2[0]];
//				lep1_pdgID		=	Muon_pdgId()[h_cand2[0]];
//				lep1_tightID	=	Muon_tightId()[h_cand2[0]];
//
//				lep2_p4			=	Tau_p4()[h_cand1[0]];
//				lep2_pt			=	Tau_pt()[h_cand1[0]];
//				lep2_eta		=	Tau_eta()[h_cand1[0]];
//				lep2_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep2_phi		=	Tau_phi()[h_cand1[0]];	
//				lep2_charge		=	Tau_charge()[h_cand1[0]];
//				lep2_pdgID		=	15 * sgn( lep2_charge );
//				lep2_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
//				lep2_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
//				lep2_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];
//			
//				lep12_dr		= deltaR( Muon_p4()[h_cand2[0]], Tau_p4()[h_cand1[0]] );
//			}
//			if ( cat2 ){
//				lep1_p4			=	Electron_p4()[h_cand2[0]];
//				lep1_pt			=	Electron_pt()[h_cand2[0]];
//				lep1_eta		=	Electron_eta()[h_cand2[0]];
//				lep1_eta_bdt	=	Electron_eta()[h_cand2[0]] * sgn( gg_eta );
//				lep1_phi		=	Electron_phi()[h_cand2[0]];	
//				lep1_charge		=	Electron_charge()[h_cand2[0]];
//				lep1_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand2[0]];
//				lep1_pdgID		=	Electron_pdgId()[h_cand2[0]];
//
//				lep2_p4			=	Tau_p4()[h_cand1[0]];
//				lep2_pt			=	Tau_pt()[h_cand1[0]];
//				lep2_eta		=	Tau_eta()[h_cand1[0]];
//				lep2_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep2_phi		=	Tau_phi()[h_cand1[0]];	
//				lep2_charge		=	Tau_charge()[h_cand1[0]];
//				lep2_pdgID		=	15 * sgn( lep2_charge );
//				lep2_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
//				lep2_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
//				lep2_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];
//
//				lep12_dr		= deltaR( Electron_p4()[h_cand2[0]], Tau_p4()[h_cand1[0]] );
//			}
//			if ( cat3 ){
//				lep1_p4			=	Tau_p4()[h_cand1[0]];
//				lep1_pt			=	Tau_pt()[h_cand1[0]];
//				lep1_eta		=	Tau_eta()[h_cand1[0]];
//				lep1_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep1_phi		=	Tau_phi()[h_cand1[0]];	
//				lep1_charge		=	Tau_charge()[h_cand1[0]];
//				lep1_pdgID		=	15 * sgn( lep1_charge );
//				lep1_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
//				lep1_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
//				lep1_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];
//
//				lep2_p4			=	Tau_p4()[h_cand2[0]];
//				lep2_pt			=	Tau_pt()[h_cand2[0]];
//				lep2_eta		=	Tau_eta()[h_cand2[0]];
//				lep2_eta_bdt	=	Tau_eta()[h_cand2[0]] * sgn( gg_eta );
//				lep2_phi		=	Tau_phi()[h_cand2[0]];	
//				lep2_charge		=	Tau_charge()[h_cand2[0]];
//				lep2_pdgID		=	15 * sgn( lep2_charge );
//				lep2_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand2[0]];	
//				lep2_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand2[0]];	
//				lep2_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand2[0]];
//
//				lep12_dr		= deltaR( Tau_p4()[h_cand2[0]], Tau_p4()[h_cand1[0]] );
//			}
//			if ( cat4 ){
//				lep1_p4			=	Muon_p4()[h_cand1[0]];
//				lep1_pt			=	Muon_pt()[h_cand1[0]];
//				lep1_eta		=	Muon_eta()[h_cand1[0]];
//				lep1_eta_bdt	=	Muon_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep1_phi		=	Muon_phi()[h_cand1[0]];	
//				lep1_charge		=	Muon_charge()[h_cand1[0]];
//				lep1_pdgID		=	Muon_pdgId()[h_cand1[0]];
//				lep1_tightID	=	Muon_tightId()[h_cand1[0]];
//
//				lep2_p4			=	Muon_p4()[h_cand2[0]];
//				lep2_pt			=	Muon_pt()[h_cand2[0]];
//				lep2_eta		=	Muon_eta()[h_cand2[0]];
//				lep2_eta_bdt	=	Muon_eta()[h_cand2[0]] * sgn( gg_eta );
//				lep2_phi		=	Muon_phi()[h_cand2[0]];	
//				lep2_charge		=	Muon_charge()[h_cand2[0]];
//				lep2_pdgID		=	Muon_pdgId()[h_cand2[0]];
//				lep2_tightID	=	Muon_tightId()[h_cand2[0]];
//
//				lep12_dr		= deltaR( Muon_p4()[h_cand2[0]], Muon_p4()[h_cand1[0]] );
//				m_Z				= ( lep1_p4 + lep2_p4 ).M();
//			}
//			if ( cat5 ){
//				lep1_p4			=	Electron_p4()[h_cand1[0]];
//				lep1_pt			=	Electron_pt()[h_cand1[0]];
//				lep1_eta		=	Electron_eta()[h_cand1[0]];
//				lep1_eta_bdt	=	Electron_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep1_phi		=	Electron_phi()[h_cand1[0]];	
//				lep1_charge		=	Electron_charge()[h_cand1[0]];
//				lep1_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand1[0]];
//				lep1_pdgID		=	Electron_pdgId()[h_cand1[0]];
//
//				lep2_p4			=	Electron_p4()[h_cand2[0]];
//				lep2_pt			=	Electron_pt()[h_cand2[0]];
//				lep2_eta		=	Electron_eta()[h_cand2[0]];
//				lep2_eta_bdt	=	Electron_eta()[h_cand2[0]] * sgn( gg_eta );
//				lep2_phi		=	Electron_phi()[h_cand2[0]];	
//				lep2_charge		=	Electron_charge()[h_cand2[0]];
//				lep2_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand2[0]];
//				lep2_pdgID		=	Electron_pdgId()[h_cand2[0]];
//
//				lep12_dr		= deltaR( Electron_p4()[h_cand2[0]], Electron_p4()[h_cand1[0]] );
//				m_Z				= ( lep1_p4 + lep2_p4 ).M();
//			}
//			if ( cat6 ){
//				lep1_p4			=	Muon_p4()[h_cand1[0]];
//				lep1_pt			=	Muon_pt()[h_cand1[0]];
//				lep1_eta		=	Muon_eta()[h_cand1[0]];
//				lep1_eta_bdt	=	Muon_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep1_phi		=	Muon_phi()[h_cand1[0]];	
//				lep1_charge		=	Muon_charge()[h_cand1[0]];
//				lep1_pdgID		=	Muon_pdgId()[h_cand1[0]];
//				lep1_tightID	=	Muon_tightId()[h_cand1[0]];
//
//				lep2_p4			=	Electron_p4()[h_cand2[0]];
//				lep2_pt			=	Electron_pt()[h_cand2[0]];
//				lep2_eta		=	Electron_eta()[h_cand2[0]];
//				lep2_eta_bdt	=	Electron_eta()[h_cand2[0]] * sgn( gg_eta );
//				lep2_phi		=	Electron_phi()[h_cand2[0]];	
//				lep2_charge		=	Electron_charge()[h_cand2[0]];
//				lep2_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand2[0]];
//				lep2_pdgID		=	Electron_pdgId()[h_cand2[0]];
//
//				lep12_dr		= deltaR( Electron_p4()[h_cand2[0]], Muon_p4()[h_cand1[0]] );
//			}
//			if ( cat7 ){
//				lep1_p4					=	Tau_p4()[h_cand1[0]];
//				lep1_pt					=	Tau_pt()[h_cand1[0]];
//				lep1_eta				=	Tau_eta()[h_cand1[0]];
//				lep1_eta_bdt			=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep1_phi				=	Tau_phi()[h_cand1[0]];	
//				lep1_charge				=	Tau_charge()[h_cand1[0]];
//				lep1_pdgID				=	15 * sgn( lep2_charge );
//				lep1_id_vs_e			=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
//				lep1_id_vs_m			=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
//				lep1_id_vs_jet			=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];
//
//				lep2_pt					=	IsoTrack_pt()[h_cand2[0]];
//				lep2_eta				=	IsoTrack_eta()[h_cand2[0]];
//				lep2_eta_bdt			=	IsoTrack_eta()[h_cand2[0]] * sgn( gg_eta );
//				lep2_phi				=	IsoTrack_phi()[h_cand2[0]];	
//				lep2_charge				=	IsoTrack_pdgId()[h_cand2[0]]/fabs(IsoTrack_pdgId()[h_cand2[0]]);
//				lep2_pdgID				=	IsoTrack_pdgId()[h_cand2[0]];
//				lep2_pfRelIso03_all 	= 	IsoTrack_pfRelIso03_all()[h_cand2[0]];;
//				lep2_pfRelIso03_chg 	= 	IsoTrack_pfRelIso03_chg()[h_cand2[0]];;
//
//				LorentzVector iso_track(IsoTrack_pt()[h_cand2[0]], IsoTrack_eta()[h_cand2[0]], IsoTrack_phi()[h_cand2[0]], 0);
//				lep2_p4 = iso_track;
//				lep12_dr		= deltaR( iso_track, Tau_p4()[h_cand1[0]] );
//			}
//
//			if ( cat8 ){
//				lep1_pt			=	Tau_pt()[h_cand1[0]];
//				lep1_eta		=	Tau_eta()[h_cand1[0]];
//				lep1_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
//				lep1_phi		=	Tau_phi()[h_cand1[0]];	
//				lep1_charge		=	Tau_charge()[h_cand1[0]];
//				lep1_pdgID		=	15 * sgn( lep1_charge );
//				lep1_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
//				lep1_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
//				lep1_id_vs_jet	=	(int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];
//			}
//
//			if ( lep1_pt > lep2_pt ){
//				max_lep_pt = lep1_pt;
//				min_lep_pt = lep2_pt;
//			}
//			else{
//				max_lep_pt = lep2_pt;
//				min_lep_pt = lep1_pt;
//			}
//
//			dPhi_MET_l	= deltaPhi( t_MET_phi, lep1_phi );
//
//			if ( category < 8 ){
//
//				if ( lep1_pt > lep2_pt ) dPhi_MET_l	= deltaPhi( t_MET_phi, lep1_phi );
//				else { 					 dPhi_MET_l	= deltaPhi( t_MET_phi, lep2_phi ); }
//
//				MET_ll_dPhi							= deltaPhi( MET_phi() , diTau_p4.phi() );
//
//				lep12_dphi							= deltaPhi( lep2_phi , lep1_phi );
//				lep12_deta							= fabs(lep2_eta - lep1_eta) ;
//				lep12_deta_bdt						= fabs(lep2_eta - lep1_eta) * sgn( gg_eta ) ;
//
//				m_tautau_vis						= (lep1_p4 + lep2_p4).M()	;
//				pt_tautau_vis						= (lep1_p4 + lep2_p4).pt()	;
//				eta_tautau_vis						= (lep1_p4 + lep2_p4).eta()	;
//				eta_tautau_vis_bdt					= eta_tautau_vis * sgn( gg_eta );
//				phi_tautau_vis						= (lep1_p4 + lep2_p4).phi()	;
//
//				gg_tt_CS							= fabs( getCosThetaStar_CS_old( (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])), diTau_p4 ) );
//
//				pt_tautauSVFitLoose					= diTau_p4.pt();
//				eta_tautauSVFitLoose				= diTau_p4.eta();
//				eta_tautauSVFitLoose_bdt			= diTau_p4.eta() * sgn( gg_eta ) ;
//				phi_tautauSVFitLoose				= diTau_p4.phi();
//				m_tautauSVFitLoose					= diTau_p4.M();
//				dR_tautauSVFitLoose					= deltaR( tau1_p4 , tau2_p4 );
//				dR_ggtautauSVFitLoose				= deltaR( (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])), diTau_p4 );
//				dPhi_tautauSVFitLoose				= deltaPhi( tau1_p4 , tau2_p4 );
//				dPhi_ggtautauSVFitLoose				= deltaPhi( (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])), diTau_p4 );
//
//				tt_hel_phys							= 	fabs( helicityCosTheta_phys( tau1_p4, tau2_p4 ) ) ;
//				gg_tt_hel_phys						= 	fabs( helicityCosTheta_phys( (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])) + diTau_p4, (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]))  ) ) ;
//				if (roll ){
//					tt_hel							= 	fabs( helicityCosTheta( diTau_p4 , tau1_p4 ) ) ;
//					gg_tt_hel						= 	fabs( helicityCosTheta( (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])) + diTau_p4, (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]))  ) ) ;
//				}
//				else{
//					tt_hel							= 	fabs( helicityCosTheta( diTau_p4 , tau2_p4 ) ) ;
//					gg_tt_hel						= 	fabs( helicityCosTheta( (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])) + diTau_p4, diTau_p4  ) ) ;
//				}
//
//				mX	= ( diTau_p4 + Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]) ).M() - ( diTau_p4.M() - mHiggs ) - ( mgg - mHiggs );
//			}
//
//
//			//remove main ZGamma bkg in 0tau2lep 
//			if ( category >3 && category < 7 ){
//				m_llg_lead	= ( lep1_p4 + lep2_p4 + Photon_p4().at(gHidx[0]) ).M();
//				m_llg_subl	= ( lep1_p4 + lep2_p4 + Photon_p4().at(gHidx[1]) ).M();
//			}
//			if ( fabs( m_llg_lead - mZ ) < mllg_window | fabs( m_llg_subl - mZ ) < mllg_window  ) continue;
//
//			if ( sel_jets.size() > 0 ){
//				jet1_pt		=	Jet_pt()[sel_jets[0]];
//				jet1_eta	=	Jet_eta()[sel_jets[0]];
//				jet1_eta_bdt=	jet1_eta* sgn( gg_eta );
//				jet1_phi	=	Jet_phi()[sel_jets[0]];
//				jet1_bTag	=	Jet_btagDeepFlavB()[sel_jets[0]];
//				jet1_id		=	Jet_jetId()[sel_jets[0]];
//			}
//			if ( sel_jets.size() > 1 ){
//				jet2_pt		=	Jet_pt()[sel_jets[1]];
//				jet2_eta	=	Jet_eta()[sel_jets[1]];
//				jet2_eta_bdt=	jet2_eta* sgn( gg_eta );
//				jet2_phi	=	Jet_phi()[sel_jets[1]];
//				jet2_bTag	=	Jet_btagDeepFlavB()[sel_jets[1]];
//				jet2_id		=	Jet_jetId()[sel_jets[1]];
//			}

			//make histograms for yields
			h_mgg->Fill( mgg, weight );
			if ( cat1 || cat2 ) h_mgg_1t1l->Fill( mgg, weight );
			if ( cat3 ) h_mgg_2t0l->Fill( mgg, weight );
			if ( cat4 || cat5 || cat6 ) h_mgg_0t2l->Fill( mgg, weight );
			if ( cat7 ) h_mgg_1t0l_iso->Fill( mgg, weight );
			if ( cat8 ) h_mgg_1t0l->Fill( mgg, weight );

			for (unsigned int j=0; j<sel_isoTracks.size(); j++){
				unsigned int i= sel_isoTracks[j];
				h_iso_pt->Fill( IsoTrack_pt().at(i) , weight );
				h_iso_dxy->Fill( IsoTrack_dxy().at(i) , weight );
				h_iso_dz->Fill( IsoTrack_dz().at(i) , weight );
				h_iso_eta->Fill( IsoTrack_eta().at(i) , weight );
				h_iso_fromPV->Fill( IsoTrack_fromPV().at(i) , weight );
				h_iso_miniPFRelIso_all->Fill( IsoTrack_miniPFRelIso_all().at(i) , weight );
				h_iso_miniPFRelIso_chg->Fill( IsoTrack_miniPFRelIso_chg().at(i) , weight );
				h_iso_pdgID->Fill( IsoTrack_pdgId().at(i) , weight );
				h_iso_pfRelIso03_all->Fill( IsoTrack_pfRelIso03_all().at(i) , weight );
				h_iso_pfRelIso03_chg->Fill( IsoTrack_pfRelIso03_chg().at(i) , weight );
				h_iso_phi->Fill( IsoTrack_phi().at(i) , weight );
		}
		h_iso_nIso->Fill( nIsoTrack() , weight );

			out_tree->Fill();

        } // Event loop
        delete file;
    } // File loop
    bar.finish();

	f1->Write();
	f1->Close();
	return 0;
}
