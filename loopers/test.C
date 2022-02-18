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

int ScanChain( TChain *ch, string proc, int year, float scale_factor = 1, bool resonant = false ) {

	TString file_name = proc + "_check" +  std::to_string(year);

	TFile* f1 = new TFile("test/" + file_name + ".root", "RECREATE");
	H1(mgg, 60, 100 , 180 );
	H1(mgg_1t0l, 60, 100 , 180 );
	H1(mgg_1t1l, 60, 100 , 180 );
	H1(mgg_2t0l, 60, 100 , 180 );
	H1(mgg_0t2l, 60, 100 , 180 );
	H1(mgg_1t0l_iso, 60, 100 , 180 );

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

	out_tree->Branch("m_Z"					,	&m_Z					, 	"m_Z/F"				);	  

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
	out_tree->Branch("lep1_id_vs_e"			,	&lep1_id_vs_e			, 	"lep1_id_vs_e/I"			);	  
	out_tree->Branch("lep1_id_vs_m"			,	&lep1_id_vs_m			,  	"lep1_id_vs_m/I"			);
	out_tree->Branch("lep1_id_vs_jet"		,	&lep1_id_vs_jet			,	"lep1_id_vs_jet/I"			);
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

	//define process ids
	if (proc.find(std::string("HH_ggWW_semileptonic")) != std::string::npos)		process_id = -4;
	if (proc.find(std::string("HH_ggWW_dileptonic")) != std::string::npos) 			process_id = -3;
	//if (proc.find(std::string("HH_ggZZ")) != std::string::npos) 				process_id = -2;
	if (proc.find(std::string("HH_ggZZ_2l2q")) != std::string::npos) 			process_id = -6;
	if (proc.find(std::string("HH_ggZZ_4l")) != std::string::npos) 				process_id = -5;
	if (proc.find(std::string("HH_ggTauTau")) != std::string::npos) 			process_id = -1;
	if (proc.find(std::string("Data")) != std::string::npos) 				process_id = 0;
	if (proc.find(std::string("ZGamma")) != std::string::npos) 				process_id = 2;
	if (proc.find(std::string("DiPhoton")) != std::string::npos) 				process_id = 3;
	if (proc.find(std::string("WGamma")) != std::string::npos) 				process_id = 4;
	if (proc.find(std::string("TTbar")) != std::string::npos) 				process_id = 5;
	if (proc.find(std::string("TTGamma")) != std::string::npos) 				process_id = 6;
	if (proc.find(std::string("TTGG")) != std::string::npos) 				process_id = 7;
	if (proc.find(std::string("GJets")) != std::string::npos)				process_id = 8;
	if (proc.find(std::string("VH")) != std::string::npos) 					process_id = 9;
	if (proc.find(std::string("ttH")) != std::string::npos) 				process_id = 10;
	if (proc.find(std::string("ggH")) != std::string::npos) 				process_id = 11;
	if (proc.find(std::string("VBFH")) != std::string::npos) 				process_id = 12;

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
            tree->LoadTree(loop_event);

            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

			clear_branches();

			//di-photon selection
			mgg = (float)(Photon_p4().at(gHidx()[0]) + Photon_p4().at(gHidx()[1]) ).M();
			if ( mgg < mgg_lower || mgg > mgg_upper ) continue;
			if ( (proc == "Data" || !resonant) && mgg > mgg_sideband_lower && mgg < mgg_sideband_upper ) continue;

			//photon selection
			if ( Photon_pt().at(gHidx()[0]) < pho_pt_cut 				|| Photon_pt().at(gHidx()[1]) < pho_pt_cut ) continue;
			if ( Photon_pt().at(gHidx()[0]) / mgg < lead_pt_mgg_cut 	|| Photon_pt().at(gHidx()[1]) / mgg < sublead_pt_mgg_cut ) continue;
			if ( Photon_mvaID().at(gHidx()[0]) < pho_idmva_cut 			|| Photon_mvaID().at(gHidx()[1]) < pho_idmva_cut ) continue;
			if ( Photon_electronVeto().at(gHidx()[0]) < pho_eveto_cut 	|| Photon_electronVeto().at(gHidx()[1]) < pho_eveto_cut ) continue;
			if ( fabs(Photon_eta().at(gHidx()[0])) > pho_eta_cut 		|| fabs(Photon_eta().at(gHidx()[1])) > pho_eta_cut ) continue;
			if ( fabs(Photon_eta().at(gHidx()[0])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx()[0])) < trans_eta_high ) continue;
			if ( fabs(Photon_eta().at(gHidx()[1])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx()[1])) < trans_eta_high ) continue;

			out_tree->Fill();

        } // Event loop
        delete file;
    } // File loop
    bar.finish();

	f1->Write();
	f1->Close();
	return 0;
}
