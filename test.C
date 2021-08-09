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

#include "NanoCORE/Nano.h"
#include "NanoCORE/Base.h"
#include "NanoCORE/tqdm.h"
#include "NanoCORE/utils.cc"
#include "parameters.h"

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

vector<unsigned int> event_vec;

int ScanChain( TChain *ch, string proc, int year, float scale_factor = 1, bool resonant = false ) {

	TString file_name = proc + "_" +  std::to_string(year);
	TFile* f1 = new TFile("outputs/" + file_name + ".root", "RECREATE");
	H1(mgg, 60, 100 , 180 );
	H1(mgg_1t0l, 60, 100 , 180 );
	H1(mgg_1t1l, 60, 100 , 180 );
	H1(mgg_2t0l, 60, 100 , 180 );
	H1(mgg_0t2l, 60, 100 , 180 );
	H1(lead_pho_pt, 60, 0 , 70 );
	H1(lead_pho_pt_1t0l, 60, 0 , 70 );
	H1(lead_pho_pt_1t1l, 60, 0 , 70 );
	H1(lead_pho_pt_2t0l, 60, 0 , 70 );
	H1(lead_pho_pt_0t2l, 60, 0 , 70 );
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

        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

            nt.GetEntry(event);
            tree->LoadTree(event);

            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

			//di-photon selection
			double mgg = (Photon_p4().at(gHidx()[0]) + Photon_p4().at(gHidx()[1]) ).M();
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

			//trigger requirements
			//if ( year == 2016 && !HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() ) continue;
			//if ( (year == 2017 || year == 2018 ) && !HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() ) continue;

			vector<unsigned int> sel_eles;
			for(unsigned int i=0; i<nElectron(); i++){
				if (Electron_pt().at(i) > ele_pt && fabs(Electron_eta().at(i)) < ele_eta && ( fabs(Electron_eta().at(i)) < trans_eta_low || fabs(Electron_eta().at(i)) > trans_eta_high ) && fabs(Electron_dxy().at(i)) < ele_dxy && fabs(Electron_dz().at(i)) < ele_dz 
				&& ( (Electron_pfRelIso03_all().at(i) < ele_pfRelIso && Electron_mvaFall17V2noIso_WP90().at(i)) || (Electron_mvaFall17V2Iso_WP90().at(i) ) ) 
				&& deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > ele_dR_pho && deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > ele_dR_pho 
				){
					sel_eles.push_back(i);
				}
			}
			
			vector<unsigned int> sel_muons;
			for(unsigned int i=0; i<nMuon(); i++){
				if (Muon_pt().at(i) > muon_pt && fabs(Muon_eta().at(i)) < muon_eta && fabs(Muon_dxy().at(i)) < muon_dxy && fabs(Muon_dz().at(i)) < muon_dz 
				&& Muon_pfRelIso03_all().at(i) < muon_pfRelIso 
				&& Muon_isGlobal().at(i) 
				&& deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > muon_dR_pho && deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > muon_dR_pho 
				){
					sel_muons.push_back(i);
				}
			}
			
			vector<unsigned int> sel_taus;
			for(unsigned int i=0; i<nTau(); i++){
				if (Tau_pt().at(i) > tau_pt && fabs(Tau_eta().at(i)) < tau_eta && Tau_idDecayModeNewDMs().at(i) && fabs(Tau_dz().at(i)) < tau_dz 
				&& Tau_idDeepTau2017v2p1VSe().at(i) >= tau_deepID_e && Tau_idDeepTau2017v2p1VSmu().at(i) >= tau_deepID_m && Tau_idDeepTau2017v2p1VSjet().at(i) >= tau_deepID_j 
				&& deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > tau_dR_pho && deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > tau_dR_pho 
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
			//require opposite sign leptons
			bool os_cut = true;
			//1tau1lep
			if ( sel_taus.size() == 1 && (sel_eles.size() + sel_muons.size() == 1) && sel_eles.size() == 1 && Tau_charge().at(sel_taus.at(0))*Electron_charge().at(sel_eles.at(0)) > 0 ) 	os_cut = false;
			if ( sel_taus.size() == 1 && (sel_eles.size() + sel_muons.size() == 1) && sel_muons.size() == 1 && Tau_charge().at(sel_taus.at(0))*Muon_charge().at(sel_muons.at(0)) > 0 ) 		os_cut = false;
			//2tau0lep
			if ( sel_taus.size() == 2 && (sel_eles.size() + sel_muons.size() == 0) && Tau_charge().at(sel_taus.at(0))*Tau_charge().at(sel_taus.at(1)) > 0 ) os_cut = false;
			//0tau2lep
			if ( sel_taus.size() == 0 && (sel_eles.size() + sel_muons.size() == 2) && sel_eles.size() == 2 && Electron_charge().at(sel_eles.at(0))*Electron_charge().at(sel_eles.at(1)) > 0 ) os_cut = false;
			if ( sel_taus.size() == 0 && (sel_eles.size() + sel_muons.size() == 2) && sel_muons.size() == 2 && Muon_charge().at(sel_muons.at(0))*Muon_charge().at(sel_muons.at(1)) > 0 ) 	os_cut = false;
			if ( sel_taus.size() == 0 && (sel_eles.size() + sel_muons.size() == 2) && sel_eles.size() == sel_muons.size() && Electron_charge().at(sel_eles.at(0))*Muon_charge().at(sel_muons.at(0)) > 0 ) os_cut = false;
			if ( !os_cut ) continue;

			float weight = 1.;
			if ( proc != "Data" ) weight = genWeight() * scale_factor;

			h_mgg->Fill( mgg, weight );
			h_lead_pho_pt->Fill( Photon_pt().at(gHidx().at(0)), weight );

			if ( sel_taus.size() == 1 && (sel_eles.size() + sel_muons.size() == 0 ) ){ 
				h_mgg_1t0l->Fill( mgg, weight );
				h_lead_pho_pt_1t0l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}
			if ( sel_taus.size() == 1 && (sel_eles.size() + sel_muons.size() == 1 ) ){
				h_mgg_1t1l->Fill( mgg, weight );
				h_lead_pho_pt_1t1l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}
			if ( sel_taus.size() == 2 && (sel_eles.size() + sel_muons.size() == 0 ) ){
				h_mgg_2t0l->Fill( mgg, weight );
				h_lead_pho_pt_2t0l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}
			if ( sel_taus.size() == 0 && (sel_eles.size() + sel_muons.size() == 2 ) ){
				h_mgg_0t2l->Fill( mgg, weight );
				h_lead_pho_pt_0t2l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}
        } // Event loop
        delete file;
    } // File loop
    bar.finish();

	f1->Write();
	f1->Close();
	return 0;
}
