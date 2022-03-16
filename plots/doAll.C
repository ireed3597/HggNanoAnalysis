{
	gROOT->ProcessLine(".L plotHistsMuon.C+");

	TString dir= "/home/users/fsetti/HHggTauTau/HggNanoAnalysis/outputs_UL/";
	vector<TString> years = { "2016", "2017", "2018" };
	vector<TString> procs = { "DiPhoton", "HH_ggTauTau" };
	vector<TString> hists = {"iso_pt","iso_dxy","iso_dz", "iso_eta", "iso_fromPV", "iso_miniPFRelIso_all", "iso_miniPFRelIso_chg", "iso_pdgID", "iso_pfRelIso03_all", "iso_pfRelIso03_chg", "iso_phi", "iso_nIso" };
	TString out_path, reReco_file, UL_file;
	for ( unsigned int i=0; i<years.size(); i++ ){
		TString year = years[i];
		for ( unsigned int j=0; j<procs.size(); j++){
			TString proc = procs[j];
			for ( unsigned int k=0; k<hists.size(); k++){
				TString hist = hists[k];
				out_path = "/home/users/fsetti/public_html/HH2ggtautau/isoTrack_variables/" + year + "/" + proc + "/" + hist ;
				reReco_file	= dir + proc + "_old_" + year + ".root";
				UL_file			= dir + proc + "_new_" + year + ".root";
				ratioHists_v2( reReco_file, UL_file, hist, out_path );
			}
		}
	} 
}
