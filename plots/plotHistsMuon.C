#include <stdlib.h>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TProfile.h"
#include "THStack.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"



void comparePlots(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", TString h2Name="name_of_h2", TString h3Name="name_of_h3", TString h4Name="name_of_h4"){

        TFile* file = TFile::Open(path1);

        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        h1_Copy->SetStats(0);

        TH1D *h2 = (TH1D*)file->Get(h2Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h2Name);
        h2_Copy->SetDirectory(0);
        h2_Copy->SetStats(0);

        TH1D *h3 = (TH1D*)file->Get(h3Name);
        TH1D *h3_Copy = (TH1D*)h3->Clone(h3Name);
        h3_Copy->SetDirectory(0);
        h3_Copy->SetStats(0);

        TH1D *h4 = (TH1D*)file->Get(h4Name);
        TH1D *h4_Copy = (TH1D*)h4->Clone(h4Name);
        h4_Copy->SetDirectory(0);
        h4_Copy->SetStats(0);

        file->Close();

	if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 &&  h3_Copy->GetEntries() != 0 && h4_Copy->GetEntries() != 0 ){
//		gStyle->SetOptStat(11);
	        h1_Copy->SetLineColor(kRed);
	        h1_Copy->SetFillColor(kRed);
	        h1_Copy->SetFillStyle(3003);
	        h2_Copy->SetLineColor(kBlue);
	        h2_Copy->SetLineWidth(3);
//	        h2_Copy->SetFillColor(kBlue);
//	        h2_Copy->SetFillStyle(3002);
	        h3_Copy->SetLineColor(kGreen);
	        h3_Copy->SetLineWidth(3);
//	        h3_Copy->SetFillColor(kGreen);
//	        h3_Copy->SetFillStyle(3002);
	        h4_Copy->SetLineColor(kBlack);
	        h4_Copy->SetFillColor(kBlack);
	        h4_Copy->SetFillStyle(3002);
	
	        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
//		c->SetLogy();
		h1_Copy->GetYaxis()->SetLabelOffset(0.005);
		h1_Copy->GetYaxis()->SetTitle("[Hz]    ");  // Define Y ..
		h1_Copy->GetXaxis()->SetLabelOffset(0.005);
//		h1_Copy->GetXaxis()->SetTitle("time [ns]");  // Define Y ..
	        h1_Copy->Draw("HIST,ERR");
	        h2_Copy->Draw("HIST,ERR,SAME");
	        h3_Copy->Draw("HIST,ERR,SAME");
	        h4_Copy->Draw("HIST,ERR,SAME");
	        TLegend *legend = new TLegend(.75, .80, .95, .95);
	        legend->AddEntry(h1_Copy, "Ch18");
	        legend->AddEntry(h2_Copy, "Ch20");
	        legend->AddEntry(h4_Copy, "Ch28");
	        legend->AddEntry(h3_Copy, "Ch21");
		legend->Draw();
	        c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Comparison_"+h1Name+".png");
	        c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Comparison_"+h1Name+".pdf");
		c->Close();
	}

}



void comparePlots_v2( TString path1, TString path2, TString path3, TString hName, TString plotName ){

	const int nBinsA = 18;
	Double_t edgesA[nBinsA+1] = { 10, 20, 40, 70, 150, 250, 400, 750, 1400, 2500, 4500, 9000, 15000, 30000, 50000, 80000, 150000, 300000, 500000 };
	// move bin edges around by 20%
	Double_t edges[nBinsA+1], edgesUp[nBinsA+1], edgesDown[nBinsA+1];
	for (unsigned int i=0; i < nBinsA+1; i++){
		edges[i]     = 1.0*(edgesA[i]);
		edgesUp[i]   = 1.2*(edgesA[i]);
		edgesDown[i] = 0.8*(edgesA[i]);
	}

        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(hName);
        TH1D *h1_Copy = (TH1D*)h1->Clone(hName);
        h1_Copy->SetDirectory(0);
        h1_Copy->SetStats(0);
	h1_Copy->SetBins(nBinsA,edges);
        file->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2 = (TH1D*)file2->Get(hName);
        TH1D *h2_Copy = (TH1D*)h2->Clone(hName);
        h2_Copy->SetDirectory(0);
        h2_Copy->SetStats(0);
	h2_Copy->SetBins(nBinsA,edgesUp);
        file2->Close();

        TFile* file3 = TFile::Open(path3);
        TH1D *h3 = (TH1D*)file3->Get(hName);
        TH1D *h3_Copy = (TH1D*)h3->Clone(hName);
        h3_Copy->SetDirectory(0);
        h3_Copy->SetStats(0);
	h3_Copy->SetBins(nBinsA,edgesDown);
        file3->Close();

        h1_Copy->SetLineColor(kRed);
        h2_Copy->SetLineColor(kBlue);
        h3_Copy->SetLineColor(kGreen);

        h1_Copy->SetLineWidth(3);
        h2_Copy->SetLineWidth(3);
        h3_Copy->SetLineWidth(3);

        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
        c->SetLogx();
        c->SetLogy();
        c->SetGridy();
	if( h2_Copy->GetMaximum() > h1_Copy->GetMaximum() && h2_Copy->GetMaximum() > h3_Copy->GetMaximum() ) h1_Copy->SetMaximum( 1.1*h2_Copy->GetMaximum() );
	if( h3_Copy->GetMaximum() > h1_Copy->GetMaximum() && h3_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h1_Copy->SetMaximum( 1.1*h3_Copy->GetMaximum() );
        h1_Copy->GetYaxis()->SetLabelOffset(0.005);
        h1_Copy->GetYaxis()->SetTitle("Error on #Delta_{t}  [ns]");  // Define Y ..
        h1_Copy->GetXaxis()->SetLabelOffset(0.005);
        h1_Copy->GetXaxis()->SetTitle("pulse Area [pVs]");  // Define Y ..
        h1_Copy->SetTitle("Error on #Delta_{t} w/ varying bin edges");  // Define Y ..
        h1_Copy->Draw();
        h2_Copy->Draw("SAME");
        h3_Copy->Draw("SAME");
        TLegend *legend = new TLegend(.75, .80, .95, .95);
        legend->AddEntry(h3_Copy, "-20%");
        legend->AddEntry(h1_Copy, "+0%");
        legend->AddEntry(h2_Copy, "+20%");
        legend->Draw();
        c->SaveAs("/home/users/fsetti/public_html/milliQan/TimeDelay/"+plotName+".png");
        c->Close();


}


void plotHist( TString path1, TString h1Name, TString outPath = "/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/SimMCratio_", bool twoD = false,  bool logy = false, bool logx = false ){
	
	gStyle->SetOptStat("eoumr");

        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        h1_Copy->SetStats(0);
	file->Close();
	if ( h1_Copy->GetEntries() == 0 ) return;
	
	TCanvas *c = new TCanvas("c", "canvas", 800, 600);
	c->cd();
	c->SetGridy();
//	c->SetGrid();
//	h1_Copy->GetXaxis()->SetTitle("time [ns]");
//	h1_Copy->GetYaxis()->SetTitle("a.u.");
	if ( logy ) c->SetLogy();
	if ( logx ) c->SetLogx();
	if ( !twoD ) h1_Copy->Draw("HIST");
	if ( twoD ) h1_Copy->Draw("colz");
	c->SaveAs(outPath+h1Name+".png");
	c->Close();

}




void plotHist_v2( TH1D *h  ){

        gStyle->SetOptStat("eoumr");
        h->SetDirectory(0);
        h->SetStats(0);
        if ( h->GetEntries() == 0 ) return;
        TString name = h->GetName();

        TCanvas *c = new TCanvas("c", "canvas", 800, 480);
        c->SetGridy();
        h->GetXaxis()->SetTitle("Run #");
        h->GetYaxis()->SetTitle("threshold [mV]");
        h->SetMarkerStyle(34);
        h->Draw();
        c->SaveAs("/home/users/fsetti/public_html/milliQan/Trigger_Efficiency/ThresholdsStudy/"+name+".png");
        c->Close();
}




void compareHists(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", TString h2Name="name_of_h2" ){

	gStyle->SetOptStat("ou");


        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);

        TH1D *h2 = (TH1D*)file->Get(h2Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h2Name);
        h2_Copy->SetDirectory(0);
//        h2_Copy->SetStats(0);
        file->Close();

        if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 ){
                h1_Copy->SetLineColor(kRed);
                h1_Copy->SetFillColor(kRed);
                h1_Copy->SetFillStyle(3003);
                h2_Copy->SetLineColor(kBlue);
                h2_Copy->SetFillColor(kBlue);
                h2_Copy->SetFillStyle(3002);
                if ( h1_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h2_Copy->SetMaximum(h1_Copy->GetMaximum()*1.1);

                TCanvas *c = new TCanvas("c", "canvas", 900, 700);
//		c->SetLogy();
                h2_Copy->GetYaxis()->SetLabelOffset(0.005);
 //               h2_Copy->GetYaxis()->SetTitle("[Hz]    ");  // Define Y ..
                h2_Copy->GetXaxis()->SetLabelOffset(0.005);
                h2_Copy->Draw("HIST,ERR");
                h1_Copy->Draw("HIST,ERR,SAME");
                TLegend *legend = new TLegend(.1, .6, .38, .9);
                legend->AddEntry(h1_Copy, "slab");
                legend->AddEntry(h2_Copy, "same layer, #mu pulse");
                legend->Draw();
                c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Comparison_"+h1Name+"_vs_"+h2Name+".png");
                c->Close();
        }

}


void ratioHists(TString path1, TString path2, TString h1Name, TString outPath ){


        TFile* file1 = TFile::Open(path1);
        TH1D *h1_Copy = (TH1D*)file1->Get(h1Name);
        TH1D *h1 = (TH1D*)h1_Copy->Clone(h1Name);
        h1->SetDirectory(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2_Copy = (TH1D*)file2->Get(h1Name);
        TH1D *h2 = (TH1D*)h2_Copy->Clone(h1Name);
        h2->SetDirectory(0);
        file2->Close();

        if ( h1->GetEntries() != 0 && h2->GetEntries() != 0 ){

		TCanvas *c = new TCanvas("c", "canvas", 800, 800);
		gStyle->SetOptStat("ou");

		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
		pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
//		pad1->SetLogy();         // Vertical grid
		pad1->SetGridx();         // Vertical grid
		pad1->SetGridy();         // Vertical grid
		pad1->Draw();             // Draw the upper pad: pad1
		pad1->cd();               // pad1 becomes the current pad
		h1->SetStats(0);          // No statistics on upper plot
		h1->SetMinimum(0.);
		if ( h1->GetMaximum() < h2->GetMaximum() ) h1->SetMaximum(2.1*h2->GetMaximum());
		h1->Scale(1.,"width");
		h2->Scale(1.,"width");
		h1->Draw();               // Draw h1
		h2->Draw("HIST,SAME");         // Draw h2 on top of h1
		h1->GetXaxis()->SetTitle("Channel #");
		h1->GetYaxis()->SetTitle("Pulse Prob.");
	        TLegend *legend = new TLegend(.75, .75, .95, .95);
                legend->AddEntry(h1, "Data");
                legend->AddEntry(h2, "MC");
                legend->Draw();

		c->cd();          // Go back to the main canvas before defining pad2
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
		pad2->SetTopMargin(0.1);
		pad2->SetBottomMargin(0.25);
		pad2->SetGridy(); // vertical grid
		pad2->Draw();
		pad2->cd();       // pad2 becomes the current pad

		TH1D *h3 = (TH1D*)h1->Clone("h3");
		h3->SetLineColor(kBlack);
		h3->SetMinimum(0.);  // Define Y ..
		h3->SetMaximum(2); // .. range
//		h3->Sumw2();
		h3->SetStats(0);      // No statistics on lower plot
		h3->Divide(h2);
		h3->Draw("ep");       // Draw the ratio plot

//		h1->SetLineColor(kBlue+1);
//		h1->SetLineWidth(3);
		h1->SetMarkerStyle(19);


		h2->SetLineColor(kRed);
		h2->SetLineWidth(3);

		h3->SetTitle(""); // Remove the ratio title


		h3->GetYaxis()->SetTitle("Data/MC   ");
		h3->GetYaxis()->SetNdivisions(505);
		h3->GetYaxis()->SetTitleSize(20);
		h3->GetYaxis()->SetTitleFont(43);
		h3->GetYaxis()->SetTitleOffset(1.55);
		h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		h3->GetYaxis()->SetLabelSize(15);

		h3->GetXaxis()->SetTitle("Channel #");
		h3->GetXaxis()->SetTitleSize(20);
		h3->GetXaxis()->SetTitleFont(43);
		h3->GetXaxis()->SetTitleOffset(4.);
		h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		h3->GetXaxis()->SetLabelSize(15);

		c->SaveAs(outPath+".png");

		c->Close();
        }

}




void ratioHists_v2 (TString path1, TString path2, TString h1Name, TString outPath ){


        TFile* file1 = TFile::Open(path1);
        TH1D *h1_Copy = (TH1D*)file1->Get(h1Name);
        TH1D *h1 = (TH1D*)h1_Copy->Clone(h1Name);
        h1->SetDirectory(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2_Copy = (TH1D*)file2->Get(h1Name);
        TH1D *h2 = (TH1D*)h2_Copy->Clone(h1Name);
        h2->SetDirectory(0);
        file2->Close();

        if ( h1->GetEntries() != 0 && h2->GetEntries() != 0 ){

		TCanvas *c = new TCanvas("c", "canvas", 800, 800);
		gStyle->SetOptStat("ou");

		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
		pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
//		pad1->SetLogy();         // Vertical grid
		pad1->SetGridx();         // Vertical grid
		pad1->SetGridy();         // Vertical grid
		pad1->Draw();             // Draw the upper pad: pad1
		pad1->cd();               // pad1 becomes the current pad
		h1->SetStats(0);          // No statistics on upper plot
		h2->SetStats(0);          // No statistics on upper plot
		h1->SetLineColor(kBlue);
		h1->SetLineWidth(3);
		h2->SetLineColor(kRed);
		h2->SetLineWidth(3);
//		h1->SetMinimum(0.11);
		if ( h1->GetMaximum() > h2->GetMaximum() ){
			h1->Draw("HIST");               // Draw h1
			h2->Draw("HIST,SAME");         // Draw h2 on top of h1
		}
		else{
			h2->Draw("HIST");               // Draw h1
			h1->Draw("HIST,SAME");         // Draw h2 on top of h1
		}
		//h1->GetXaxis()->SetTitle("chan #");
    TLegend *legend = new TLegend(.75, .75, .95, .95);
		legend->AddEntry(h1, "ReReco");
		legend->AddEntry(h2, "20UL");
		legend->Draw();

		c->cd();          // Go back to the main canvas before defining pad2
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
		pad2->SetTopMargin(0.1);
		pad2->SetBottomMargin(0.2);
		pad2->SetGridy(); // vertical grid
		pad2->Draw();
		pad2->cd();       // pad2 becomes the current pad

		TH1D *h3 = (TH1D*)h1->Clone("h3");
		h3->SetMinimum(0);  // Define Y ..
		h3->SetMaximum(2); // .. range
//		h3->Sumw2();
		h3->SetStats(0);      // No statistics on lower plot
		h3->Divide(h2);
		h3->Draw("ep");       // Draw the ratio plot

//		h1->SetLineColor(kBlue+1);
//		h1->SetLineWidth(3);
		h1->SetMarkerStyle(19);
		h3->SetLineColor(kBlack);
		h3->SetMarkerStyle(19);

		h2->SetLineColor(kRed);
		h2->SetLineWidth(3);

		h3->SetTitle(""); // Remove the ratio title

		h3->GetYaxis()->SetTitle("ReReco / UL");
		h3->GetYaxis()->SetNdivisions(505);
		h3->GetYaxis()->SetTitleSize(20);
		h3->GetYaxis()->SetTitleFont(43);
		h3->GetYaxis()->SetTitleOffset(1.55);
		h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		h3->GetYaxis()->SetLabelSize(15);

		//h3->GetXaxis()->SetTitle("ReReco / UL");
		h3->GetXaxis()->SetTitleSize(20);
		h3->GetXaxis()->SetTitleFont(43);
		h3->GetXaxis()->SetTitleOffset(4.);
		h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		h3->GetXaxis()->SetLabelSize(15);

		c->SaveAs(outPath+".png");

		c->Close();
        }
}


void compareHists_v2(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString path2="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", bool logy = false ){

//	gStyle->SetOptStat("uo");


        TFile* file1 = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file1->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        h1_Copy->SetStats(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2 = (TH1D*)file2->Get(h1Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h1Name);
        h2_Copy->SetDirectory(0);
        h2_Copy->SetStats(0);
        file2->Close();
//	h1_Copy->Scale(h2_Copy->GetEntries()/h1_Copy->GetEntries());

        if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 ){
                h1_Copy->SetLineColor(kRed);
                h1_Copy->SetLineWidth(3);
                h1_Copy->SetMarkerStyle(19);
//                h1_Copy->SetFillColor(kRed);
//                h1_Copy->SetFillStyle(3003);
                h2_Copy->SetLineColor(kBlue);
                h2_Copy->SetLineWidth(5);
                h2_Copy->SetMarkerStyle(3);
//                h2_Copy->SetFillColor(kBlue);
//                h2_Copy->SetFillStyle(3002);
//		h1_Copy->Scale(1.,"width");
//		h2_Copy->Scale(1.,"width");
//                if ( h1_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h2_Copy->SetMaximum(h1_Copy->GetMaximum()*2.1);

                TCanvas *c = new TCanvas("c", "canvas", 900, 700);
		if ( logy ) c->SetLogy();
		c->SetLogx();
		c->SetGridy();
                h1_Copy->GetYaxis()->SetLabelOffset(0.005);
                h1_Copy->GetYaxis()->SetTitle("RMS");  // Define Y ..
                h1_Copy->GetXaxis()->SetLabelOffset(0.005);
                h1_Copy->GetXaxis()->SetTitle("nPE");  // Define Y ..
                h1_Copy->Draw("HIST");
                h2_Copy->Draw("HIST,SAME");
	        TLegend *legend = new TLegend(.75, .65, .95, .85);
//                TLegend *legend = new TLegend(.1, .8, .38, .95);
                legend->AddEntry(h1_Copy, "Data");
                legend->AddEntry(h2_Copy, "MC");
                legend->Draw();
                c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/021219/Comparison_"+h1Name+".png");
                c->Close();
        }

}


void compareHists_v3(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString path2="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", TString h2Name="name_of_h2", bool logy = false ){

	gStyle->SetOptStat("uo");


        TFile* file1 = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file1->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2 = (TH1D*)file2->Get(h2Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h2Name);
        h2_Copy->SetDirectory(0);
        file2->Close();
//	h1_Copy->Scale(h2_Copy->GetEntries()/h1_Copy->GetEntries());

        if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 ){
                h1_Copy->SetLineColor(kRed);
                h1_Copy->SetLineWidth(3);
//                h1_Copy->SetFillColor(kRed);
//                h1_Copy->SetFillStyle(3003);
                h2_Copy->SetLineColor(kBlue);
                h2_Copy->SetLineWidth(5);
//                h2_Copy->SetFillColor(kBlue);
//                h2_Copy->SetFillStyle(3002);
                if ( h1_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h2_Copy->SetMaximum(h1_Copy->GetMaximum()*1.1);

                TCanvas *c = new TCanvas("c", "canvas", 900, 700);
		if ( logy ) c->SetLogy();
                h2_Copy->GetYaxis()->SetLabelOffset(0.005);
//                h2_Copy->GetYaxis()->SetTitle("PE");  // Define Y ..
                h2_Copy->GetXaxis()->SetLabelOffset(0.005);
                h2_Copy->Draw("HIST");
                h1_Copy->Draw("HIST,SAME");
	        TLegend *legend = new TLegend(.75, .65, .95, .85);
//                TLegend *legend = new TLegend(.1, .8, .38, .95);
                legend->AddEntry(h1_Copy, "New");
                legend->AddEntry(h2_Copy, "Old");
                legend->Draw();
                c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/021219/Comparison_"+h1Name+".png");
                c->Close();
        }

}


void compareProfiles(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString path2="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString hName="name_of_h" ){

	gStyle->SetEndErrorSize(5);

        TFile* file1 = TFile::Open(path1);
        TProfile *h1 = (TProfile*)file1->Get(hName);
        TProfile *h1_Copy = (TProfile*)h1->Clone(hName);
        h1_Copy->SetDirectory(0);
        h1_Copy->SetStats(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TProfile *h2 = (TProfile*)file2->Get(hName);
        TProfile *h2_Copy = (TProfile*)h2->Clone(hName);
        h2_Copy->SetDirectory(0);
        h2_Copy->SetStats(0);
        file2->Close();

        if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 ){
                h1_Copy->SetLineWidth(3);
                h2_Copy->SetLineWidth(3);
                h1_Copy->SetLineColor(kRed);
                h2_Copy->SetLineColor(kBlue);
//                if ( h1_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h2_Copy->SetMaximum(h1_Copy->GetMaximum()*1.1);

                TCanvas *c = new TCanvas("c", "canvas", 900, 700);
		c->SetGrid();
		c->SetLogx();
                h1_Copy->Draw("E1");
                h2_Copy->Draw("E1,SAME");
//                TLegend *legend = new TLegend(.1, .8, .38, .95);
        	TLegend *legend = new TLegend(.75, .80, .95, .95);
                legend->AddEntry(h1_Copy, "Data");
                legend->AddEntry(h2_Copy, "MC");
                legend->Draw();
                c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/021219/Comparison_"+hName+".png");
                c->Close();
        }

}




void stackHists (TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", TString h2Name="name_of_h2" , TString h3Name="name_of_h3" ){


        THStack *hs = new THStack("hist","");

        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        TH1D *h2 = (TH1D*)file->Get(h2Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h2Name);
        h2_Copy->SetDirectory(0);
        h2_Copy->SetStats(0);
        TH1D *h3 = (TH1D*)file->Get(h3Name);
        TH1D *h3_Copy = (TH1D*)h3->Clone(h3Name);
        h3_Copy->SetDirectory(0);
        h3_Copy->SetStats(0);
        file->Close();

        h1_Copy->SetLineColor(kRed);
        h1_Copy->SetFillColor(kRed);
        h1_Copy->SetFillStyle(3003);
        h2_Copy->SetLineColor(kBlue);
        h2_Copy->SetFillColor(kBlue);
        h2_Copy->SetFillStyle(3002);
        h3_Copy->SetLineColor(kGreen);
        h3_Copy->SetFillColor(kGreen);
        h3_Copy->SetFillStyle(3001);

	hs->Add(h1_Copy);
	hs->Add(h2_Copy);
	hs->Add(h3_Copy);

        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
//	c->SetLogy();
//        hs->GetYaxis()->SetLabelOffset(0.005);
//        hs->GetYaxis()->SetTitle("[Hz]    ");  // Define Y ..
//        hs->GetXaxis()->SetLabelOffset(0.005);
        hs->Draw("HIST,ERR");
        TLegend *legend = new TLegend(.75, .80, .95, .95);
        legend->AddEntry(h1_Copy, "Layer1");
        legend->AddEntry(h2_Copy, "Layer2");
        legend->AddEntry(h3_Copy, "Layer3");
        legend->Draw();
        c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Comparison_Stack_"+h1Name+".png");
        c->Close();


}
