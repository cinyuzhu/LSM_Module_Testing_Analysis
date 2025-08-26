#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"


// // C++ utils used by AutoZRange / VectorsToTH2F
#include <limits>
#include <algorithm>
#include "TString.h"  // TString, Form (since you use them)
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include "TInterpreter.h"  // add this include

// // // ... inside MakeComposite(...) BEFORE binding branches:
// // gInterpreter->GenerateDictionary("std::vector<std::vector<int> >, std::vector<std::vector<float> >", "vector");
#include <vector>
#ifdef __CLING__
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#endif


// TH2F* VectorsToTH2F(const std::vector<int>& x, const std::vector<int>& y, const std::vector<float>& val,const char* name = "h2"){
//     if (x.size() != y.size() || y.size() != val.size()) {
//         std::cerr << "VectorsToTH2F: size mismatch x=" << x.size()
//           << " y=" << y.size() << " val=" << val.size() << std::endl;
//         return nullptr;
//     }
//     if (x.empty()) {
//         return new TH2F(name, name, 1, -0.5, 0.5, 1, -0.5, 0.5);
//     }

//     int min_x = std::numeric_limits<int>::max();
//     int max_x = std::numeric_limits<int>::min();
//     int min_y = std::numeric_limits<int>::max();
//     int max_y = std::numeric_limits<int>::min();

//     for (size_t i = 0; i < x.size(); ++i) {
//         min_x = std::min(min_x, x[i]);
//         max_x = std::max(max_x, x[i]);
//         min_y = std::min(min_y, y[i]);
//         max_y = std::max(max_y, y[i]);
//     }

//     const int nx = max_x - min_x + 1;
//     const int ny = max_y - min_y + 1;

//     static int hcounter = 0;
//     TString hname = Form("%s_%d", name, hcounter++);
//     // Match your legacy binning: [min-1, max] with nx bins, and same for Y
//     TH2F* h2 = new TH2F(hname, hname, nx, min_x - 1, max_x, ny, min_y - 1, max_y);

//     for (size_t i = 0; i < x.size(); ++i) {
//         // Direct bin indexing like your old code
//         const int bx = (x[i] - min_x + 1);
//         const int by = (y[i] - min_y + 1);
//         if (bx >= 1 && bx <= nx && by >= 1 && by <= ny) {
//             h2->SetBinContent(bx, by, h2->GetBinContent(bx, by) + val[i]);
//         }
//     }
//     return h2;
// }

TH2F* VectorsToTH2F(const std::vector<int>& x,
    const std::vector<int>& y,
    const std::vector<float>& val,
    const char* name="h2")
{
// ...
int min_x = *std::min_element(x.begin(), x.end());
int max_x = *std::max_element(x.begin(), x.end());
int min_y = *std::min_element(y.begin(), y.end());
int max_y = *std::max_element(y.begin(), y.end());

// Edges at half-integers so each integer pixel is the bin CENTER
TH2F* h2 = new TH2F(name, name,
        max_x - min_x + 1, min_x - 0.5, max_x + 0.5,
        max_y - min_y + 1, min_y - 0.5, max_y + 0.5);

// Fill at the *pixel center* coordinates
for (size_t i = 0; i < x.size(); ++i)
h2->Fill(x[i], y[i], val[i]);

return h2;
}


double compute_skewness(TH1D* h) {
    double mean = h->GetMean();
    double sigma = h->GetRMS();
    if (sigma == 0) return 0.0;

    double skew = 0.0;
    double norm = 0.0;

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double x = h->GetBinCenter(i);
        double y = h->GetBinContent(i);
        skew += y * std::pow((x - mean) / sigma, 3);
        norm += y;
    }

    return (norm > 0) ? skew / norm : 0.0;
}



void MakeComposite_Fe55_Clusters(const char* filelist = "test_list.txt", const int image_number = 4) {
    TChain chain("clustersRec");

    std::ifstream fin(filelist);
    std::string filename;
    while (std::getline(fin, filename)) {
        std::cout << "Adding file: " << filename << std::endl;
        chain.Add(filename.c_str());
    }

    // Set up branch pointers
    std::vector<std::vector<int>>* pixels_x = nullptr;
    std::vector<std::vector<int>>* pixels_y = nullptr;
    std::vector<std::vector<float>>* pixels_E = nullptr;
    std::vector<float>* PosX = nullptr;
    std::vector<float>* PosY = nullptr;
    // std::vector<float>* valid = nullptr;
    std::vector<float>*  cluster_rmsxy= nullptr;
    std::vector<float>*  cluster_rmsx= nullptr;
    std::vector<float>*  cluster_rmsy= nullptr;
    std::vector<float>* cluster_E = nullptr;

    chain.SetBranchAddress("pixels_x", &pixels_x);
    chain.SetBranchAddress("pixels_y", &pixels_y);
    chain.SetBranchAddress("pixels_E", &pixels_E);
    chain.SetBranchAddress("PosX", &PosX);
    chain.SetBranchAddress("PosY", &PosY);
    // chain.SetBranchAddress("PosY", &valid);
    chain.SetBranchAddress("wSTD_XY", &cluster_rmsxy);
    chain.SetBranchAddress("fSTD_X", &cluster_rmsx);
    chain.SetBranchAddress("fSTD_Y", &cluster_rmsy);

    chain.SetBranchAddress("Energy", &cluster_E);

    const int nbins = 401;
    const double range = 4;
    const int bins_per_pix = 25;
    const double efact = 3.8;
    int total_clusters = 0;
    int nAccepted = 0;
    int nAccepted_front = 0;
    int nAccepted_back = 0;
    int nAccepted_En = 0;
    int nColumns = 0;
    int nRows = 0;
    double min_energy = 4.5;
    double max_energy = 7.0;


    if(image_number==4){
        nColumns = 6400;
        nRows = 1600;
    } 
    else if(image_number==5){
        nColumns = 1280;
        nRows = 320;
    } 
    else if(image_number==6){
        nColumns = 640;
        nRows = 1600;
    }
    else if(image_number==7){
        nColumns = 640;
        nRows = 1600;
    }
    
    


    double sigma_xy_front_cut = 1.0;
    Long64_t nEntries = chain.GetEntries();

    TH2D* hcomp_front = new TH2D("hcomp_front", "Composite Cluster Front Events", nbins, -range, range, nbins, -range, range);
    TH2D* hcomp_back = new TH2D("hcomp_back", "Composite Cluster Back Events", nbins, -range, range, nbins, -range, range);
    TH1F* hPosX = new TH1F("hPosX", "PosX", 20, 0, nColumns);
    TH1F* hPosY = new TH1F("hPosY", "PosY", 20, 0, nRows);
    TH2F* hPosXY = new TH2F("hPosXY", "PosXY", (int)(nColumns/100), 0, nColumns, (int)(nRows/100), 0, nRows);
    TH1F* hSigmaxy = new TH1F("hSigmaxy", "Sigmaxy", 20, 0, 2.0);
    TH1F* hSigmax = new TH1F("hSigmax", "Sigmax", 20, 0, 2.0);
    TH1F* hSigmay = new TH1F("hSigmay", "Sigmay", 20, 0, 2.0);
    TH1F* hClusterE = new TH1F("hClusterE", "hClusterE", 50, min_energy, max_energy);
    // TH2F* hposx_sigmaxy = new TH2F("hposx_sigmaxy", "hposx_sigmaxy", 640, 0, 640, 20, 0, 2.0);
    TH2F* hposx_sigmaxy = new TH2F("hposx_sigmaxy", "posx_sigmaxy", nColumns, 0, nColumns, 20, 0, 2.0);
    TH2F* henergy_sigmaxy = new TH2F("hposx_sigmaxy", "posx_sigmaxy", 50, min_energy, max_energy, 20, 0, 2.0);


    

    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        chain.GetEntry(entry);
        int nClusters = PosX->size();
        total_clusters += nClusters;
        std::cout << "Found " << nClusters << " clusters in entry " << entry << std::endl;


        for (size_t n = 0; n < nClusters; n++) {
            double cur_cluster_E = cluster_E->at(n); // Energy of the cluster
            double cur_cluster_x = PosX->at(n); // Mean x position of the cluster
            double cur_cluster_y = PosY->at(n); // Mean y position of the cluster
            // double cur_cluster_rmsxy = cluster_rmsxy->at(n); // Standard deviation in xy for the cluster
            double cur_cluster_rmsx = cluster_rmsx->at(n); // Standard deviation in xy for the cluster
            double cur_cluster_rmsy = cluster_rmsy->at(n); // Standard deviation in xy for the cluster
            double cur_cluster_rmsxy = std::sqrt((cur_cluster_rmsx*cur_cluster_rmsx + cur_cluster_rmsy*cur_cluster_rmsy)/2.0); // Combined standard deviation in xy for the cluster

            const auto& cur_pixel_x = pixels_x->at(n); // Pixel x coordinates for the cluster (vector of ints) 
            const auto& cur_pixel_y = pixels_y->at(n); // Pixel y coordinates for the cluster (vector of ints)
            const auto& cur_pixel_energy = pixels_E->at(n); // Pixel energies for the cluster (vector of floats)
            double borq = 0;

            double cluster_rms_temp_variable = 0; // Temporary variable to store rms value depending of which image we are processing. This defines the front/back cut
            
            
            if (cur_cluster_E > min_energy && cur_cluster_E < max_energy ) {
                // hClusterE->Fill(cur_cluster_E);
                hSigmaxy->Fill(cur_cluster_rmsxy);
                hSigmax->Fill(cur_cluster_rmsx);
                hSigmay->Fill(cur_cluster_rmsy);
                hposx_sigmaxy->Fill(cur_cluster_x, cur_cluster_rmsxy,cur_cluster_E);
                hPosX->Fill(cur_cluster_x);
                hPosY->Fill(cur_cluster_y);
                hPosXY->Fill(cur_cluster_x, cur_cluster_y, cur_cluster_E);
                henergy_sigmaxy->Fill(cur_cluster_E, cur_cluster_rmsxy);
                // hClusterE->Fill(cur_cluster_E);
                nAccepted_En++;
            } 

            if(image_number==4){
                cluster_rms_temp_variable = cur_cluster_rmsxy;
            } 
            else if(image_number==5){
                cluster_rms_temp_variable = cur_cluster_rmsxy;
                
            } 
            else if(image_number==6){
                cluster_rms_temp_variable = cur_cluster_rmsx;

            }
            else if(image_number==7){
                cluster_rms_temp_variable = cur_cluster_rmsy;
            }

            // Selection for front side cluster events
            if (cur_cluster_E > min_energy && cur_cluster_E < max_energy && cluster_rms_temp_variable < sigma_xy_front_cut) {
                hClusterE->Fill(cur_cluster_E);
                TH2F* hview = VectorsToTH2F(cur_pixel_x, cur_pixel_y, cur_pixel_energy, "h2_cluster");
                if (!hview) continue;

                for (int i = 1; i <= nbins; i++) {
                    for (int j = 1; j <= nbins; j++) {
                        const double gx = cur_cluster_x + hcomp_front->GetXaxis()->GetBinLowEdge(i);
                        const double gy = cur_cluster_y + hcomp_front->GetYaxis()->GetBinLowEdge(j);
                        const double bc = hview->GetBinContent(hview->FindBin(gx, gy)) / (bins_per_pix * bins_per_pix);
                        if (std::abs(hcomp_front->GetXaxis()->GetBinCenter(i)) > 3.0 &&
                            std::abs(hcomp_front->GetYaxis()->GetBinCenter(j)) > 3.0) {
                            borq += bc;
                        }
                    }
                }

                // if (std::abs(borq) < 25.0) { // selection on the maximum charge in the corner pixels

                    // Accumulate into composite, centered on (mean_x,mean_y)
                    for (int i = 1; i <= nbins; ++i) {
                        for (int j = 1; j <= nbins; ++j) {
                            const double gx = cur_cluster_x + hcomp_front->GetXaxis()->GetBinLowEdge(i);
                            const double gy = cur_cluster_y + hcomp_front->GetYaxis()->GetBinLowEdge(j);
                            const double bc = hview->GetBinContent(hview->FindBin(gx, gy))/ (bins_per_pix * bins_per_pix);
                            hcomp_front->SetBinContent(i, j, hcomp_front->GetBinContent(i, j) + bc);
                        }
                    }
                    nAccepted_front++;
                // }

                if ((n % 10) == 0) {
                    std::cout << "Processed cluster " << (n + 1) << " / " << nClusters << std::endl;
                }
                hview->Delete();
            } 


            // Selection for back side cluster events
            if (cur_cluster_E > min_energy && cur_cluster_E < max_energy && cluster_rms_temp_variable > sigma_xy_front_cut) {
                TH2F* hview = VectorsToTH2F(cur_pixel_x, cur_pixel_y, cur_pixel_energy, "h2_cluster");
                if (!hview) continue;

                for (int i = 1; i <= nbins; i++) {
                    for (int j = 1; j <= nbins; j++) {
                        const double gx = cur_cluster_x + hcomp_back->GetXaxis()->GetBinLowEdge(i);
                        const double gy = cur_cluster_y + hcomp_back->GetYaxis()->GetBinLowEdge(j);
                        const double bc = hview->GetBinContent(hview->FindBin(gx, gy)) / (bins_per_pix * bins_per_pix);
                        if (std::abs(hcomp_back->GetXaxis()->GetBinCenter(i)) > 3.0 &&
                            std::abs(hcomp_back->GetYaxis()->GetBinCenter(j)) > 3.0) {
                            borq += bc;
                        }
                    }
                }

                if (std::abs(borq) < 25.0) {

                    // Accumulate into composite, centered on (mean_x,mean_y)
                    for (int i = 1; i <= nbins; ++i) {
                        for (int j = 1; j <= nbins; ++j) {
                            const double gx = cur_cluster_x + hcomp_back->GetXaxis()->GetBinLowEdge(i);
                            const double gy = cur_cluster_y + hcomp_back->GetYaxis()->GetBinLowEdge(j);
                            const double bc = hview->GetBinContent(hview->FindBin(gx, gy))/ (bins_per_pix * bins_per_pix);
                            hcomp_back->SetBinContent(i, j, hcomp_back->GetBinContent(i, j) + bc);
                        }
                    }
                    nAccepted_back++;
                }

                if ((n % 10) == 0) {
                    std::cout << "Processed cluster " << (n + 1) << " / " << nClusters << std::endl;
                }
                hview->Delete();

            } 
            
        }
    }

    std::cout << "--------------------------------------"  << std::endl;
    std::cout << "Total clusters processed in "<<nEntries<<" processed images: " << total_clusters << std::endl;
    std::cout << "Accepted clusters passing energy cut: " << nAccepted_En << std::endl;
    std::cout << "Accepted clusters passing front-side event cut: " << nAccepted_front << std::endl;
    std::cout << "Accepted clusters passing back-side event cut: " << nAccepted_back << std::endl;
    hcomp_front->Scale(1. / nAccepted_front);
    hcomp_back->Scale(1. / nAccepted_back);
    std::cout << "Total integral for front-side events [e-]: " << hcomp_front->Integral() * 1000 / 3.8 << std::endl;
    std::cout << "Total integral for back-side events [e-]: " << hcomp_back->Integral() * 1000 / 3.8 << std::endl;
    std::cout << "--------------------------------------"  << std::endl;
    
    // Print statistical variables
    int mid_bin = nbins / 2 + 1;
    std::cout << "Center bin X FRONT: " << hcomp_front->GetXaxis()->GetBinCenter(mid_bin) << std::endl;
    std::cout << "Center bin Y FRONT: " << hcomp_front->GetYaxis()->GetBinCenter(mid_bin) << std::endl;
    std::cout << "Center bin X BACK: " << hcomp_back->GetXaxis()->GetBinCenter(mid_bin) << std::endl;
    std::cout << "Center bin Y BACK: " << hcomp_back->GetYaxis()->GetBinCenter(mid_bin) << std::endl;


    //---------------------------------------------//
    //-------- Analysis for Back-Side events -------//
    //---------------------------------------------//

    gStyle->SetOptStat(0);
    // gStyle->SetNumberContours(255);
    // gStyle->SetPalette(kBird); // nice perceptual palette; switch to kViridis if you have it

    // Canvas with two panels: left = heatmap, right = projections
    TCanvas* cBack = new TCanvas("cBack", "Back-side composite + projections", 1200, 550);
    cBack->Divide(2,1);

    // ---------- LEFT: composite heatmap ----------
    cBack->cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.12);
    gPad->SetLogz();

    hcomp_back->SetTitle("Back-side Composite Cluster;#Delta x [pixels];#Delta y [pixels]");
    hcomp_back->Draw("COLZ");

    // (optional) quick geometry printout
    // {
    // auto ax = hcomp_back->GetXaxis();
    // double xmin = ax->GetXmin(), xmax = ax->GetXmax();
    // double bw   = ax->GetBinWidth(1);
    // int mid     = nbins/2 + 1;
    // std::cout<<std::setprecision(12)
    //         <<"xmin="<<xmin<<" xmax="<<xmax
    //         <<" bw="<<bw<<" mid_center="<<ax->GetBinCenter(mid)<<"\n";
    // }

    // ---------- RIGHT: edge projections (normalized & styled) ----------
    cBack->cd(2);
    gPad->SetGridx(true);
    gPad->SetGridy(true);
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.12);

    // Project charge onto x/y along the first pixel from -4 to -3 (bins 51..100)
    TH1D* hbelow = hcomp_back->ProjectionX("hbelow", 51, 100); hbelow->SetTitle("ProjX (Y: -4 #rightarrow -3 px)");
    TH1D* hleft  = hcomp_back->ProjectionY("hleft",  51, 100); hleft->SetTitle("ProjY (X: -4 #rightarrow -3 px)");

    // Print statistical variables of projections (pre-normalization)
    std::cout << "--------------------------------------\n";
    std::cout << "BACK-SIDE EVENTS CTI METRICS\n";

    const double ePerADU = 3.8;
    const double hleft_int_e  = hleft->Integral()*1000.0/ePerADU;
    const double hbelow_int_e = hbelow->Integral()*1000.0/ePerADU;

    std::cout << "hleft  integral [e-]: " << hleft_int_e  << std::endl;
    std::cout << "Mean RMS Skewness: " << hleft->GetMean() << " "
            << hleft->GetRMS() << " " << hleft->GetSkewness() << std::endl;

    std::cout << "hbelow integral [e-]: " << hbelow_int_e << std::endl;
    std::cout << "Mean RMS Skewness: " << hbelow->GetMean() << " "
            << hbelow->GetRMS() << " " << hbelow->GetSkewness() << std::endl;
    std::cout << "--------------------------------------\n";

    // Style
    hbelow->SetLineColor(kBlack);
    hbelow->SetLineWidth(3);
    hbelow->SetFillColorAlpha(kBlack, 0.08);
    hbelow->SetMarkerStyle(20);
    hbelow->SetMarkerColor(kBlack);

    hleft->SetLineColor(kBlue+1);
    hleft->SetLineWidth(3);
    hleft->SetFillColorAlpha(kBlue+1, 0.08);
    hleft->SetMarkerStyle(20);
    hleft->SetMarkerColor(kBlue+1);

    // Normalize to unit area for a fair overlay
    if (hbelow->Integral() > 0) hbelow->Scale(1.0 / hbelow->Integral());
    if (hleft->Integral()  > 0) hleft->Scale(1.0  / hleft->Integral());

    // Axis cosmetics (use X-projection to define frame)
    hbelow->SetTitle(";Pixel offset [px];Normalized counts");
    hbelow->GetXaxis()->SetTitleSize(0.05);
    hbelow->GetYaxis()->SetTitleSize(0.05);
    hbelow->GetXaxis()->SetLabelSize(0.045);
    hbelow->GetYaxis()->SetLabelSize(0.045);
    hbelow->GetYaxis()->SetNdivisions(505);

    // Common y-range with 10% headroom
    double ymax = std::max(hbelow->GetMaximum(), hleft->GetMaximum());
    hbelow->SetMaximum(ymax * 1.15);
    hbelow->SetMinimum(0.0);

    // Draw overlay
    hbelow->Draw("HIST");        // frame + first curve
    hleft->Draw("HIST SAME");    // second curve

    // Legend
    auto leg = new TLegend(0.58782, 0.721, 0.8869, 0.880);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.025);
    leg->AddEntry(hbelow, "ProjX (Y: -4 #rightarrow -3 px)", "l");
    leg->AddEntry(hleft,  "ProjY (X: -4 #rightarrow -3 px)", "l");
    leg->Draw();

    // On-plot stats box
    TPaveText* pave = new TPaveText(0.15, 0.72, 0.55, 0.88, "NDC");
    pave->SetFillColor(0);
    pave->SetFillStyle(0);
    pave->SetBorderSize(0);
    pave->SetTextFont(42);
    pave->SetTextSize(0.04);
    // pave->AddText(Form("hleft integral [e^-]:  %.1f", hleft_int_e));
    // pave->AddText(Form("hbelow integral [e^-]: %.1f", hbelow_int_e));
    pave->Draw();

    // (optional) save the canvas
    // cBack->SaveAs("back_composite_projections.png");


    // //---------------------------------------------//
    // //--------Analysis for Front Side events-------//
    // //---------------------------------------------//

    // //Get charge for central pixel and pixels around it in composite
    // hcomp_front->SetAxisRange(-0.5,0.5-0.001,"X");
    // hcomp_front->SetAxisRange(-0.5,0.5-0.001,"Y");
    // double centrp = hcomp_front->Integral();
    
    // hcomp_front->SetAxisRange(-0.5,0.5-0.001,"X");
    // hcomp_front->SetAxisRange(0.5,1.5-0.001,"Y");
    // double abovep = hcomp_front->Integral();
    
    // hcomp_front->SetAxisRange(-0.5,0.5-0.001,"X");
    // hcomp_front->SetAxisRange(-1.5,-0.5-0.001,"Y");
    // double belowp = hcomp_front->Integral();
    
    // hcomp_front->SetAxisRange(-0.5,0.5-0.001,"Y");
    // hcomp_front->SetAxisRange(0.5,1.5-0.001,"X");
    // double rightp = hcomp_front->Integral();
    
    // hcomp_front->SetAxisRange(-0.5,0.5-0.001,"Y");
    // hcomp_front->SetAxisRange(-1.5,-0.5-0.001,"X");
    // double leftp = hcomp_front->Integral();
    
    // //Print fractions of pixel contents relative to central pixel

    // std::cout << "--------------------------------------"  << std::endl;
    // std::cout << "FRONT-SIDE EVENTS CTI METRICS"  << std::endl;
    // std::cout << "Above fraction: " << abovep/centrp << std::endl;
    // std::cout << "Below fraction: " << belowp/centrp << std::endl;
    // std::cout << "Right fraction: " << rightp/centrp << std::endl;
    // std::cout << "Left fraction: " << leftp/centrp << std::endl;
    // std::cout << "--------------------------------------"  << std::endl;

    
    // hcomp_front->SetAxisRange(-4,4-0.001,"X");
    // hcomp_front->SetAxisRange(-4,4-0.001,"Y");
    // TCanvas* cview_front = new TCanvas();

    // hcomp_front->Draw("COLZ");
    // cview_front->SetLogz();

    //---------------------------------------------//
    //-------- Analysis for Front-Side events ------//
    //---------------------------------------------//

    // helper: integrate a rectangle without mutating axis ranges
    auto SumRect = [&](TH2D* h, double xlo, double xhi, double ylo, double yhi) -> double {
        auto ax = h->GetXaxis();
        auto ay = h->GetYaxis();
        const double eps = 1e-6;
        int x1 = ax->FindBin(xlo + eps), x2 = ax->FindBin(xhi - eps);
        int y1 = ay->FindBin(ylo + eps), y2 = ay->FindBin(yhi - eps);
        return h->Integral(x1, x2, y1, y2);
    };
    
    // center & neighbors
    const double centrp = SumRect(hcomp_front, -0.5,  0.5, -0.5,  0.5);
    const double abovep = SumRect(hcomp_front, -0.5,  1.5,  0.5,  2.5);
    const double belowp = SumRect(hcomp_front, -0.5,  0.5, -1.5, -0.5);
    const double rightp = SumRect(hcomp_front,  1.5,  -0.5, 2.5,  0.5);
    const double leftp  = SumRect(hcomp_front, -1.5, -0.5, -0.5,  0.5);
    
    const double fAbove = (centrp>0)? abovep/centrp : 0.0;
    const double fBelow = (centrp>0)? belowp/centrp : 0.0;
    const double fRight = (centrp>0)? rightp/centrp : 0.0;
    const double fLeft  = (centrp>0)? leftp /centrp : 0.0;
    
    // clean 2-panel canvas
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(255);
    gStyle->SetPalette(kBird);
    
    TCanvas* cFront = new TCanvas("cFront","Front-side composite + CTI fractions", 1200, 550);
    cFront->Divide(2,1);
    
    // -------- LEFT: composite with color-coded boxes (no text) --------
    cFront->cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.12);
    gPad->SetLogz();
    
    hcomp_front->SetTitle("Front-side Composite;#Delta x [pixels];#Delta y [pixels]");
    hcomp_front->GetXaxis()->SetRangeUser(-4,4);
    hcomp_front->GetYaxis()->SetRangeUser(-4,4);
    hcomp_front->Draw("COLZ");
    
    // only draw boxes, color-coded
    auto drawBoxOnly = [&](double x1,double y1,double x2,double y2, Color_t lc){
        TBox* b = new TBox(x1,y1,x2,y2);
        b->SetLineColor(lc);
        b->SetLineWidth(3);
        b->SetFillStyle(0); // transparent
        b->Draw("same");
        return b; // return pointer for legend prototypes
    };
    
    TObject* boxCenter = drawBoxOnly(-0.5,-0.5, 0.5, 0.5, kBlack);
    TObject* boxAbove  = drawBoxOnly(-0.5, 1.5, 0.5, 2.5, kBlue+1);
    TObject* boxBelow  = drawBoxOnly(-0.5,-1.5, 0.5,-0.5, kRed+1);
    TObject* boxRight  = drawBoxOnly( 1.5,-0.5, 2.5, 0.5, kGreen+2);
    TObject* boxLeft   = drawBoxOnly(-1.5,-0.5,-0.5, 0.5, kMagenta);
    
    // -------- RIGHT: color-coded bars + legend --------
    cFront->cd(2);
    gPad->SetGridx(true);
    gPad->SetGridy(true);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.12);
    
    // frame only (no automatic bars)
    TH1F* hFrame = new TH1F("hFracFrame","Neighbor Fractions;Region;Fraction (w.r.t. center)", 5, 0.5, 5.5);
    hFrame->GetXaxis()->SetBinLabel(1,"Above");
    hFrame->GetXaxis()->SetBinLabel(2,"Below");
    hFrame->GetXaxis()->SetBinLabel(3,"Right");
    hFrame->GetXaxis()->SetBinLabel(4,"Left");
    hFrame->GetXaxis()->SetBinLabel(5,"Center=1.0");
    
    double yMaxFrac = std::max(1.0, std::max({fAbove, fBelow, fRight, fLeft}));
    hFrame->SetMaximum(yMaxFrac * 1.25);
    hFrame->SetMinimum(0.0);
    hFrame->GetXaxis()->SetLabelSize(0.05);
    hFrame->GetYaxis()->SetLabelSize(0.045);
    hFrame->GetXaxis()->SetTitleSize(0.05);
    hFrame->GetYaxis()->SetTitleSize(0.05);
    hFrame->Draw("AXIS");  // axes only
    
    // draw a dashed reference at 1.0
    TLine* ref = new TLine(0.5, 1.0, 5.5, 1.0);
    ref->SetLineStyle(2);
    ref->SetLineColor(kGray+2);
    ref->Draw("same");
    
    // helper to paint one colored bar (using a TBox so each bar can have its own color)
    auto drawBar = [&](int bin, double y, Color_t col)->TBox* {
        double x1 = hFrame->GetXaxis()->GetBinLowEdge(bin);
        double x2 = hFrame->GetXaxis()->GetBinUpEdge(bin);
        TBox* b = new TBox(x1, 0.0, x2, y);
        b->SetFillColorAlpha(col, 0.35);
        b->SetLineColor(col);
        b->SetLineWidth(3);
        b->Draw("same");
        return b;
    };
    
    TBox* bAbove  = drawBar(1, fAbove,  kBlue+1);
    TBox* bBelow  = drawBar(2, fBelow,  kRed+1);
    TBox* bRight  = drawBar(3, fRight,  kGreen+2);
    TBox* bLeft   = drawBar(4, fLeft,   kMagenta);
    TBox* bCenter = drawBar(5, 1.0,     kBlack);
    
    // legend mirrors the box colors; no extra text on the heatmap
    auto legF = new TLegend(0.205, 0.537, 0.504, 0.738);
    legF->SetBorderSize(0);
    legF->SetFillStyle(0);
    legF->AddEntry(bCenter, Form("Center charge = %.0f", centrp), "f");
    legF->AddEntry(bAbove,  Form("Above  = %.3f", fAbove), "f");
    legF->AddEntry(bBelow,  Form("Below  = %.3f", fBelow), "f");
    legF->AddEntry(bRight,  Form("Right  = %.3f", fRight), "f");
    legF->AddEntry(bLeft,   Form("Left   = %.3f", fLeft ), "f");
    legF->Draw();
  

    //------------------------------------------------------------------//
    // Fit the cluster energy distribution with a double-Gaussian model
    //------------------------------------------------------------------//

    TCanvas* c_cluster_E = new TCanvas("c_cluster_E","Cluster Energy (double-Gaussian)", 900, 600);
    c_cluster_E->cd();
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.12);
    // gPad->SetGridx(true);
    // gPad->SetGridy(true);

    // Style like plt.bar(..., color="gray", edgecolor="black", alpha=0.6)
    gPad->SetGridx(); gPad->SetGridy();

    hClusterE->SetTitle("Cluster Energy Distribution;Energy [keV];Counts");
   // base filled bars
    hClusterE->SetFillColorAlpha(kGray+1, 0.6);
    hClusterE->SetLineColor(kGray+1);   // ignored; we'll draw borders separately
    hClusterE->SetBarWidth(0.95);
    hClusterE->SetBarOffset(0.025);
    hClusterE->Draw("BAR");

    // border-only clone
    auto hOutline = (TH1F*)hClusterE->Clone("hClusterE_outline");
    hOutline->SetFillStyle(0);          // no fill
    hOutline->SetLineColor(kBlack);     // black border
    hOutline->SetLineWidth(2);
    hOutline->Draw("BAR SAME");
    

    // Histogram already drawn: hClusterE->Draw("HIST");

    // Global range (for the total, if you want it)
    const double ene_cut_low  = hClusterE->GetXaxis()->GetXmin();
    const double ene_cut_high = hClusterE->GetXaxis()->GetXmax();

    // --- per-peak fit ranges ---
    const double g1_lo = 5.3, g1_hi = 6.1;  // Kα window (tweak as you like)
    const double g2_lo = 6.0, g2_hi = 6.8;  // Kβ window

    // ---- fit g1 only in [g1_lo, g1_hi] ----
    TF1* g1 = new TF1("g1","gaus", g1_lo, g1_hi);
    g1->SetParameters(hClusterE->GetMaximum(), 5.90, 0.10); // [A, mu, sigma]
    g1->SetParLimits(2, 0.03, 0.40);                        // keep sigma from collapsing
    hClusterE->Fit(g1, "RQ0S");

    // ---- fit g2 only in [g2_lo, g2_hi] ----
    TF1* g2 = new TF1("g2","gaus", g2_lo, g2_hi);
    g2->SetParameters(hClusterE->GetMaximum(), 6.40, 0.12);
    g2->SetParLimits(2, 0.03, 0.40);
    hClusterE->Fit(g2, "RQ0S");

    // ---- (optional) build total on full range using those seeds, then refit ----
    TF1* f2g = new TF1("f2g","gaus(0)+gaus(3)", ene_cut_low, ene_cut_high);
    f2g->SetParameters(g1->GetParameter(0), g1->GetParameter(1), g1->GetParameter(2),
                    g2->GetParameter(0), g2->GetParameter(1), g2->GetParameter(2));
    f2g->SetParLimits(2, 0.03, 0.40);
    f2g->SetParLimits(5, 0.03, 0.40);
    // hClusterE->Fit(f2g, "RQ0S");  // or skip this line if you only want separate fits

    // ---- draw components (each only over its own range) ----
    g1->SetLineColor(kRed+1);   g1->SetLineStyle(kSolid); g1->SetLineWidth(3); g1->SetNpx(600);
    g2->SetLineColor(kBlue+1);  g2->SetLineStyle(kSolid); g2->SetLineWidth(3); g2->SetNpx(600);
    g1->Draw("SAME");
    g2->Draw("SAME");
    // If you also want the total across the whole axis:
    // f2g->SetLineColor(kBlack); f2g->SetLineWidth(2); f2g->Draw("SAME");
    // If you want the total fit too, uncomment:
    // f2g->SetLineColor(kBlack); f2g->SetLineWidth(2); f2g->Draw("SAME");

    // --- Extract params & errors (prefer combined fit if available, else per-peak) ---
    bool have_total = (f2g != nullptr);

    double A1, mu1, s1, A2, mu2, s2;
    double A1_err, mu1_err, s1_err, A2_err, mu2_err, s2_err;

    A1 = g1->GetParameter(0);   mu1 = g1->GetParameter(1);   s1 = g1->GetParameter(2);
    A2 = g2->GetParameter(0);   mu2 = g2->GetParameter(1);   s2 = g2->GetParameter(2);
    A1_err = g1->GetParError(0); mu1_err = g1->GetParError(1); s1_err = g1->GetParError(2);
    A2_err = g2->GetParError(0); mu2_err = g2->GetParError(1); s2_err = g2->GetParError(2);

    // --- Manual chi2/ndf over the union of the two windows, using yfit = g1+g2 ---
    auto inG1 = [&](double x){ return (x >= g1_lo && x <= g1_hi); };
    auto inG2 = [&](double x){ return (x >= g2_lo && x <= g2_hi); };

    double chi2 = 0.0;
    int npts = 0;
    for (int ib = 1; ib <= hClusterE->GetNbinsX(); ++ib) {
        const double x = hClusterE->GetXaxis()->GetBinCenter(ib);
        if (!(inG1(x) || inG2(x))) continue;          // only evaluate in the fitted windows
        const double y    = hClusterE->GetBinContent(ib);
        const double yfit = g1->Eval(x) + g2->Eval(x); // sum model (works even if f2g not used)
        const double var  = y + 1e-6;                  // avoid 0 division
        chi2 += (y - yfit)*(y - yfit) / var;
        ++npts;
    }

    // 6 parameters total (3 per peak)
    const int npar = 6;
    const int ndf  = std::max(1, npts - npar);
    const double chi2_ndf = chi2 / ndf;
    const double prob = TMath::Prob(chi2, ndf);  // same tail prob as 1 - CDF in Python

    TLegend* legend = new TLegend(0.620267, 0.62087, 0.8942, 0.88);
    legend->SetBorderSize(1);
    legend->SetFillStyle(1001);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
    legend->SetTextAlign(12); // left-align
    legend->AddEntry(
        g1,
        Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}",
            mu1, mu1_err, s1, s1_err),
        "l"
    );
    legend->AddEntry(
        g2,
        Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}",
            mu2, mu2_err, s2, s2_err),
        "l"
    );
    // gray box for Data (use the histogram itself)
    legend->AddEntry(hClusterE, "Data", "f");
    legend->Draw();

    // Optional: match Python ticks/text sizes a bit
    hClusterE->GetXaxis()->SetLabelSize(0.045);
    hClusterE->GetYaxis()->SetLabelSize(0.045);
    hClusterE->GetXaxis()->SetTitleSize(0.05);
    hClusterE->GetYaxis()->SetTitleSize(0.05);

    // (Optional) save
    // c_cluster_E->SaveAs("energy_double_gauss_style.png");

    // Console summary (like your printouts)
    std::cout << "Double-Gaussian fit (Python-style bounds)\n"
            << "  A1=" << A1 << " +/- " << A1_err
            << "  mu1=" << mu1 << " +/- " << mu1_err
            << ", sigma1=" << s1 << " +/- " << s1_err << "\n"
            << "  A2=" << A2 << " +/- " << A2_err
            << "  mu2=" << mu2 << " +/- " << mu2_err
            << ", sigma2=" << s2 << " +/- " << s2_err << "\n"
            << "  chi2/ndf=" << chi2_ndf << ",  p=" << prob << "\n";

    c_cluster_E->SaveAs("cluster_energy_double_gauss.png");

    //---------------------------------------------//
    //--------      Additional Plots        -------//
    //---------------------------------------------//

    TCanvas* c_Posx_vs_sigma_xy = new TCanvas();
    c_Posx_vs_sigma_xy->cd();
    hposx_sigmaxy->SetTitle("Mean x vs #sigma_{xy}");
    hposx_sigmaxy->GetXaxis()->SetTitle("Mean x [column number]");
    hposx_sigmaxy->GetYaxis()->SetTitle("#sigma_{xy}[pixels]");
    hposx_sigmaxy->GetZaxis()->SetTitle("Counts");
    // hposx_sigmaxy->SetMarkersize(0.5);
    hposx_sigmaxy->SetMarkerStyle(20);
    hposx_sigmaxy->SetMarkerColor(kBlue);
    hposx_sigmaxy->SetLineColor(kBlue);
    hposx_sigmaxy->Draw("COLZ");
    c_Posx_vs_sigma_xy->SetLogz();

    TCanvas* c_Energy_vs_sigma_xy = new TCanvas();
    c_Energy_vs_sigma_xy->cd();
    henergy_sigmaxy->SetTitle("Cluster Energy vs #sigma_{xy}");
    henergy_sigmaxy->GetXaxis()->SetTitle("Cluster Energy [keV]");
    henergy_sigmaxy->GetYaxis()->SetTitle("#sigma_{xy} [pixels]");
    henergy_sigmaxy->GetZaxis()->SetTitle("Counts");
    // henergy_sigmaxy->SetMarkersize(0.5);
    henergy_sigmaxy->SetMarkerStyle(20);
    // henergy_sigmaxy->SetMarkerColor(kRed);
    // henergy_sigmaxy->SetLineColor(kRed);
    henergy_sigmaxy->Draw("COLZ");
    // c_Energy_vs_sigma_xy->SetLogz();

    // Create the canvas for cluster sigmaxy distribution
    TCanvas* c_sigmaxy = new TCanvas();
    c_sigmaxy->cd();
    hSigmaxy->SetTitle("Cluster #sigma_{xy} Distribution");
    hSigmaxy->GetXaxis()->SetTitle("#sigma_{xy} [pixels]");
    hSigmaxy->GetYaxis()->SetTitle("Counts");
    hSigmaxy->SetLineColor(kRed);
    hSigmaxy->SetFillColor(kRed);
    hSigmaxy->SetFillStyle(3001);
    hSigmaxy->Draw("hist"); 

    // Create the canvas for cluster sigmax distribution
    TCanvas* c_sigmax = new TCanvas();
    c_sigmax->cd();
    hSigmax->SetTitle("Cluster #sigma_{x} Distribution");
    hSigmax->GetXaxis()->SetTitle("#sigma_{x} [pixels]");
    hSigmax->GetYaxis()->SetTitle("Counts");
    hSigmax->SetLineColor(kRed);
    hSigmax->SetFillColor(kRed);
    hSigmax->SetFillStyle(3001);
    hSigmax->Draw("hist");

    // Create the canvas for cluster sigmay distribution
    TCanvas* c_sigmay = new TCanvas();
    c_sigmay->cd();
    hSigmay->SetTitle("Cluster #sigma_{y} Distribution");
    hSigmay->GetXaxis()->SetTitle("#sigma_{y} [pixels]");
    hSigmay->GetYaxis()->SetTitle("Counts");
    hSigmay->SetLineColor(kRed);
    hSigmay->SetFillColor(kRed);
    hSigmay->SetFillStyle(3001);
    hSigmay->Draw("hist");


    // create the canvas for cluster position distributions
    TCanvas* c_PosXY = new TCanvas();
    c_PosXY->cd();
    hPosXY->SetTitle("Cluster Position Distribution");
    hPosXY->GetXaxis()->SetTitle("X Position [column number]");
    hPosXY->GetYaxis()->SetTitle("Y Position [row number]");
    // hPosXY->SetMarkerStyle(20);
    // hPosXY->SetMarkerColor(kGreen);
    // hPosXY->SetLineColor(kGreen);
    hPosXY->Draw("BOXCOLZ");   
    // hPosX->Draw("hist");

    // create the canvas for cluster position distributions
    TCanvas* c_PosX = new TCanvas();
    c_PosX->cd();
    hPosX->SetTitle("Cluster Position Distribution");
    hPosX->GetXaxis()->SetTitle("X Position [column number]");
    hPosX->GetYaxis()->SetTitle("Counts");
    // hPosXY->SetMarkerStyle(20);
    // hPosXY->SetMarkerColor(kGreen);
    // hPosXY->SetLineColor(kGreen);
    // hPosXY->Draw("BOXCOLZ");   
    hPosX->Draw("hist");

    // create the canvas for cluster position distributions
    TCanvas* c_PosY = new TCanvas();
    c_PosY->cd();
    hPosY->SetTitle("Cluster Position Distribution");
    hPosY->GetXaxis()->SetTitle("Y Position [rows number]");
    hPosY->GetYaxis()->SetTitle("Counts");
    // hPosXY->SetMarkerStyle(20);
    // hPosXY->SetMarkerColor(kGreen);
    // hPosXY->SetLineColor(kGreen);
    // hPosXY->Draw("BOXCOLZ");   
    hPosY->Draw("hist");

    // Draw empty axes (without content)

    // Normalize marker size relative to the maximum bin content
    // double maxContent = hPosXY->GetMaximum();
    // double minSize = 0.2;
    // double maxSize = 2.0;
    // hPosXY->Draw("axis");

    // // Loop over bins
    // for (int ix = 1; ix <= hPosXY->GetNbinsX(); ix++) {
    //     for (int iy = 1; iy <= hPosXY->GetNbinsY(); iy++) {
    //         double content = hPosXY->GetBinContent(ix, iy);
    //         if (content <= 0) continue;

    //         double x = hPosXY->GetXaxis()->GetBinCenter(ix);
    //         double y = hPosXY->GetYaxis()->GetBinCenter(iy);

    //         // Scale marker size ∝ sqrt(content) so area reflects energy
    //         double norm = content / maxContent;
    //         double size = minSize + (maxSize - minSize) * sqrt(norm);

    //         TMarker *m = new TMarker(x, y, 20); // circle marker
    //         m->SetMarkerSize(size);
    //         m->SetMarkerColor(kRed+1);
    //         m->Draw("same");
    //     }
    // }
    // c_PosXY->SetLogz();


    

    // Create the canvas
    // TCanvas* c = new TCanvas("c", "Composite Cluster with Projections", 1000, 1000);
    // gStyle->SetOptStat(0);
    // c->Divide(2, 2);
    // c->cd(1); hcomp->Draw("COLZ");
    // c->cd(2); hPosX->Draw();
    // c->cd(3); hPosY->Draw();
}
