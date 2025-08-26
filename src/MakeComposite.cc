// // MakeComposite.cc — adapted to new ROOT format with aliases, front/back mode, and external file loading
// // Author: Diego 
// // IMportant note, thre is only one cluster per entry in every file, 

// #include "TH2F.h"
// #include "TH1D.h"
// #include "TFile.h"
// #include "TTree.h"
// #include "TChain.h"
// #include "TCanvas.h"
// #include "TSystemDirectory.h"
// // ROOT drawing helpers
// #include <TStyle.h>
// #include <TPaveText.h>
// #include <TLatex.h>
// #include <TLine.h>
// #include <TMath.h>
// #include <TH1.h>
// #include <TH2.h>
// #include <TCanvas.h>
// #include "TLegend.h"   // <-- Add this!


// // C++ utils used by AutoZRange / VectorsToTH2F
// #include <limits>
// #include <algorithm>
// #include "TString.h"  // TString, Form (since you use them)
// #include <vector>
// #include <iostream>
// #include <cmath>
// #include <string>
// #include <sstream>
// #include <fstream>
// // #include "TInterpreter.h"  // add this include

// // // ... inside MakeComposite(...) BEFORE binding branches:
// // gInterpreter->GenerateDictionary("std::vector<std::vector<int> >, std::vector<std::vector<float> >", "vector");
// #include <vector>
// #ifdef __CLING__
// #pragma link C++ class std::vector<int>+;
// #pragma link C++ class std::vector<float>+;
// #pragma link C++ class std::vector<std::vector<int> >+;
// #pragma link C++ class std::vector<std::vector<float> >+;
// #endif



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

// // Robust Z-range from quantiles (ignores zeros by default)
// void AutoZRange(TH2* h, double qlow=0.01, double qhigh=0.995, bool ignoreZeros=true) {
//     std::vector<double> vals;
//     vals.reserve(h->GetNbinsX()*h->GetNbinsY());
//     for (int ix=1; ix<=h->GetNbinsX(); ++ix) {
//         for (int iy=1; iy<=h->GetNbinsY(); ++iy) {
//             double c = h->GetBinContent(ix, iy);
//             if (ignoreZeros && c<=0) continue;
//             vals.push_back(c);
//         }
//     }
//     if (vals.size()<16) return; // too few nonzero bins; skip
//     std::sort(vals.begin(), vals.end());
//     auto q = [&](double f)->double{
//         if (f<=0) return vals.front();
//         if (f>=1) return vals.back();
//         const double idx = f*(vals.size()-1);
//         const size_t i0 = (size_t)std::floor(idx);
//         const size_t i1 = std::min(i0+1, vals.size()-1);
//         const double t = idx - i0;
//         return vals[i0]*(1.0-t) + vals[i1]*t;
//     };
//     const double zmin = q(qlow);
//     const double zmax = q(qhigh);
//     if (zmax>zmin) { h->SetMinimum(zmin); h->SetMaximum(zmax); }
// }

// void SetStandardAliases(TTree* tree) {
//     tree->SetAlias("cluster_E", "Energy");    // keV
//     tree->SetAlias("cluster_x", "PosX");
//     tree->SetAlias("cluster_y", "PosY");
//     tree->SetAlias("cluster_id_alias", "cluster_id");
//     tree->SetAlias("rmsxy", "wSTD_XY");       // or "STD_XY" if that’s your intended metric
//     tree->SetAlias("valid", "has_seed==1");
// }

// void MakeComposite(const char* listfile, const char* selection = "", int front = 1, double efact = 0.0029) {
//     TChain* tree = new TChain("clustersRec");
//     std::cout<<tree->GetEntries() << " entries in chain before adding files."<<std::endl;
//     std::ifstream infile(listfile);
//     std::string line;
//     int nadded = 0;
//     while (std::getline(infile, line)) {
//         // trim spaces and CR
//         line.erase(0, line.find_first_not_of(" \t\r"));
//         if (line.empty() || line[0]=='#') continue;
//         line.erase(line.find_last_not_of(" \t\r") + 1);
//         std::cout << "Adding file: " << line << std::endl;
//         nadded += tree->Add(line.c_str());
//     }
//     if (nadded <= 0) { std::cerr << "ERROR: no files added from list.\n"; return; }
//     std::cout<<tree->GetEntries() << " entries in chain after adding files."<<std::endl;


//     SetStandardAliases(tree);

//     std::vector<std::vector<int>>*   pixels_x = nullptr; // each entry has a vector of pixel x coordinates for the reconstructed cluster
//     std::vector<std::vector<int>>*   pixels_y = nullptr; // each entry has a vector of pixel y coordinates for the reconstructed cluster
//     std::vector<std::vector<float>>* pixels_E = nullptr; // each entry has a vector of pixel energies in keV for the reconstructed cluster
//     std::vector<float>*              PosX     = nullptr; // each entry has a 1 element vector with the mean x position of the cluster
//     std::vector<float>*              PosY     = nullptr; // each entry has a 1 element vector with the mean y position of the cluster
//     std::vector<int>*                Npix     = nullptr;
//     std::vector<float>*              wSTD_X   = nullptr; // each entry has a 1 element vectore with the standard deviation in x for the cluster
//     std::vector<float>*              wSTD_Y   = nullptr; // each entry has a 1 element vectore with the standard deviation in x for the cluster
//     std::vector<float>*              wSTD_XY   = nullptr; // each entry has a 1 element vectore with the standard deviation in xy for the cluster
//     std::vector<int>*                cluster_id = nullptr; 

//     tree->SetBranchAddress("pixels_x", &pixels_x);
//     tree->SetBranchAddress("pixels_y", &pixels_y);
//     tree->SetBranchAddress("pixels_E", &pixels_E);
//     tree->SetBranchAddress("PosX", &PosX);
//     tree->SetBranchAddress("PosY", &PosY);
//     tree->SetBranchAddress("wSTD_X", &wSTD_X);
//     tree->SetBranchAddress("wSTD_Y", &wSTD_Y);
//     tree->SetBranchAddress("cluster_id", &cluster_id);
//     tree->SetBranchAddress("wSTD_XY", &wSTD_XY);

    

//     //Get the entries that pass the selection in array v1
//     // Selection stage
//     tree->Draw("Entry$", selection, "GOFF");
//     Long64_t nsel = tree->GetSelectedRows();
//     std::cout << "Selected rows: " << nsel << " with selection: " << selection << std::endl;
//     Double_t* v1 = tree->GetV1();  // contains entry numbers
//     Int_t counter=0;
//     Int_t total = tree->GetSelectedRows();
//     std::vector<Long64_t> entries;    

//     std::cout << "Total entries selected = " << entries.size() << std::endl;
//     for (size_t i = 0; i < std::min<size_t>(entries.size(), 10); ++i)
//     std::cout << "entries[" << i << "] = " << entries[i] << std::endl;

    
    
    
//     // Set up composite histogram parameters
//     const int    nbins = 400;
//     const double bmin  = -4.0;
//     const double bmax  =  4.0;
//     const double bins_per_pix = 50.0; // for normalization
//     int ncl = 0;
//     double sum_skew = 0;

//     // Create histograms
//     TH2F* hcomp = new TH2F("composite", "composite", nbins, bmin, bmax, nbins, bmin, bmax);
//     TH1D* projX = new TH1D("projX", "Projection X", nbins, bmin, bmax);
//     TH1D* projY = new TH1D("projY", "Projection Y", nbins, bmin, bmax);
//     TH1D* hleft  = new TH1D("hleft",  "", nbins, bmin, bmax);
//     TH1D* hbelow = new TH1D("hbelow", "", nbins, bmin, bmax);
//     TCanvas* cview = new TCanvas("cview", "Composite View", 1200, 400);
//     cview->Divide(3, 1);

//     // Here We Process Each Entry
//     std::cout << "-------------------------" << std::endl;
//     std::cout << "---Processing entries ---" << std::endl;
//     std::cout << "-------------------------" << std::endl;
//     // while(counter<total && counter>=0){
//     for (Long64_t counter = 0; counter < nsel; ++counter) {
//         Long64_t entry = static_cast<Long64_t>(v1[counter]);  // ⚠️ This cast is critical

//         if (tree->GetEntry(entry) <= 0) {
//             std::cerr << "MakeComposite: GetEntry(" << entry << ") <= 0, skipping\n";
//             continue;
//         }
//     // for (size_t ie = 0; ie < entries.size(); ++ie) {
//     // for (size_t ie = 0; ie < 5; ++ie) {
//         // const Long64_t e = entries[ie];
//         // if (tree->GetEntry(v1[counter]) <= 0) {
//         //     std::cerr << "MakeComposite: GetEntry(" << counter << ") <= 0, skipping\n";
//         //     continue;
//         // }
//         // Long64_t entry = tree->GetEntryNumber(counter);
//         // if (tree->GetEntry(entry) <= 0) {
//         //     std::cerr << "MakeComposite: GetEntry(" << entry << ") <= 0, skipping\n";
//         //     continue;
//         // }
        

//         // Per-cluster loop
//         // int cluster_count = cluster_id ? cluster_id->size() : 0;
//         std::cout << "Entry " << counter << "  Nclusters = " << counter << "\n";
//         const auto& cur_pixel_x = pixels_x->at(0);
//         const auto& cur_pixel_y = pixels_y->at(0);
//         const auto& cur_pixel_energy = pixels_E->at(0);

//         if (cur_pixel_x.size() != cur_pixel_y.size() || cur_pixel_y.size() != cur_pixel_energy.size()) {
//             std::cerr << "MakeComposite: size mismatch for cluster " << counter
//                         << ": pixels_x=" << cur_pixel_x.size()
//                         << " pixels_y=" << cur_pixel_y.size()
//                         << " pixels_E=" << cur_pixel_energy.size() << "\n";
//             continue;
//         }
//         if (cur_pixel_x.empty()) {
//             std::cerr << "MakeComposite: empty pixel vectors for cluster " << cluster_id->at(0) << "\n";
//             continue;
//         }

//         // Tight image of this cluster in pixel coordinates
//         TH2F* hview = VectorsToTH2F(cur_pixel_x, cur_pixel_y, cur_pixel_energy, "h2_cluster");
//         if (!hview) continue;

//         const double mean_x = PosX->at(0);
//         const double mean_y = PosY->at(0);
        
//         std::cout<< "Cluster " << cluster_id->at(0) << ": mean_x=" << mean_x
//                     << ", mean_y=" << mean_y << "\n";

//         // Corner-charge veto (same idea as legacy): integrate outside 3×3 pixel square
//         double borq = 0.0;
//         for (int ix = 1; ix <= nbins; ++ix) {
//             for (int iy = 1; iy <= nbins; ++iy) {
//                 const double gx = mean_x + hcomp->GetXaxis()->GetBinLowEdge(ix);
//                 const double gy = mean_y + hcomp->GetYaxis()->GetBinLowEdge(iy);

//                 // const double gx = hcomp->GetXaxis()->GetBinCenter(ix) - mean_x;
//                 // const double gy = hcomp->GetYaxis()->GetBinCenter(iy) - mean_y;
//                 const double bc = hview->GetBinContent(hview->FindBin(gx, gy)) / (bins_per_pix * bins_per_pix);
//                 if (std::abs(hcomp->GetXaxis()->GetBinCenter(ix)) > 3.0 &&
//                     std::abs(hcomp->GetYaxis()->GetBinCenter(iy)) > 3.0) {
//                     borq += bc;
//                 }
//             }
//         }

//         if (std::abs(borq) < 25.0) {
//             // Accumulate into composite, centered on (mean_x,mean_y)
//             for (int ix = 1; ix <= nbins; ++ix) {
//                 for (int iy = 1; iy <= nbins; ++iy) {
//                     const double gx = mean_x + hcomp->GetXaxis()->GetBinLowEdge(ix);
//                     const double gy = mean_y + hcomp->GetYaxis()->GetBinLowEdge(iy);
//                     // const double gx = hcomp->GetXaxis()->GetBinCenter(ix) - mean_x;
//                     // const double gy = hcomp->GetYaxis()->GetBinCenter(iy) - mean_y;
//                     const double bc = hview->GetBinContent(hview->FindBin(gx, gy)) / (bins_per_pix * bins_per_pix);
//                     hcomp->SetBinContent(ix, iy, hcomp->GetBinContent(ix, iy) + bc);
//                 }
//             }
//             ++ncl;
//         }
//         if ((counter % 10) == 0) {
//             std::cout << "Processed entry " << (counter + 1) << " / " << entries.size() << std::endl;
//         }
//         counter++;
//         delete hview;
//     }

//     if (ncl > 0) {
//         hcomp->Scale(1.0 / ncl);
//         projX->Add(hcomp->ProjectionX());
//         projY->Add(hcomp->ProjectionY());
//     }

//     std::cout << "Composite built from " << ncl << " clusters.\n";
//     // Keep your legacy print (keV→e−): efact kept for compatibility, though pixels_E are already keV
//     std::cout << "Total integral: " << hcomp->Integral() * 1000.0 / 3.8 << " electrons\n";

//     // Style
//     gStyle->SetOptStat(0);
//     gStyle->SetNumberContours(100);
//     gStyle->SetPalette(kBird);

   

//     // Project charge in left and below regions
//     for (int ix = 1; ix <= hcomp->GetNbinsX(); ++ix) {
//         for (int iy = 1; iy <= hcomp->GetNbinsY(); ++iy) {
//             double gx = hcomp->GetXaxis()->GetBinCenter(ix);
//             double gy = hcomp->GetYaxis()->GetBinCenter(iy);
//             double val = hcomp->GetBinContent(ix, iy);

//             // Left region
//             if (gx >= -4 && gx <= -3)
//                 hleft->Fill(gy, val);

//             // Below region
//             if (gy >= -4 && gy <= -3)
//                 hbelow->Fill(gx, val);
//         }
//     }

//     std::cout << "Left integral: " << hleft->Integral() * 1000.0 / 3.8 << " electrons\n";
//     std::cout << "Below integral: " << hbelow->Integral() * 1000.0 / 3.8 << " electrons\n";


//     // // Normalize to unit area
//     // hleft->Scale(1.0 / hleft->Integral());
//     // hbelow->Scale(1.0 / hbelow->Integral());

//     // Compute stats
//     auto print_stats = [](TH1D* h, const TString& label) {
//         double mean = h->GetMean();
//         double rms = h->GetRMS();
//         double skew = h->GetSkewness();
//         printf("%s stats → Mean: %.3f | RMS: %.3f | Skewness: %.3f\n",
//             label.Data(), mean, rms, skew);
//         return std::make_tuple(mean, rms, skew);
//     };

//     auto [meanL, rmsL, skewL] = print_stats(hleft,  "Left ");
//     auto [meanB, rmsB, skewB] = print_stats(hbelow, "Below");

//     // Draw everything
//     TCanvas* c = new TCanvas("c", "", 1400, 600);
//     c->Divide(2,1);

//     c->cd(1);
//     hcomp->Draw("colz");

//     // Rectangle overlay (left and below regions)
//     TLine* boxL1 = new TLine(-3, -4, -3,  4);
//     TLine* boxL2 = new TLine(-4, -4, -4,  4);
//     TLine* boxL3 = new TLine(-4, -4, -3, -4);
//     TLine* boxL4 = new TLine(-4,  4, -3,  4);
//     TLine* boxB1 = new TLine(-4, -3,  4, -3);
//     TLine* boxB2 = new TLine(-4, -4,  4, -4);
//     TLine* boxB3 = new TLine(-4, -4, -4, -3);
//     TLine* boxB4 = new TLine( 4, -4,  4, -3);

//     for (auto ln : {boxL1, boxL2, boxL3, boxL4, boxB1, boxB2, boxB3, boxB4}) {
//         ln->SetLineColor(kGreen+2);
//         ln->SetLineStyle(2);
//         ln->Draw();
//     }

//     // Text labels
//     TLatex latex;
//     latex.SetTextSize(0.04);
//     latex.SetTextColor(kGreen+2);
//     latex.SetTextAlign(21);
//     latex.DrawLatex(-3.5,  4.2, "left");
//     latex.DrawLatex( 4.2, -3.5, "below");

//     // --- Stats Box
//     TPaveText* stats = new TPaveText(0.15, 0.65, 0.55, 0.85, "NDC");
//     stats->SetFillColor(kGreen-5);
//     stats->SetTextColor(kBlack);
//     stats->SetTextAlign(12);
//     stats->AddText("     Mean     RMS   Skewness");
//     stats->AddText(Form("Left   %.2f     %.2f     %.2f", meanL, rmsL, skewL));
//     stats->AddText(Form("Below  %.2f     %.2f     %.2f", meanB, rmsB, skewB));
//     stats->Draw();

//     // Plot projections
//     c->cd(2);
//     hleft->SetLineColor(kRed);
//     hleft->SetLineWidth(2);
//     hleft->SetTitle("hleft vs hbelow;Position [pix];Normalized charge");
//     hleft->Draw("hist");

//     hbelow->SetLineColor(kBlack);
//     hbelow->SetLineWidth(2);
//     hbelow->Draw("hist same");

//     TLegend* leg = new TLegend(0.6, 0.75, 0.85, 0.88);
//     leg->AddEntry(hleft,  "hleft",  "l");
//     leg->AddEntry(hbelow, "hbelow", "l");
//     leg->Draw();

//     c->SaveAs("composite_projection.pdf");

//     // // ---- styling
//     // gStyle->SetOptStat(0);
//     // gStyle->SetNumberContours(100);
//     // gStyle->SetPalette(kBird); // nice perceptual palette

//     // // titles & axes
//     // hcomp->SetTitle(Form("Composite (%zu entries, %d clusters);#Delta x [pix];#Delta y [pix]",
//     //                     entries.size(), ncl));
//     // hcomp->GetZaxis()->SetTitle("Charge (keV/pixel)");
//     // hcomp->GetZaxis()->SetTitleOffset(1.2);
//     // hcomp->GetXaxis()->CenterTitle(true);
//     // hcomp->GetYaxis()->CenterTitle(true);

//     // // auto Z range (robust to outliers)
//     // AutoZRange(hcomp, 0.01, 0.995, /*ignoreZeros=*/true);

//     // // draw on pad 1
//     // cview->cd(1);
//     // gPad->SetRightMargin(0.14);
//     // gPad->SetLogz(front!=0);
//     // hcomp->Draw("COLZ");

//     // // crosshairs at (0,0) and a 1-pixel bounding box
//     // TLine ln; ln.SetLineStyle(2); ln.SetLineColor(kGray+2);
//     // ln.DrawLine(0.0, hcomp->GetYaxis()->GetXmin(), 0.0, hcomp->GetYaxis()->GetXmax());
//     // ln.DrawLine(hcomp->GetXaxis()->GetXmin(), 0.0, hcomp->GetXaxis()->GetXmax(), 0.0);

//     // TLine box; box.SetLineColor(kBlack); box.SetLineStyle(1);
//     // box.DrawLine(-0.5, -0.5,  0.5, -0.5);
//     // box.DrawLine( 0.5, -0.5,  0.5,  0.5);
//     // box.DrawLine( 0.5,  0.5, -0.5,  0.5);
//     // box.DrawLine(-0.5,  0.5, -0.5, -0.5);

//     // // annotate total charge
//     // TLatex tx; tx.SetNDC();
//     // tx.SetTextSize(0.035);
//     // tx.DrawLatex(0.16, 0.92, Form("Integral: %.1f e^{-}", hcomp->Integral()*1000.0/3.8));

//     // // projections / front-back behavior
//     // if (front) {
//     //     // zoom to +/-4 pixels view already set by axis limits; nothing else
//     // } else {
//     //     cview->cd(2);
//     //     gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
//     //     projX->SetTitle(";#Delta x [pix];Charge (keV/pixel)");
//     //     projX->SetLineWidth(2);
//     //     projX->Draw("hist");

//     //     cview->cd(3);
//     //     gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
//     //     projY->SetTitle(";#Delta y [pix];Charge (keV/pixel)");
//     //     projY->SetLineWidth(2);
//     //     projY->Draw("hist");
//     // }

//     // // force refresh and save an image (works with -b -q too)
//     // cview->Modified(); cview->Update();
//     // const char* tag = front ? "front" : "back";
//     // cview->SaveAs(Form("composite_%s.pdf", tag));
//     // // also keep the histogram for later reuse if you like:
//     // // TFile fout("composite.root","RECREATE"); hcomp->Write(); projX->Write(); projY->Write(); fout.Close();

//     // tree->ResetBranchAddresses();




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

int MakeComposite(const char* filelist = "test_list.txt", double min_energy = 4.5) {
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
    std::vector<float>* cluster_E = nullptr;

    chain.SetBranchAddress("pixels_x", &pixels_x);
    chain.SetBranchAddress("pixels_y", &pixels_y);
    chain.SetBranchAddress("pixels_E", &pixels_E);
    chain.SetBranchAddress("PosX", &PosX);
    chain.SetBranchAddress("PosY", &PosY);
    chain.SetBranchAddress("Energy", &cluster_E);

    if (chain.GetEntries() != 1) {
        std::cerr << "Expected exactly one entry in clustersRec tree, found: " << chain.GetEntries() << std::endl;
        return 1;
    }

    chain.GetEntry(0);
    int nClusters = PosX->size();
    std::cout << "Found " << nClusters << " clusters in entry 0\n";

    const int nbins = 400; // Number of bins for the composite histogram
    const double range = 4.0;
    TH2D* hcomp = new TH2D("hcomp", "Composite Cluster", nbins, -range, range, nbins, -range, range);
    TH1D* hleft = new TH1D("hleft", "Left Projection", nbins, -range, range);
    TH1D* hbelow = new TH1D("hbelow", "Below Projection", nbins, -range, range);

    int nAccepted = 0;
    for (size_t i = 0; i < PosX->size(); ++i) {
        if (cluster_E->at(i) < min_energy)
            continue;

        double cx = PosX->at(i);
        double cy = PosY->at(i);
        const auto& px = pixels_x->at(i);
        const auto& py = pixels_y->at(i);
        const auto& pe = pixels_E->at(i);

        for (size_t j = 0; j < px.size(); ++j) {
            double dx = px[j] - cx;
            double dy = py[j] - cy;
            double val = pe[j];
            hcomp->Fill(dx, dy, val);
        }

        ++nAccepted;
    }

    std::cout << "Composite built from " << nAccepted << " clusters\n";

    for (int i = 1; i <= nbins; ++i) {
        for (int j = 1; j <= nbins; ++j) {
            double content = hcomp->GetBinContent(i, j);
            hleft->AddBinContent(i, content);
            hbelow->AddBinContent(j, content);
        }
    }

    hleft->Scale(1.0 / hleft->Integral());
    hbelow->Scale(1.0 / hbelow->Integral());

    // Create the canvas
    TCanvas* c = new TCanvas("c", "Composite Cluster with Projections", 1000, 1000);
    gStyle->SetOptStat(0);
    c->Divide(2, 2);
    c->cd(2); hbelow->Draw();
    c->cd(3); hleft->Draw();
    c->cd(1); hcomp->Draw("COLZ");

    // Stats
    double mean_l = hleft->GetMean(), rms_l = hleft->GetRMS(), skew_l = compute_skewness(hleft);
    double mean_b = hbelow->GetMean(), rms_b = hbelow->GetRMS(), skew_b = compute_skewness(hbelow);
    double int_l = hleft->Integral(), int_b = hbelow->Integral(), int_tot = hcomp->Integral();

    c->cd(4);
    TPaveText* pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    pt->AddText(Form("Composite built from %d clusters", nAccepted));
    pt->AddText(Form("Total integral: %.1f electrons", int_tot));
    pt->AddText(Form("Left integral: %.1f electrons", int_l));
    pt->AddText(Form("Below integral: %.1f electrons", int_b));
    pt->AddText(Form("Left  stats → Mean: %.3f | RMS: %.3f | Skewness: %.3f", mean_l, rms_l, skew_l));
    pt->AddText(Form("Below stats → Mean: %.3f | RMS: %.3f | Skewness: %.3f", mean_b, rms_b, skew_b));
    pt->Draw();

    c->SaveAs("composite_projection.pdf");
    std::cout << "Composite saved to composite_projection.pdf\n";

    return 0;
}


// void MakeComposite_Fe55_Clusters(const char* filelist = "test_list.txt", double min_energy = 4.5){
//     TChain chain("clustersRec");
//     std::ifstream fin(filelist);
//     std::string filename;
//     while (std::getline(fin, filename)) {
//         std::cout << "Adding file: " << filename << std::endl;
//         chain.Add(filename.c_str());
//     }
//     // Set up branch pointers for the main variables
//     // There is an interesting point about the format of this data, 
//     // There is only a single entry in the TTree, which contains all clusters.
//     // So we can thing of a single entry as a full CCD image
//     std::vector<std::vector<int>>* pixels_x = nullptr;
//     std::vector<std::vector<int>>* pixels_y = nullptr;
//     std::vector<std::vector<float>>* pixels_E = nullptr;
//     std::vector<float>* PosX = nullptr;
//     std::vector<float>* PosY = nullptr;
//     std::vector<float>* cluster_E = nullptr;
//     chain.SetBranchAddress("pixels_x", &pixels_x);
//     chain.SetBranchAddress("pixels_y", &pixels_y);
//     chain.SetBranchAddress("pixels_E", &pixels_E);
//     chain.SetBranchAddress("PosX", &PosX);
//     chain.SetBranchAddress("PosY", &PosY);
//     chain.SetBranchAddress("Energy", &cluster_E);
//     if (chain.GetEntries() != 1) {
//         std::cerr << "Expected exactly one entry in clustersRec tree, found: " << chain.GetEntries() << std::endl;
//         // break;
//     }
//     chain.GetEntry(0);
//     int nClusters = PosX->size();
//     std::cout << "Found " << nClusters << " clusters in entry 0\n";

//     const int nbins = 21; // Number of bins for the composite histogram
//     const double range = 4.5;
//     const int bins_per_pix = 50; // for normalization
//     const double efact = 3.8; // energy factor for conversion to electrons
//     TH2D* hcomp = new TH2D("hcomp", "Composite Cluster", nbins, -range, range, nbins, -range, range);
//     // TH1D* hleft = new TH1D("hleft", "Left Projection", nbins, -range, range);
//     // TH1D* hbelow = new TH1D("hbelow", "Below Projection", nbins, -range, range);

//     int nAccepted = 0;
//     for (size_t n = 0; n < nClusters; n++) {
//         if (cluster_E->at(n) < min_energy)
//             continue;

//         double cur_cluster_x = PosX->at(n);
//         double cur_cluster_y = PosY->at(n);
//         const auto& cur_pixel_x = pixels_x->at(n);
//         const auto& cur_pixel_y = pixels_y->at(n);
//         const auto& cur_pixel_energy = pixels_E->at(n);
//         double borq = 0;


//         TH2F* hview = VectorsToTH2F(cur_pixel_x, cur_pixel_y, cur_pixel_energy, "h2_cluster");
//         if (!hview) continue;

//         //Not the must efficient or best way to do this, but add the charge in the four corner pixels
//         //Used below to remove other events appearing within the 8x8 window and contaminating the cluster
//         for(int i=1; i<=nbins; i++){
//             for(int j=1; j<=nbins; j++){
//                 // double bc = hview->GetBinContent(hview->FindBin(mx+hcomp->GetXaxis()->GetBinLowEdge(i),my+hcomp->GetYaxis()->GetBinLowEdge(j)))/50./50.; //divide by bin size to keep integral constant
//                 // if(TMath::Abs(hcomp->GetXaxis()->GetBinCenter(i)) > 3 && TMath::Abs(hcomp->GetYaxis()->GetBinCenter(j)) > 3)
//                 //     borq+=bc;
//                 const double gx = cur_cluster_x + hcomp->GetXaxis()->GetBinCenter(i);
//                 const double gy = cur_cluster_y + hcomp->GetYaxis()->GetBinCenter(j);
//                 const double bc = hview->GetBinContent(hview->FindBin(gx, gy)) / (bins_per_pix * bins_per_pix);
//                 if (std::abs(hcomp->GetXaxis()->GetBinCenter(i)) > 3.0 &&
//                     std::abs(hcomp->GetYaxis()->GetBinCenter(j)) > 3.0) {
//                     borq += bc;
//                 }
//             }
//         }

//         if (std::abs(borq) < 25.0) {//selection on the maximum charge in the corner pixels
//             // Accumulate into composite, centered on (mean_x,mean_y)
//             // for (int i = 1; i <= nbins; ++i) {
//             //     for (int j = 1; j <= nbins; ++j) {
//             //         const double gx = cur_cluster_x + hcomp->GetXaxis()->GetBinCenter(i);
//             //         const double gy = cur_cluster_y + hcomp->GetYaxis()->GetBinCenter(j);
//             //         const double bc = hview->GetBinContent(hview->FindBin(gx, gy))/ (bins_per_pix * bins_per_pix);
//             //         hcomp->SetBinContent(i, j, hcomp->GetBinContent(i, j) + bc);
//             //     }
//             // }
//             for (size_t k = 0; k < cur_pixel_x.size(); ++k) {
//                 double gx = cur_pixel_x.at(k) - cur_cluster_x;  // Δx relative to cluster center
//                 double gy = cur_pixel_y.at(k) - cur_cluster_y;  // Δy relative to cluster center
//                 double energy = cur_pixel_energy.at(k); // / (bins_per_pix * bins_per_pix);  // Normalize
            
//                 hcomp->Fill(gx, gy, energy);  // Fill composite histogram
//             }
//             nAccepted++; //count number of clusters to make composite
//         }

//         if ((n % 10) == 0) {
//             std::cout << "Processed cluster " << (n + 1) << " / " << nClusters << std::endl;
//         }

//         hview->Delete();

//     }

//     int mid_bin = nbins / 2 + 1;
//     std::cout << "Center bin X: " << hcomp->GetXaxis()->GetBinCenter(mid_bin) << std::endl;
//     std::cout << "Center bin Y: " << hcomp->GetYaxis()->GetBinCenter(mid_bin) << std::endl;

//     hcomp->Scale(1./nAccepted); //scale by number of clusters to average composite
//     //Print integrated, average cluster charge in electrons
//     std::cout << "Total integral: " << hcomp->Integral()*efact*1000/3.8 << std::endl;
//     // std::cout << "Found " << nClusters << " clusters in entry 0\n";
//     hcomp->Draw("COLZ");

// }