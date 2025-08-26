// MakeComposite.cc
// Updated for TTree with single entry containing all clusters

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./MakeComposite <input_root_file>" << std::endl;
        return 1;
    }

    const char* input_file = argv[1];

    TFile* f = TFile::Open(input_file);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open file: " << input_file << std::endl;
        return 1;
    }

    TTree* tree = (TTree*) f->Get("clustersRec");
    if (!tree) {
        std::cerr << "TTree 'clustersRec' not found in file." << std::endl;
        return 1;
    }

    // Set up branches
    std::vector<std::vector<int>>*   pixels_x = nullptr;
    std::vector<std::vector<int>>*   pixels_y = nullptr;
    std::vector<std::vector<float>>* pixels_E = nullptr;
    std::vector<float>*              PosX     = nullptr;
    std::vector<float>*              PosY     = nullptr;

    tree->SetBranchAddress("pixels_x", &pixels_x);
    tree->SetBranchAddress("pixels_y", &pixels_y);
    tree->SetBranchAddress("pixels_E", &pixels_E);
    tree->SetBranchAddress("PosX", &PosX);
    tree->SetBranchAddress("PosY", &PosY);

    // Create composite histogram
    const int nbins = 101;
    const double range = 50.5;
    TH2D* hcomp = new TH2D("hcomp", "Composite Cluster;#Delta x;#Delta y", nbins, -range, range, nbins, -range, range);

    // Load the single entry
    if (tree->GetEntries() != 1 || tree->GetEntry(0) <= 0) {
        std::cerr << "Expected one entry in tree. Found: " << tree->GetEntries() << std::endl;
        return 1;
    }

    int ncl = 0;
    int n_clusters = PosX->size();
    std::cout << "Processing " << n_clusters << " clusters..." << std::endl;

    for (int ic = 0; ic < n_clusters; ++ic) {
        double mx = PosX->at(ic);
        double my = PosY->at(ic);
        double borq = 0;

        std::vector<int>& px = pixels_x->at(ic);
        std::vector<int>& py = pixels_y->at(ic);
        std::vector<float>& pe = pixels_E->at(ic);

        if (px.size() != py.size() || px.size() != pe.size()) {
            std::cerr << "Cluster " << ic << " pixel vector size mismatch. Skipping.\n";
            continue;
        }

        for (size_t ip = 0; ip < px.size(); ++ip) {
            double dx = px[ip] - mx;
            double dy = py[ip] - my;
            double e = pe[ip];

            if (std::abs(dx) > 3 && std::abs(dy) > 3) borq += e;
        }

        if (std::abs(borq) < 25) {
            for (size_t ip = 0; ip < px.size(); ++ip) {
                double dx = px[ip] - mx;
                double dy = py[ip] - my;
                double e = pe[ip];
                hcomp->Fill(dx, dy, e);
            }
            ncl++;
            std::cout << "Added cluster " << ic << " with PosX = " << mx << ", PosY = " << my << std::endl;
        }
    }

    std::cout << "Composite built from " << ncl << " clusters.\n";

    // Draw and save
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "", 800, 700);
    hcomp->Draw("colz");
    c->SaveAs("composite_cluster.pdf");

    delete f;
    return 0;
}
