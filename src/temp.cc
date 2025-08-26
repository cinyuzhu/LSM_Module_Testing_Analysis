// Skipper CCD cluster reader in the style of the legacy composite script (LSM - AEC style)
// - Opens a ROOT file (or wildcard) into a TChain
// - Reads per-entry, per-cluster data from the 'clustersRec' tree
// - Fields read: cluster_id, Energy, PosX, PosY, pixels_x, pixels_y, pixels_E,
//                Npix, wSTD_X, wSTD_Y, STD_XY
//
// Build:
//   g++ -std=c++17 read_clusters_like_old.cc `root-config --cflags --libs` -o read_clusters_like_old
//
// Run (all entries):
//   ./read_clusters_like_old /path/to/file.root
//
// With selection (same style as TTree::Draw):
//   ./read_clusters_like_old "/path/to/*.root" "Energy>0 && Npix>3"
//
// Custom tree name:
//   ./read_clusters_like_old /path/to/file.root "" "clustersRec"
//
// Notes:
// - This mirrors the structure of your ArraysToTH2F / MakeComposite approach,
//   but binds to std::vector / std::vector<std::vector<...>> in the new schema.
// - The VectorsToTH2F helper is a drop-in analog to ArraysToTH2F for quick inspection.

#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

// Old method analog: turn per-cluster vectors into a small TH2F around the cluster
// Includes only cluster pixels (no extra border here to keep it minimal)
static TH2F* VectorsToTH2F(const std::vector<int>& x,
                           const std::vector<int>& y,
                           const std::vector<float>& val,
                           const char* name = "h2")
{
    if (x.size() != y.size() || y.size() != val.size()) {
        ::Error("VectorsToTH2F", "Size mismatch x=%zu y=%zu val=%zu",
                x.size(), y.size(), val.size());
        return nullptr;
    }
    if (x.empty()) {
        return new TH2F(name, name, 1, -0.5, 0.5, 1, -0.5, 0.5);
    }

    int min_x = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int min_y = std::numeric_limits<int>::max();
    int max_y = std::numeric_limits<int>::min();

    for (size_t i = 0; i < x.size(); ++i) {
        min_x = std::min(min_x, x[i]);
        max_x = std::max(max_x, x[i]);
        min_y = std::min(min_y, y[i]);
        max_y = std::max(max_y, y[i]);
    }

    // Bin edges mirror the legacy style: one bin per integer pixel, tight bounding box
    const int nx = (max_x - min_x + 1);
    const int ny = (max_y - min_y + 1);

    // Unique name per call to avoid directory clashes
    static int hcounter = 0;
    TString hname = Form("%s_%d", name, hcounter++);
    TH2F* h2 = new TH2F(hname, hname, nx, min_x - 0.5, max_x + 0.5,
                                     ny, min_y - 0.5, max_y + 0.5);

    for (size_t i = 0; i < x.size(); ++i) {
        const int bx = h2->GetXaxis()->FindBin(x[i]);
        const int by = h2->GetYaxis()->FindBin(y[i]);
        // If the same (x,y) appears twice (shouldn't), we add
        h2->SetBinContent(bx, by, h2->GetBinContent(bx, by) + val[i]);
    }
    return h2;
}

// Adapted from your MakeComposite pattern: iterate over selected entries,
// and for each entry, iterate clusters and print the requested fields.
// Selection is a standard TTree::Draw selection string (e.g., "Energy>0 && Npix>3").
// tree_name defaults to "clustersRec".
void DumpSelectedClusters(TChain* c,
                          const TString& selection = "",
                          size_t max_clusters_to_print_per_entry = 5)
{
    // Bind branches (new schema)
    std::vector<int>*                cluster_id = nullptr;
    std::vector<float>*              Energy     = nullptr;
    std::vector<float>*              PosX       = nullptr;
    std::vector<float>*              PosY       = nullptr;
    std::vector<std::vector<int>>*   pixels_x   = nullptr;
    std::vector<std::vector<int>>*   pixels_y   = nullptr;
    std::vector<std::vector<float>>* pixels_E   = nullptr;
    std::vector<int>*                Npix       = nullptr;
    std::vector<float>*              wSTD_X     = nullptr;
    std::vector<float>*              wSTD_Y     = nullptr;
    std::vector<float>*              STD_XY     = nullptr;

    auto require = [&](const char* bname, void* addr) {
        if (!c->GetBranch(bname)) {
            ::Error("DumpSelectedClusters", "Missing required branch '%s'", bname);
            throw std::runtime_error(Form("Missing required branch '%s'", bname));
        }
        if (c->SetBranchAddress(bname, addr) != 0) {
            ::Error("DumpSelectedClusters", "Failed to SetBranchAddress '%s'", bname);
            throw std::runtime_error(Form("Failed to SetBranchAddress '%s'", bname));
        }
    };

    require("cluster_id", &cluster_id);
    require("Energy",     &Energy);
    require("PosX",       &PosX);
    require("PosY",       &PosY);
    require("pixels_x",   &pixels_x);
    require("pixels_y",   &pixels_y);
    require("pixels_E",   &pixels_E);
    require("Npix",       &Npix);
    require("wSTD_X",     &wSTD_X);
    require("wSTD_Y",     &wSTD_Y);
    require("STD_XY",     &STD_XY);

    // Build list of entries matching the selection (same pattern you use)
    Long64_t nsel = 0;
    std::vector<Long64_t> selected_entries;
    if (selection.Length() > 0) {
        c->Draw("Entry$", selection.Data(), "GOFF");
        nsel = c->GetSelectedRows();
        Double_t* v1 = c->GetV1();
        selected_entries.reserve(nsel);
        for (Long64_t i = 0; i < nsel; ++i) selected_entries.push_back(static_cast<Long64_t>(v1[i]));
    } else {
        // No selection: process all entries
        nsel = c->GetEntries();
        selected_entries.reserve(nsel);
        for (Long64_t i = 0; i < nsel; ++i) selected_entries.push_back(i);
    }

    std::cout << "Selected entries: " << nsel << "\n";

    // Loop over selected entries
    for (Long64_t idx = 0; idx < nsel; ++idx) {
        const Long64_t entry = selected_entries[idx];
        if (c->GetEntry(entry) <= 0) {
            ::Warning("DumpSelectedClusters", "GetEntry(%lld) returned <=0, skipping", entry);
            continue;
        }

        const size_t ncl = cluster_id->size();
        auto check = [&](const char* nm, size_t sz) {
            if (sz != ncl) {
                throw std::runtime_error(Form("Size mismatch at entry %lld: '%s'=%zu vs cluster_id=%zu",
                                              entry, nm, sz, ncl));
            }
        };
        check("Energy",   Energy->size());
        check("PosX",     PosX->size());
        check("PosY",     PosY->size());
        check("Npix",     Npix->size());
        check("wSTD_X",   wSTD_X->size());
        check("wSTD_Y",   wSTD_Y->size());
        check("STD_XY",   STD_XY->size());
        check("pixels_x", pixels_x->size());
        check("pixels_y", pixels_y->size());
        check("pixels_E", pixels_E->size());

        std::cout << "\nEntry " << entry << "  Nclusters = " << ncl << "\n";

        const size_t print_n = std::min(max_clusters_to_print_per_entry, ncl);
        for (size_t i = 0; i < print_n; ++i) {
            const auto& xs = pixels_x->at(i);
            const auto& ys = pixels_y->at(i);
            const auto& Es = pixels_E->at(i);

            const size_t nx = xs.size();
            const size_t ny = ys.size();
            const size_t nE = Es.size();
            const int npix  = (i < Npix->size()) ? Npix->at(i) : -1;

            if (!(nx == ny && ny == nE && nE == static_cast<size_t>(npix))) {
                ::Warning("DumpSelectedClusters",
                          "entry %lld cluster %zu: pixel sizes nx=%zu ny=%zu nE=%zu Npix=%d",
                          entry, i, nx, ny, nE, npix);
            }

            std::cout << "  i=" << std::setw(5) << i
                      << "  id=" << std::setw(6) << cluster_id->at(i)
                      << "  E(keV)=" << std::setw(10) << Energy->at(i)
                      << "  Pos=(" << PosX->at(i) << "," << PosY->at(i) << ")"
                      << "  Npix=" << npix
                      << "  wSTD=(" << wSTD_X->at(i) << "," << wSTD_Y->at(i) << ")"
                      << "  STD_XY=" << STD_XY->at(i)
                      << "\n";

            // Show first few pixels (x,y,E)
            const size_t kmax = std::min<size_t>(Es.size(), 5);
            std::cout << "    pixels [first " << kmax << "]: ";
            for (size_t k = 0; k < kmax; ++k) {
                std::cout << "(" << xs[k] << "," << ys[k] << "," << Es[k] << ") ";
            }
            if (Es.size() > kmax) std::cout << "...";
            std::cout << "\n";

            // Optional quick visual like your ArraysToTH2F (disabled by default)
            // TH2F* hview = VectorsToTH2F(xs, ys, Es, "h2_cluster");
            // if (hview) { new TCanvas(); hview->Draw("COLZ"); }
        }

        if (ncl > print_n) {
            std::cout << "  ... (" << (ncl - print_n) << " more clusters in this entry)\n";
        }
        if ((idx % 100) == 0) {
            std::cout << "Processed " << (idx + 1) << " / " << nsel << " selected entries\n";
        }
    }

    c->ResetBranchAddresses();
}

int main(int argc, char** argv)
{
    // Args:
    //   1) input (file or wildcard for TChain::Add)
    //   2) [optional] selection string (TTree::Draw syntax)
    //   3) [optional] tree name (default "clustersRec")
    //
    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " <input.root|pattern> [selection] [tree_name]\n\n"
                  << "Examples:\n"
                  << "  " << argv[0] << " /data/run123.root\n"
                  << "  " << argv[0] << " \"/data/run*.root\" \"Energy>0 && Npix>3\"\n"
                  << "  " << argv[0] << " /data/run123.root \"\" Clusters\n";
        return 1;
    }

    const TString inpat   = argv[1];
    const TString select  = (argc >= 3) ? argv[2] : "";
    const TString treename= (argc >= 4) ? argv[3] : "clustersRec";

    // Build chain in your legacy style
    TChain chain(treename);
    const int nadded = chain.Add(inpat);
    if (nadded <= 0) {
        std::cerr << "ERROR: No files matched/added for pattern: " << inpat << "\n";
        return 2;
    }
    std::cout << "Added " << nadded << " file(s) to chain '" << treename << "'.\n";

    try {
        DumpSelectedClusters(&chain, select, /*max_clusters_to_print_per_entry=*/5);
    } catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 3;
    }
    return 0;
}
