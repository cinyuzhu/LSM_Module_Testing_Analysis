#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>

void CheckClustersPerEntry(const char* filename = "/Users/diegovenegasvargas/Documents/LSM_Module_Testing/Production_Modules/2025-08-07-DM01_DM02_DM03_Low_Temp/Analysis/Image_4/panaSKImg_clustersRec_avg_Image_4_Low_Temp_106_20250807_101152_8_ext1.root", 
    const char* treename = "clustersRec") {
TFile* file = TFile::Open(filename);
if (!file || file->IsZombie()) {
std::cerr << "Error opening file.\n";
return;
}

TTree* tree = (TTree*)file->Get(treename);
if (!tree) {
std::cerr << "Error getting tree.\n";
return;
}

std::vector<float>* PosX = nullptr;
std::vector<float>* PosY = nullptr;
std::vector<int>* Npix = nullptr;

tree->SetBranchAddress("PosX", &PosX);
tree->SetBranchAddress("PosY", &PosY);
tree->SetBranchAddress("Npix", &Npix);

Long64_t nentries = tree->GetEntries();
std::cout << "Total entries: " << nentries << std::endl;
std::cout << "Checking clusters per entry...\n";
for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    if (PosX && PosY) {
        std::cout << "Entry " << i << ": PosX size = " << PosX->size() 
                  << ", PosY size = " << PosY->size() 
                  << ", Npix size = " << (Npix ? Npix->size() : 0) << std::endl;
    } else {
        std::cout << "Entry " << i << ": Missing PosX or PosY data.\n";
    }
}
std::cout << "Checking first 20 entries:\n";

// for (Long64_t i = 0; i < nentries; ++i) {
// tree->GetEntry(i);
// std::cout << "Entry " << i<< ", PosX->size() = " << PosX->size() << std::endl;
// << ", PosY->size() = " << PosY->size() << std::endl;
// }

file->Close();
}
