// MakeComposite_Fe55_Clusters.C
// ROOT >=6 recommended

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TMath.h"

// --- CLING vector dictionaries (if running as macro) ---
#ifdef __CLING__
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#endif

// ------------------------ helpers ------------------------
TH2F* VectorsToTH2F(const std::vector<int>& x,
                    const std::vector<int>& y,
                    const std::vector<float>& val,
                    const char* name="h2")
{
  if (x.empty()) return new TH2F(name,name,1,-0.5,0.5,1,-0.5,0.5);
  int min_x = *std::min_element(x.begin(), x.end());
  int max_x = *std::max_element(x.begin(), x.end());
  int min_y = *std::min_element(y.begin(), y.end());
  int max_y = *std::max_element(y.begin(), y.end());
  TH2F* h2 = new TH2F(name, name,
                      max_x - min_x + 1, min_x - 0.5, max_x + 0.5,
                      max_y - min_y + 1, min_y - 0.5, max_y + 0.5);
  for (size_t i=0;i<x.size();++i) h2->Fill(x[i], y[i], val[i]);
  return h2;
}

static double SumRect(const TH2D* h, double xlo, double xhi, double ylo, double yhi) {
  const double eps=1e-6;
  auto ax=h->GetXaxis(), ay=h->GetYaxis();
  int x1=ax->FindBin(xlo+eps), x2=ax->FindBin(xhi-eps);
  int y1=ay->FindBin(ylo+eps), y2=ay->FindBin(yhi-eps);
  return const_cast<TH2D*>(h)->Integral(x1,x2,y1,y2);
}

struct ExtResults {
  // bookkeeping
  int    nFiles=0;
  long   totalClusters=0;
  int    nAcceptedEnergy=0, nAcceptedFront=0, nAcceptedBack=0;

  // composites & spectra
  TH2D*  hBack=nullptr;
  TH2D*  hFront=nullptr;
  TH1F*  hEnergy=nullptr;
  TH1D*  projX=nullptr; // ProjectionX over Y:-4→-3
  TH1D*  projY=nullptr; // ProjectionY over X:-4→-3

  // CTI fractions (front)
  double centerCharge=0, fAbove=0, fBelow=0, fRight=0, fLeft=0;

  // back-side projection metrics (pre-normalization, electrons)
  double leftInt_e=0, belowInt_e=0;
  double leftMean=0, leftRMS=0, leftSkew=0;
  double belowMean=0, belowRMS=0, belowSkew=0;

  // energy fit (two per-peak Gaussians)
  double mu1=0, mu1_err=0, s1=0, s1_err=0;
  double mu2=0, mu2_err=0, s2=0, s2_err=0;
  double chi2=0, p=0; int ndf=0;

  // integrals (e-)
  double frontIntegral_e=0, backIntegral_e=0;
};

static void styled_axes() {
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);
  gStyle->SetPalette(kBird);
}

// ---------------------- main routine ----------------------
void MakeComposite_Fe55_Clusters_Grouped(const char* filelist = "test_list.txt",
                                 const char* moduleUsed = "DMXX",
                                 int image_number = 4)
{
  styled_axes();

  // ---------- read file list & group by extension ----------
  std::ifstream fin(filelist);
  if (!fin) { std::cerr<<"[ERROR] cannot open "<<filelist<<"\n"; return; }
  std::vector<std::string> byExt[5]; // 1..4
  std::string line;
  while (std::getline(fin,line)) {
    if (line.empty()) continue;
    for (int e=1;e<=4;++e) {
      std::string tag = std::string("_ext")+std::to_string(e)+".root";
      if (line.rfind(tag)!=std::string::npos || line.find(tag)!=std::string::npos) {
        byExt[e].push_back(line);
        break;
      }
    }
  }

  // ---------- image geometry choices ----------
  int nColumns=6400, nRows=1600;
  if      (image_number==5){ nColumns=1280; nRows=320; }
  else if (image_number==6){ nColumns=640;  nRows=1600;}
  else if (image_number==7){ nColumns=640;  nRows=1600;}

  // ---------- per-extension results ----------
  ExtResults R[5];

  // ---------- base output dirs ----------
  TString baseTag = Form("Module_%s_Image_%d", moduleUsed, image_number);
  gSystem->mkdir(baseTag, kTRUE);

  // loop over extensions
  for (int ext=1; ext<=4; ++ext) {
    if (byExt[ext].empty()) continue;
    R[ext].nFiles = (int)byExt[ext].size();

    // per-ext output dir
    TString extDir = Form("%s/ext%d", baseTag.Data(), ext);
    gSystem->mkdir(extDir, kTRUE);

    // chain
    TChain ch("clustersRec");
    for (auto& f : byExt[ext]) ch.Add(f.c_str());

    // branches
    std::vector<std::vector<int>>   *pixels_x=nullptr,*pixels_y=nullptr;
    std::vector<std::vector<float>> *pixels_E=nullptr;
    std::vector<float> *PosX=nullptr,*PosY=nullptr,*cluster_E=nullptr;
    std::vector<float> *rmsxy=nullptr,*rmsx=nullptr,*rmsy=nullptr;

    ch.SetBranchAddress("pixels_x", &pixels_x);
    ch.SetBranchAddress("pixels_y", &pixels_y);
    ch.SetBranchAddress("pixels_E", &pixels_E);
    ch.SetBranchAddress("PosX", &PosX);
    ch.SetBranchAddress("PosY", &PosY);
    ch.SetBranchAddress("Energy", &cluster_E);
    ch.SetBranchAddress("wSTD_XY", &rmsxy);
    ch.SetBranchAddress("fSTD_X",  &rmsx);
    ch.SetBranchAddress("fSTD_Y",  &rmsy);

    // histos
    const int nbins=401; const double range=4; const int bins_per_pix=25;
    TH2D* hBack  = new TH2D(Form("hcomp_back_ext%d",ext),"Back Composite;#Delta x [pix];#Delta y [pix]", nbins,-range,range, nbins,-range,range);
    TH2D* hFront = new TH2D(Form("hcomp_front_ext%d",ext),"Front Composite;#Delta x [pix];#Delta y [pix]", nbins,-range,range, nbins,-range,range);
    TH1F* hE     = new TH1F(Form("hClusterE_ext%d",ext), "Cluster Energy;Energy [keV];Counts", 50, 4.5, 7.0);

    // selection
    double minE=4.5, maxE=7.0;
    double sigma_xy_front_cut=1.0;

    Long64_t nEnt = ch.GetEntries();
    for (Long64_t ie=0; ie<nEnt; ++ie) {
      ch.GetEntry(ie);
      const int N = PosX->size();
      R[ext].totalClusters += N;

      for (int n=0;n<N;++n) {
        const double E = (*cluster_E)[n];
        const double x = (*PosX)[n];
        const double y = (*PosY)[n];
        const double sx = (*rmsx)[n];
        const double sy = (*rmsy)[n];
        const double sxy = std::sqrt(0.5*(sx*sx + sy*sy));

        const auto& px = (*pixels_x)[n];
        const auto& py = (*pixels_y)[n];
        const auto& pe = (*pixels_E)[n];

        if (E>minE && E<maxE) {
          R[ext].nAcceptedEnergy++;
        }

        // choose front/back discriminator per image_number
        double discr = sxy;
        if      (image_number==6) discr = sx;
        else if (image_number==7) discr = sy;

        if (E>minE && E<maxE && discr < sigma_xy_front_cut) {
          hE->Fill(E);
          TH2F* hv = VectorsToTH2F(px,py,pe,"_tmpF");
          for (int i=1;i<=nbins;++i){
            for (int j=1;j<=nbins;++j){
              const double gx = x + hFront->GetXaxis()->GetBinLowEdge(i);
              const double gy = y + hFront->GetYaxis()->GetBinLowEdge(j);
              const double bc = hv->GetBinContent(hv->FindBin(gx,gy))/(bins_per_pix*bins_per_pix);
              hFront->AddBinContent(hFront->GetBin(i,j), bc);
            }
          }
          R[ext].nAcceptedFront++;
          delete hv;
        }
        if (E>minE && E<maxE && discr >= sigma_xy_front_cut) {
          TH2F* hv = VectorsToTH2F(px,py,pe,"_tmpB");
          for (int i=1;i<=nbins;++i){
            for (int j=1;j<=nbins;++j){
              const double gx = x + hBack->GetXaxis()->GetBinLowEdge(i);
              const double gy = y + hBack->GetYaxis()->GetBinLowEdge(j);
              const double bc = hv->GetBinContent(hv->FindBin(gx,gy))/(bins_per_pix*bins_per_pix);
              hBack->AddBinContent(hBack->GetBin(i,j), bc);
            }
          }
          R[ext].nAcceptedBack++;
          delete hv;
        }
      }
    }

    // normalize composites
    if (R[ext].nAcceptedFront>0) hFront->Scale(1.0/R[ext].nAcceptedFront);
    if (R[ext].nAcceptedBack >0) hBack ->Scale(1.0/R[ext].nAcceptedBack);

    // integrals in electrons
    const double ePerADU=3.8;
    R[ext].frontIntegral_e = hFront->Integral()*1000.0/ePerADU;
    R[ext].backIntegral_e  = hBack ->Integral()*1000.0/ePerADU;

    // ---------- back-side plots & projections ----------
    {
      TCanvas* cBack = new TCanvas(Form("cBack_ext%d",ext),"Back-side", 1100, 480);
      cBack->Divide(2,1);

      // left: heatmap
      cBack->cd(1);
      gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
      gPad->SetLogz();
      hBack->GetXaxis()->SetRangeUser(-4,4);
      hBack->GetYaxis()->SetRangeUser(-4,4);
      hBack->Draw("COLZ");

      // right: projections
      cBack->cd(2);
      gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
      // Y:-4→-3 == bins 51..100 (given nbins=401 over [-4,4])
      TH1D* hPX = hBack->ProjectionX(Form("hbelow_ext%d",ext), 51,100);
      TH1D* hPY = hBack->ProjectionY(Form("hleft_ext%d",ext),  51,100);

      // --- print/record metrics BEFORE scaling (as requested) ---
      R[ext].leftInt_e  = hPY->Integral()*1000.0/ePerADU;
      R[ext].belowInt_e = hPX->Integral()*1000.0/ePerADU;
      R[ext].leftMean = hPY->GetMean();  R[ext].leftRMS = hPY->GetRMS();  R[ext].leftSkew = hPY->GetSkewness();
      R[ext].belowMean= hPX->GetMean();  R[ext].belowRMS= hPX->GetRMS();  R[ext].belowSkew= hPX->GetSkewness();

      // style & normalize for overlay
      hPX->SetLineColor(kBlack);   hPX->SetLineWidth(3);  hPX->SetFillColorAlpha(kBlack,0.10);
      hPY->SetLineColor(kBlue+1);  hPY->SetLineWidth(3);  hPY->SetFillColorAlpha(kBlue+1,0.10);

      if (hPX->Integral()>0) hPX->Scale(1.0/hPX->Integral());
      if (hPY->Integral()>0) hPY->Scale(1.0/hPY->Integral());

      hPX->SetTitle(";Pixel offset [pix];Normalized counts");
      hPX->GetYaxis()->SetNdivisions(505);
      hPX->SetMinimum(0.0);
      hPX->SetMaximum(1.15*std::max(hPX->GetMaximum(),hPY->GetMaximum()));
      hPX->Draw("HIST");
      hPY->Draw("HIST SAME");

      auto leg = new TLegend(0.57,0.72,0.89,0.88);
      leg->AddEntry(hPX,"ProjX (Y: -4 #rightarrow -3 px)","l");
      leg->AddEntry(hPY,"ProjY (X: -4 #rightarrow -3 px)","l");
      leg->Draw();

      cBack->SaveAs(Form("%s/backside_ext%d.png", extDir.Data(), ext));

      // store for compare panel
      R[ext].hBack = (TH2D*)hBack->Clone(Form("hBack_keep_ext%d",ext));
      R[ext].projX = (TH1D*)hPX->Clone(Form("hPX_keep_ext%d",ext));
      R[ext].projY = (TH1D*)hPY->Clone(Form("hPY_keep_ext%d",ext));
    }

    // ---------- front-side composite + CTI fractions ----------
    {
      // fractions from fixed boxes (center +/-0.5 etc.)
      const double centrp = SumRect(hFront, -0.5,  0.5, -0.5,  0.5);
      const double abovep = SumRect(hFront, -0.5,  1.5,  0.5,  2.5);
      const double belowp = SumRect(hFront, -0.5,  0.5, -1.5, -0.5);
      const double rightp = SumRect(hFront,  0.5,  1.5, -0.5,  0.5);
      const double leftp  = SumRect(hFront, -1.5, -0.5, -0.5,  0.5);

      R[ext].centerCharge = centrp;
      R[ext].fAbove = (centrp>0)? abovep/centrp : 0;
      R[ext].fBelow = (centrp>0)? belowp/centrp : 0;
      R[ext].fRight = (centrp>0)? rightp/centrp : 0;
      R[ext].fLeft  = (centrp>0)? leftp /centrp : 0;

      TCanvas* cF = new TCanvas(Form("cFront_ext%d",ext),"Front-side", 1100, 480);
      cF->Divide(2,1);

      // left: heatmap with squares
      cF->cd(1);
      gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
      gPad->SetLogz();
      hFront->GetXaxis()->SetRangeUser(-4,4);
      hFront->GetYaxis()->SetRangeUser(-4,4);
      hFront->SetTitle("Front-side Composite;#Delta x [pix];#Delta y [pix]");
      hFront->Draw("COLZ");

      auto drawBox = [&](double x1,double y1,double x2,double y2, Color_t lc){
        TBox* b = new TBox(x1,y1,x2,y2);
        b->SetLineColor(lc); b->SetLineWidth(3); b->SetFillStyle(0); b->Draw("same");
      };
      drawBox(-0.5,-0.5, 0.5, 0.5, kBlack);
      drawBox(-0.5, 0.5, 0.5, 1.5, kBlue+1);
      drawBox(-0.5,-1.5, 0.5,-0.5, kRed+1);
      drawBox( 0.5,-0.5, 1.5, 0.5, kGreen+2);
      drawBox(-1.5,-0.5,-0.5, 0.5, kMagenta);

      // right: colored bars
      cF->cd(2);
      gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.17);
      TH1F* frame = new TH1F("frame",";Region;Fraction (w.r.t. center)",5,0.5,5.5);
      frame->GetXaxis()->SetBinLabel(1,"Above");
      frame->GetXaxis()->SetBinLabel(2,"Below");
      frame->GetXaxis()->SetBinLabel(3,"Right");
      frame->GetXaxis()->SetBinLabel(4,"Left");
      frame->GetXaxis()->SetBinLabel(5,"Center=1.0");
      frame->SetMinimum(0.0);
      frame->SetMaximum(1.25*std::max(1.0, std::max({R[ext].fAbove,R[ext].fBelow,R[ext].fRight,R[ext].fLeft})));
      frame->Draw("AXIS");

      TLine* ref = new TLine(0.5,1.0,5.5,1.0); ref->SetLineStyle(2); ref->SetLineColor(kGray+2); ref->Draw();

      auto bar = [&](int bin,double y,Color_t c){
        double x1=frame->GetXaxis()->GetBinLowEdge(bin), x2=frame->GetXaxis()->GetBinUpEdge(bin);
        TBox* b=new TBox(x1,0.0,x2,y); b->SetFillColorAlpha(c,0.35); b->SetLineColor(c); b->SetLineWidth(3); b->Draw("same");
      };
      bar(1,R[ext].fAbove,kBlue+1);
      bar(2,R[ext].fBelow,kRed+1);
      bar(3,R[ext].fRight,kGreen+2);
      bar(4,R[ext].fLeft ,kMagenta);
      bar(5,1.0,kGray+1);

      auto leg = new TLegend(0.18,0.58,0.55,0.88);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->AddEntry((TObject*)0,Form("Center charge = %.0f", R[ext].centerCharge), "");
      leg->AddEntry((TObject*)0,Form("Above = %.3f", R[ext].fAbove),  "");
      leg->AddEntry((TObject*)0,Form("Below = %.3f", R[ext].fBelow),  "");
      leg->AddEntry((TObject*)0,Form("Right = %.3f", R[ext].fRight),  "");
      leg->AddEntry((TObject*)0,Form("Left  = %.3f", R[ext].fLeft ),  "");
      leg->Draw();

      cF->SaveAs(Form("%s/frontside_ext%d.png", extDir.Data(), ext));

      R[ext].hFront = (TH2D*)hFront->Clone(Form("hFront_keep_ext%d",ext));
    }

    // ---------- energy spectrum (front-side) with two Gaussians ----------
    {
      TCanvas* cE = new TCanvas(Form("cE_ext%d",ext),"Energy", 900,600);
      gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

      // bar style: gray filled with black outline
      hE->SetFillColorAlpha(kGray+1,0.60);
      hE->SetLineColor(kGray+2);
      hE->SetBarWidth(0.95);
      hE->SetBarOffset(0.025);
      hE->Draw("BAR");
      auto outline = (TH1F*)hE->Clone(Form("hE_out_ext%d",ext));
      outline->SetFillStyle(0); outline->SetLineColor(kBlack); outline->SetLineWidth(2);
      outline->Draw("BAR SAME");

      // per-peak windows & fits
      const double g1_lo=5.3, g1_hi=6.1;
      const double g2_lo=6.0, g2_hi=6.8;

      TF1* g1 = new TF1(Form("g1_ext%d",ext),"gaus", g1_lo, g1_hi);
      g1->SetParameters(hE->GetMaximum(),5.90,0.10); g1->SetParLimits(2,0.03,0.40);
      hE->Fit(g1,"RQ0S");

      TF1* g2 = new TF1(Form("g2_ext%d",ext),"gaus", g2_lo, g2_hi);
      g2->SetParameters(hE->GetMaximum()/3.0,6.40,0.12); g2->SetParLimits(2,0.03,0.40);
      hE->Fit(g2,"RQ0S");

      g1->SetLineColor(kRed+1);  g1->SetLineWidth(3);
      g2->SetLineColor(kBlue+1); g2->SetLineWidth(3);
      g1->Draw("SAME"); g2->Draw("SAME");

      // store params
      R[ext].mu1 = g1->GetParameter(1); R[ext].mu1_err = g1->GetParError(1);
      R[ext].s1  = g1->GetParameter(2); R[ext].s1_err  = g1->GetParError(2);
      R[ext].mu2 = g2->GetParameter(1); R[ext].mu2_err = g2->GetParError(1);
      R[ext].s2  = g2->GetParameter(2); R[ext].s2_err  = g2->GetParError(2);

      // chi2/ndf over union windows using sum of components
      auto inG1=[&](double x){return x>=g1_lo && x<=g1_hi;};
      auto inG2=[&](double x){return x>=g2_lo && x<=g2_hi;};
      double chi2=0; int npts=0;
      for (int ib=1; ib<=hE->GetNbinsX(); ++ib){
        double x=hE->GetXaxis()->GetBinCenter(ib);
        if (!(inG1(x)||inG2(x))) continue;
        double y=hE->GetBinContent(ib);
        double yfit=g1->Eval(x)+g2->Eval(x);
        double var=y+1e-6;
        chi2 += (y-yfit)*(y-yfit)/var;
        ++npts;
      }
      R[ext].ndf = std::max(1,npts-6);
      R[ext].chi2=chi2; R[ext].p = TMath::Prob(chi2,R[ext].ndf);

      // legend in μ/σ format like your Python
      auto leg = new TLegend(0.60,0.63,0.89,0.88);
      leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetTextFont(42); leg->SetTextSize(0.035);
      leg->AddEntry(g1,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}",
                            R[ext].mu1,R[ext].mu1_err,R[ext].s1,R[ext].s1_err),"l");
      leg->AddEntry(g2,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}",
                            R[ext].mu2,R[ext].mu2_err,R[ext].s2,R[ext].s2_err),"l");
      leg->AddEntry(hE,"Data","f");
      leg->Draw();

      cE->SaveAs(Form("%s/energy_ext%d.png", extDir.Data(), ext));

      R[ext].hEnergy = (TH1F*)hE->Clone(Form("hE_keep_ext%d",ext));
    }

    // ---------- per-extension log (detailed) ----------
    {
      std::ofstream l(Form("%s/extension.log", Form("%s/ext%d",baseTag.Data(),ext)));
      l << "---------------- Extension ext" << ext << " ----------------\n";
      l << "Files combined: " << R[ext].nFiles << "\n";
      l << "Total clusters: " << R[ext].totalClusters << "\n";
      l << "Accepted (energy window): " << R[ext].nAcceptedEnergy << "\n";
      l << "Accepted (front-side):    " << R[ext].nAcceptedFront  << "\n";
      l << "Accepted (back-side):     " << R[ext].nAcceptedBack   << "\n";
      l << "Front integral [e-]: " << R[ext].frontIntegral_e << "\n";
      l << "Back  integral [e-]: " << R[ext].backIntegral_e  << "\n";
      l << "CTI fractions (front, w.r.t center)\n";
      l << "  center="<<R[ext].centerCharge<<"  above="<<R[ext].fAbove
        << "  below="<<R[ext].fBelow<<"  right="<<R[ext].fRight<<"  left="<<R[ext].fLeft<<"\n";
      l << "Back-side projection metrics (pre-normalization)\n";
      l << "  hleft integral [e-]: "  << R[ext].leftInt_e  << "\n";
      l << "  hleft  mean/RMS/skew: " << R[ext].leftMean  << "  " << R[ext].leftRMS  << "  " << R[ext].leftSkew  << "\n";
      l << "  hbelow integral [e-]: " << R[ext].belowInt_e << "\n";
      l << "  hbelow mean/RMS/skew: " << R[ext].belowMean << "  " << R[ext].belowRMS << "  " << R[ext].belowSkew << "\n";
      l << "Energy double-Gaussian (per-peak windows)\n";
      l << "  mu1="<<R[ext].mu1<<" +/- "<<R[ext].mu1_err<<"   sigma1="<<R[ext].s1<<" +/- "<<R[ext].s1_err<<"\n";
      l << "  mu2="<<R[ext].mu2<<" +/- "<<R[ext].mu2_err<<"   sigma2="<<R[ext].s2<<" +/- "<<R[ext].s2_err<<"\n";
      l << "  chi2/ndf="<<R[ext].chi2<<"/"<<R[ext].ndf<<"   p="<<R[ext].p<<"\n";
    }

  } // end per-extension loop

  // -------------------- comparison canvases --------------------
  gSystem->mkdir(Form("%s/compare", baseTag.Data()), kTRUE);

  // global z range for back/front
  double zminB=1e30, zmaxB=0, zminF=1e30, zmaxF=0;
  double ymaxProj=0, ymaxEnergy=0;
  for (int e=1;e<=4;++e){
    if (!R[e].hBack) continue;
    zmaxB = std::max(zmaxB, R[e].hBack->GetMaximum());
    if (R[e].hBack->GetMinimum(0)>0) zminB = std::min(zminB, R[e].hBack->GetMinimum(0));
    zmaxF = std::max(zmaxF, R[e].hFront->GetMaximum());
    if (R[e].hFront->GetMinimum(0)>0) zminF = std::min(zminF, R[e].hFront->GetMinimum(0));
    if (R[e].projX) ymaxProj = std::max({ymaxProj, R[e].projX->GetMaximum(), R[e].projY->GetMaximum()});
    if (R[e].hEnergy) ymaxEnergy = std::max(ymaxEnergy, (double)R[e].hEnergy->GetMaximum());
  }
  if (!(zminB<zmaxB)) { zminB=1e-6; zmaxB=1; }
  if (!(zminF<zmaxF)) { zminF=1e-6; zmaxF=1; }

  // ---- (1) Back-side comparison 2x2 ----
  {
    TCanvas* c = new TCanvas("compare_back","Back 2x2", 1400, 1000);
    c->Divide(2,2);
    int pad=1;
    for (int e=1;e<=4;++e){
      if (!R[e].hBack) continue;
      c->cd(pad++);

      // split each tile horizontally
      TPad* left = new TPad(Form("bL_%d",e),"",0.0,0.0,0.58,1.0); left->Draw();
      TPad* right= new TPad(Form("bR_%d",e),"",0.58,0.0,1.0,1.0); right->Draw();

      left->cd(); left->SetRightMargin(0.16); left->SetLeftMargin(0.12); left->SetBottomMargin(0.12); left->SetLogz();
      TH2D* h = (TH2D*)R[e].hBack->Clone();
      h->GetXaxis()->SetRangeUser(-4,4); h->GetYaxis()->SetRangeUser(-4,4);
      h->SetMinimum(zminB); h->SetMaximum(zmaxB);
      h->SetTitle("Back-side Composite;#Delta x [pix];#Delta y [pix]");
      h->Draw("COLZ");

      right->cd(); right->SetGridx(); right->SetGridy(); right->SetLeftMargin(0.16); right->SetBottomMargin(0.16);
      TH1D* px = (TH1D*)R[e].projX->Clone(); TH1D* py=(TH1D*)R[e].projY->Clone();
      px->SetMaximum(1.15*ymaxProj); px->SetMinimum(0.0);
      px->SetTitle(";Pixel offset [pix];Normalized counts");
      px->Draw("HIST"); py->Draw("HIST SAME");
    }
    c->SaveAs(Form("%s/compare/compare_back_2x2.png", baseTag.Data()));
  }

  // ---- (2) Front-side comparison 2x2 ----
  {
    TCanvas* c = new TCanvas("compare_front","Front 2x2", 1400, 1000);
    c->Divide(2,2);
    int pad=1;
    for (int e=1;e<=4;++e){
      if (!R[e].hFront) continue;
      c->cd(pad++);
      TPad* left = new TPad(Form("fL_%d",e),"",0.0,0.0,0.58,1.0); left->Draw();
      TPad* right= new TPad(Form("fR_%d",e),"",0.58,0.0,1.0,1.0); right->Draw();

      left->cd(); left->SetRightMargin(0.16); left->SetLeftMargin(0.12); left->SetBottomMargin(0.12); left->SetLogz();
      TH2D* h=(TH2D*)R[e].hFront->Clone();
      h->GetXaxis()->SetRangeUser(-4,4); h->GetYaxis()->SetRangeUser(-4,4);
      h->SetMinimum(zminF); h->SetMaximum(zmaxF);
      h->SetTitle("Front-side Composite;#Delta x [pix];#Delta y [pix]");
      h->Draw("COLZ");
      auto box=[&](double x1,double y1,double x2,double y2,Color_t c){ TBox*b=new TBox(x1,y1,x2,y2); b->SetLineColor(c); b->SetLineWidth(3); b->SetFillStyle(0); b->Draw("same"); };
      box(-0.5,-0.5,0.5,0.5,kBlack);
      box(-0.5,0.5,0.5,1.5,kBlue+1);
      box(-0.5,-1.5,0.5,-0.5,kRed+1);
      box(0.5,-0.5,1.5,0.5,kGreen+2);
      box(-1.5,-0.5,-0.5,0.5,kMagenta);

      right->cd(); right->SetGridx(); right->SetGridy(); right->SetLeftMargin(0.16); right->SetBottomMargin(0.18);
      TH1F* frame=new TH1F("frameF",";Region;Fraction (w.r.t. center)",5,0.5,5.5);
      frame->GetXaxis()->SetBinLabel(1,"Above"); frame->GetXaxis()->SetBinLabel(2,"Below");
      frame->GetXaxis()->SetBinLabel(3,"Right"); frame->GetXaxis()->SetBinLabel(4,"Left");
      frame->GetXaxis()->SetBinLabel(5,"Center=1.0");
      double ymax = 1.25*std::max(1.0, std::max({R[e].fAbove,R[e].fBelow,R[e].fRight,R[e].fLeft}));
      frame->SetMinimum(0.0); frame->SetMaximum(ymax); frame->Draw("AXIS");
      TLine* ref=new TLine(0.5,1.0,5.5,1.0); ref->SetLineStyle(2); ref->SetLineColor(kGray+2); ref->Draw();
      auto bar=[&](int bin,double y,Color_t c){ double x1=frame->GetXaxis()->GetBinLowEdge(bin),x2=frame->GetXaxis()->GetBinUpEdge(bin); TBox*b=new TBox(x1,0.0,x2,y); b->SetFillColorAlpha(c,0.35); b->SetLineColor(c); b->SetLineWidth(3); b->Draw("same"); };
      bar(1,R[e].fAbove,kBlue+1); bar(2,R[e].fBelow,kRed+1); bar(3,R[e].fRight,kGreen+2); bar(4,R[e].fLeft,kMagenta); bar(5,1.0,kGray+1);
    }
    c->SaveAs(Form("%s/compare/compare_front_2x2.png", baseTag.Data()));
  }

  // ---- (3) Energy spectra comparison 2x2 ----
  {
    TCanvas* c = new TCanvas("compare_energy","Energy 2x2", 1400, 1000);
    c->Divide(2,2);
    int pad=1;
    for (int e=1;e<=4;++e){
      if (!R[e].hEnergy) continue;
      c->cd(pad++); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
      TH1F* h=(TH1F*)R[e].hEnergy->Clone();
      h->SetTitle("Cluster Energy;Energy [keV];Counts");
      h->SetFillColorAlpha(kGray+1,0.60); h->SetLineColor(kGray+2); h->SetBarWidth(0.95); h->SetBarOffset(0.025);
      h->SetMaximum(1.15*ymaxEnergy);
      h->Draw("BAR");
      auto outline=(TH1F*)h->Clone(); outline->SetFillStyle(0); outline->SetLineColor(kBlack); outline->SetLineWidth(2); outline->Draw("BAR SAME");

      // refit quick per-tile to draw components (use stored means as seeds if available)
      TF1* g1=new TF1(Form("cEg1_%d",e),"gaus",5.3,6.1); g1->SetParameters(h->GetMaximum(), (R[e].mu1?R[e].mu1:5.9), std::max(0.08, R[e].s1)); g1->SetParLimits(2,0.03,0.40); h->Fit(g1,"RQ0S");
      TF1* g2=new TF1(Form("cEg2_%d",e),"gaus",6.0,6.8); g2->SetParameters(h->GetMaximum()/3.0,(R[e].mu2?R[e].mu2:6.4), std::max(0.08, R[e].s2)); g2->SetParLimits(2,0.03,0.40); h->Fit(g2,"RQ0S");
      g1->SetLineColor(kRed+1); g1->SetLineWidth(3); g2->SetLineColor(kBlue+1); g2->SetLineWidth(3);
      g1->Draw("SAME"); g2->Draw("SAME");

      auto leg=new TLegend(0.58,0.63,0.89,0.88);
      leg->AddEntry(g1,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu1,R[e].mu1_err,R[e].s1,R[e].s1_err),"l");
      leg->AddEntry(g2,Form("#splitline{#mu = %.3f #pm %.3f keV}{#sigma = %.3f #pm %.3f keV}", R[e].mu2,R[e].mu2_err,R[e].s2,R[e].s2_err),"l");
      leg->AddEntry(h,"Data","f"); leg->Draw();
    }
    c->SaveAs(Form("%s/compare/compare_energy_2x2.png", baseTag.Data()));
  }

  // -------------------- summary log --------------------
  {
    std::ofstream log(Form("%s/summary.log", baseTag.Data()));
    log << "Summary for " << baseTag << "  (Module: " << moduleUsed << ")\n";
    log << "================================================================\n";
    for (int e=1; e<=4; ++e){
      if (R[e].nFiles==0) continue;
      log << "----------------------------------------------------------------\n";
      log << "Extension: ext" << e << "\n";
      log << "Files combined: " << R[e].nFiles << "\n";
      log << "Total clusters: " << R[e].totalClusters << "\n";
      log << "Accepted (energy window): " << R[e].nAcceptedEnergy << "\n";
      log << "Accepted (front-side):    " << R[e].nAcceptedFront  << "\n";
      log << "Accepted (back-side):     " << R[e].nAcceptedBack   << "\n";
      log << "Front integral [e-]: " << R[e].frontIntegral_e << "\n";
      log << "Back  integral [e-]: " << R[e].backIntegral_e  << "\n";
      log << "CTI fractions for front events (w.r.t. center)\n";
      log << "  center charge: " << R[e].centerCharge << "\n";
      log << "  above=" << R[e].fAbove
          << " below=" << R[e].fBelow
          << " right=" << R[e].fRight
          << " left="  << R[e].fLeft  << "\n";
      log << "Back-side projection metrics (pre-normalization)\n";
      log << "  hleft  integral [e-]: " << R[e].leftInt_e  << "\n";
      log << "  hleft  mean/RMS/skew: " << R[e].leftMean   << " " << R[e].leftRMS   << " " << R[e].leftSkew   << "\n";
      log << "  hbelow integral [e-]: " << R[e].belowInt_e << "\n";
      log << "  hbelow mean/RMS/skew: " << R[e].belowMean  << " " << R[e].belowRMS  << " " << R[e].belowSkew  << "\n";
      log << "Energy fit (two separate Gaussians in per-peak windows)\n";
      log << "  mu1=" << R[e].mu1 << " +/- " << R[e].mu1_err
          << "  sigma1=" << R[e].s1 << " +/- " << R[e].s1_err << "\n";
      log << "  mu2=" << R[e].mu2 << " +/- " << R[e].mu2_err
          << "  sigma2=" << R[e].s2 << " +/- " << R[e].s2_err << "\n";
      log << "  chi2/ndf=" << R[e].chi2 << "/" << R[e].ndf
          << "  p=" << R[e].p << "\n";
    }
  }

  std::cout << "\n[OK] Outputs under: " << baseTag << "/\n";
}
