#include "interface/CompLooperBase.h"
#include "RatioPlotter.h"
#include <iostream>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "ROOT/RDataFrame.hxx"

//NOTE: Only accepts 2022 skimmed ntuples
class CompLooper : public CompLooperBase {
  public:
    CompLooper(const char *name, const char *channel, const char *label1, const char *label2);
    ~CompLooper();

    void SetNorm(bool norm=true){_norm = norm;}
		void SetMode(std::string mode){_mode = mode;}
    void SetMakePlots(bool val=true);
    void SetMakeRatios(bool val=true);
    void SetPlotFiletype(std::string ft=".png");
    void Loop();
  private:
    std::string FindFile(const char *label);
    static constexpr const char *CHANNEL = "eemm";
    bool _makePlots = false, _makeRatios = false;
    bool _norm = false;
    std::string _label1, _label2;
    std::string _filetype = ".png", _mode = "UPDATE";
};

CompLooper::CompLooper(const char *name, const char *channel, const char *label1, const char *label2)
    : _label1(label1), _label2(label2),
      CompLooperBase(name, channel, 
          ((std::string)"slimmed/" + label1 + ".root").c_str(), 
          ((std::string)"slimmed/" + label2 + ".root").c_str()
      ){
};

CompLooper::~CompLooper(){}

void CompLooper::SetMakePlots(bool val){_makePlots = val;}
void CompLooper::SetMakeRatios(bool val){_makeRatios = val;}
void CompLooper::SetPlotFiletype(std::string ft){
  if (ft.size() == 0) return;
  if (ft[0] != '.') ft = "." + ft;
  _filetype = ft;
}

void CompLooper::Loop(){
	if (_infile1 == nullptr || _infile2 == nullptr){
		std::cerr << "There was a problem loading one of the files" << std::endl;
		return;
	}
  Init();

  // Options
  bool plot = false; //override for _makePlots
  bool ratioplot = false; //override for _makeRatios

  // Setup
  TFile *histout = new TFile(("out/histout_" + _name + ".root").c_str(), _mode.c_str());
  if (!histout->IsOpen()) return;
  histout->cd();
	TCanvas *c = new TCanvas("c", "Histogram Canvas");
  c->cd();
  bool makeDirs = (plot || ratioplot || _makePlots || _makeRatios);
	
  // Optionally creating plot directories
  std::string dirname = "plots/" + _name;
  if (makeDirs){
    DIR *plotdir = opendir("plots");
    if (!plotdir) mkdir("plots", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);

    plotdir = opendir(dirname.c_str());
    if (!plotdir) mkdir(dirname.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);
  }

  // Set values based on options
  std::string l1, l2, l3, l4;
  if (_channel == "eeee"){
    l1 = "e1"; l2 = "e2"; l3 = "e3"; l4 = "e4";
  }
  else if (_channel == "eemm"){
    l1 = "e1"; l2 = "e2"; l3 = "m1"; l4 = "m2";
  }
  else if (_channel == "mmmm"){
    l1 = "m1"; l2 = "m2"; l3 = "m3"; l4 = "m4";
  }
  
  // Binning
  std::vector<Double_t> InvMass_4l_binning   = {100.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0};
  std::vector<Double_t> InvMass_pair_binning = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120.};
  std::vector<Double_t> LepPt_binning        = {0., 30., 50., 100., 200.};

  // Declaring Histograms
  TH1F *InvMass4l_1 = new TH1F("InvMass4l_1", TString::Format("4-Lepton Invariant Mass (%s)", _label1.c_str()), InvMass_4l_binning.size()-1, &InvMass_4l_binning[0]);
  TH1F *Z2Mass_1 = new TH1F("Z2Mass_1", TString::Format("Z2Mass (%s)", _label1.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass12_1 = new TH1F("InvMass12_1", TString::Format("Primary Lepton Pair Invariant Mass (%s)", _label1.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34_1 = new TH1F("InvMass34_1", TString::Format("Secondary Lepton Pair Invariant Mass (%s)", _label1.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy_1 = new TH1F("LepEnergy_1", TString::Format("Lepton Energy (%s)", _label1.c_str()), 60, 0, 600);
  TH1F *LepPt_1 = new TH1F("LepPt_1", TString::Format("Lepton Transverse Momentum (%s)", _label1.c_str()), LepPt_binning.size()-1, &LepPt_binning[0]);
  SetTitles(InvMass4l_1, "m_{4l} [GeV]");
  SetTitles(Z2Mass_1, "m_{Z_{2}} [GeV]");
  SetTitles(InvMass12_1, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34_1, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy_1, "E_{l} [GeV]");
  SetTitles(LepPt_1, "p_{l,T} [GeV]");

  TH1F *InvMass4l_2 = new TH1F("InvMass4l_2", TString::Format("4-Lepton Invariant Mass (%s)", _label2.c_str()), InvMass_4l_binning.size()-1, &InvMass_4l_binning[0]);
  TH1F *Z2Mass_2 = new TH1F("Z2Mass_2", TString::Format("Z2Mass (%s)", _label2.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass12_2 = new TH1F("InvMass12_2", TString::Format("Primary Lepton Pair Invariant Mass (%s)", _label2.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34_2 = new TH1F("InvMass34_2", TString::Format("Secondary Lepton Pair Invariant Mass (%s)", _label2.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy_2 = new TH1F("LepEnergy_2", TString::Format("Lepton Energy (%s)", _label2.c_str()), 60, 0, 600);
  TH1F *LepPt_2 = new TH1F("LepPt_2", TString::Format("Lepton Transverse Momentum (%s)", _label2.c_str()), LepPt_binning.size()-1, &LepPt_binning[0]);
  SetTitles(InvMass4l_2, "m_{4l} [GeV]");
  SetTitles(Z2Mass_1, "m_{Z_{2}} [GeV]");
  SetTitles(InvMass12_2, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34_2, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy_2, "E_{l} [GeV]");
  SetTitles(LepPt_2, "p_{l,T} [GeV]");

  std::string cutstr1 = _isT1MC? "genWeight" : "";
  _ntuple1->Draw("Mass >> InvMass4l_1", cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%s_%s_Mass >> Z2Mass_1", l3.c_str(), l4.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sEnergy >> LepEnergy_1", l1.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sEnergy >>+LepEnergy_1", l2.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sEnergy >>+LepEnergy_1", l3.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sEnergy >>+LepEnergy_1", l4.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sPt >> LepPt_1", l1.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sPt >>+LepPt_1", l2.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sPt >>+LepPt_1", l3.c_str()), cutstr1.c_str(), "goff");
  _ntuple1->Draw(TString::Format("%sPt >>+LepPt_1", l4.c_str()), cutstr1.c_str(), "goff");

  std::string cutstr2 = _isT2MC? "genWeight" : "";
  _ntuple2->Draw("Mass >> InvMass4l_2", cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%s_%s_Mass >> Z2Mass_2", l3.c_str(), l4.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sEnergy >> LepEnergy_2", l1.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sEnergy >>+LepEnergy_2", l2.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sEnergy >>+LepEnergy_2", l3.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sEnergy >>+LepEnergy_2", l4.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sPt >> LepPt_2", l1.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sPt >>+LepPt_2", l2.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sPt >>+LepPt_2", l3.c_str()), cutstr2.c_str(), "goff");
  _ntuple2->Draw(TString::Format("%sPt >>+LepPt_2", l4.c_str()), cutstr2.c_str(), "goff");

  Long64_t nentries1 = _ntuple1->GetEntries();
  Long64_t nentries2 = _ntuple2->GetEntries();
  Long64_t nentries = nentries1 > nentries2? nentries1 : nentries2;

  ROOT::Math::PtEtaPhiEVector lp1, ln1, lp2, ln2;
  Float_t l1Pt, l1Eta, l1Phi, l1Energy;
  Float_t l2Pt, l2Eta, l2Phi, l2Energy;
  Float_t l3Pt, l3Eta, l3Phi, l3Energy;
  Float_t l4Pt, l4Eta, l4Phi, l4Energy;
  Int_t l1PdgId, l2PdgId, l3PdgId, l4PdgId;
  std::cout << std::endl << "Begin looping over " << nentries << " entries..." << std::endl;
  for (unsigned int i=0; i<nentries; i++){
    if (_ntuple1->GetEntry(i)){
      if (_channel == "eeee"){
        l1Pt = e1Pt1; l1Eta = e1Eta1; l1Phi = e1Phi1; l1Energy = e1Energy1;
        l2Pt = e2Pt1; l2Eta = e2Eta1; l2Phi = e2Phi1; l2Energy = e2Energy1;
        l3Pt = e3Pt1; l3Eta = e3Eta1; l3Phi = e3Phi1; l3Energy = e3Energy1;
        l4Pt = e4Pt1; l4Eta = e4Eta1; l4Phi = e4Phi1; l4Energy = e4Energy1;
        l1PdgId = e1PdgId1; l2PdgId = e2PdgId1;
        l3PdgId = e3PdgId1; l4PdgId = e4PdgId1;
      }
      else if (_channel == "eemm"){
        l1Pt = e1Pt1; l1Eta = e1Eta1; l1Phi = e1Phi1; l1Energy = e1Energy1;
        l2Pt = e2Pt1; l2Eta = e2Eta1; l2Phi = e2Phi1; l2Energy = e2Energy1;
        l3Pt = m1Pt1; l3Eta = m1Eta1; l3Phi = m1Phi1; l3Energy = m1Energy1;
        l4Pt = m2Pt1; l4Eta = m2Eta1; l4Phi = m2Phi1; l4Energy = m2Energy1;
        l1PdgId = e1PdgId1; l2PdgId = e2PdgId1;
        l3PdgId = m1PdgId1; l4PdgId = m2PdgId1;
      }
      else if (_channel == "mmmm"){
        l1Pt = m1Pt1; l1Eta = m1Eta1; l1Phi = m1Phi1; l1Energy = m1Energy1;
        l2Pt = m2Pt1; l2Eta = m2Eta1; l2Phi = m2Phi1; l2Energy = m2Energy1;
        l3Pt = m3Pt1; l3Eta = m3Eta1; l3Phi = m3Phi1; l3Energy = m3Energy1;
        l4Pt = m4Pt1; l4Eta = m4Eta1; l4Phi = m4Phi1; l4Energy = m4Energy1;
        l1PdgId = m1PdgId1; l2PdgId = m2PdgId1;
        l3PdgId = m3PdgId1; l4PdgId = m4PdgId1;
      }
      // Set electron pairs
      if (l1PdgId > 0){
        lp1 = ROOT::Math::PtEtaPhiEVector(l1Pt, l1Eta, l1Phi, l1Energy);
        ln1 = ROOT::Math::PtEtaPhiEVector(l2Pt, l2Eta, l2Phi, l2Energy);
      }
      else{
        ln1 = ROOT::Math::PtEtaPhiEVector(l1Pt, l1Eta, l1Phi, l1Energy);
        lp1 = ROOT::Math::PtEtaPhiEVector(l2Pt, l2Eta, l2Phi, l2Energy);
      }
      // Set muon pairs
      if (l3PdgId > 0){
        lp2 = ROOT::Math::PtEtaPhiEVector(l3Pt, l3Eta, l3Phi, l3Energy);
        ln2 = ROOT::Math::PtEtaPhiEVector(l4Pt, l4Eta, l4Phi, l4Energy);
      }
      else{
        ln2 = ROOT::Math::PtEtaPhiEVector(l3Pt, l3Eta, l3Phi, l3Energy);
        lp2 = ROOT::Math::PtEtaPhiEVector(l4Pt, l4Eta, l4Phi, l4Energy);
      }
      // Set primary/secondary pairs
      if (_channel == "eemm" && std::fabs((lp2+ln2).M() - Z_MASS) < std::fabs((lp1+ln1).M() - Z_MASS) ){
        ROOT::Math::PtEtaPhiEVector lp1_copy = lp1, ln1_copy = ln1;
        lp1 = lp2;
        ln1 = ln2;
        lp2 = lp1_copy;
        ln2 = ln1_copy;
      }
      // Sort by Pt
      std::vector<ROOT::Math::PtEtaPhiEVector> leptons{lp1, ln1, lp2, ln2};
      std::sort(leptons.begin(), leptons.end(), 
        [](const ROOT::Math::PtEtaPhiEVector &v1, const ROOT::Math::PtEtaPhiEVector &v2){
          return v1.Pt() > v2.Pt();
        }
      );

      Float_t weight = _isT1MC? genWeight1 : 1.;

      // Inv Mass (12, 34)
      InvMass12_1->Fill((lp1+ln1).M(), weight);
      InvMass34_1->Fill((lp2+ln2).M(), weight);

      // PolCosTheta
      //PolCosTheta_12->Fill(GetPolCosTheta(lp1, ln1), weight);
      //PolCosTheta_34->Fill(GetPolCosTheta(lp2, ln2), weight);
    }
    if (_ntuple2->GetEntry(i)){
      if (_channel == "eeee"){
        l1Pt = e1Pt2; l1Eta = e1Eta2; l1Phi = e1Phi2; l1Energy = e1Energy2;
        l2Pt = e2Pt2; l2Eta = e2Eta2; l2Phi = e2Phi2; l2Energy = e2Energy2;
        l3Pt = e3Pt2; l3Eta = e3Eta2; l3Phi = e3Phi2; l3Energy = e3Energy2;
        l4Pt = e4Pt2; l4Eta = e4Eta2; l4Phi = e4Phi2; l4Energy = e4Energy2;
        l1PdgId = e1PdgId2; l2PdgId = e2PdgId2;
        l3PdgId = e3PdgId2; l4PdgId = e4PdgId2;
      }
      else if (_channel == "eemm"){
        l1Pt = e1Pt2; l1Eta = e1Eta2; l1Phi = e1Phi2; l1Energy = e1Energy2;
        l2Pt = e2Pt2; l2Eta = e2Eta2; l2Phi = e2Phi2; l2Energy = e2Energy2;
        l3Pt = m1Pt2; l3Eta = m1Eta2; l3Phi = m1Phi2; l3Energy = m1Energy2;
        l4Pt = m2Pt2; l4Eta = m2Eta2; l4Phi = m2Phi2; l4Energy = m2Energy2;
        l1PdgId = e1PdgId2; l2PdgId = e2PdgId2;
        l3PdgId = m1PdgId2; l4PdgId = m2PdgId2;
      }
      else if (_channel == "mmmm"){
        l1Pt = m1Pt2; l1Eta = m1Eta2; l1Phi = m1Phi2; l1Energy = m1Energy2;
        l2Pt = m2Pt2; l2Eta = m2Eta2; l2Phi = m2Phi2; l2Energy = m2Energy2;
        l3Pt = m3Pt2; l3Eta = m3Eta2; l3Phi = m3Phi2; l3Energy = m3Energy2;
        l4Pt = m4Pt2; l4Eta = m4Eta2; l4Phi = m4Phi2; l4Energy = m4Energy2;
        l1PdgId = m1PdgId2; l2PdgId = m2PdgId2;
        l3PdgId = m3PdgId2; l4PdgId = m4PdgId2;
      }
      // Set electron pairs
      if (l1PdgId > 0){
        lp1 = ROOT::Math::PtEtaPhiEVector(l1Pt, l1Eta, l1Phi, l1Energy);
        ln1 = ROOT::Math::PtEtaPhiEVector(l2Pt, l2Eta, l2Phi, l2Energy);
      }
      else{
        ln1 = ROOT::Math::PtEtaPhiEVector(l1Pt, l1Eta, l1Phi, l1Energy);
        lp1 = ROOT::Math::PtEtaPhiEVector(l2Pt, l2Eta, l2Phi, l2Energy);
      }
      // Set muon pairs
      if (l3PdgId > 0){
        lp2 = ROOT::Math::PtEtaPhiEVector(l3Pt, l3Eta, l3Phi, l3Energy);
        ln2 = ROOT::Math::PtEtaPhiEVector(l4Pt, l4Eta, l4Phi, l4Energy);
      }
      else{
        ln2 = ROOT::Math::PtEtaPhiEVector(l3Pt, l3Eta, l3Phi, l3Energy);
        lp2 = ROOT::Math::PtEtaPhiEVector(l4Pt, l4Eta, l4Phi, l4Energy);
      }
      // Set primary/secondary pairs
      if ( std::abs((lp2+ln2).M() - Z_MASS) < std::abs((lp1+ln1).M() - Z_MASS) ){
        ROOT::Math::PtEtaPhiEVector lp1_copy = lp1, ln1_copy = ln1;
        lp1 = lp2;
        ln1 = ln2;
        lp2 = lp1_copy;
        ln2 = ln1_copy;
      }
      // Sort by Pt
      std::vector<ROOT::Math::PtEtaPhiEVector> leptons{lp1, ln1, lp2, ln2};
      std::sort(leptons.begin(), leptons.end(), 
        [](const ROOT::Math::PtEtaPhiEVector &v1, const ROOT::Math::PtEtaPhiEVector &v2){
          return v1.Pt() > v2.Pt();
        }
      );

      Float_t weight = _isT2MC? genWeight2 : 1.;

      // Inv Mass (12, 34)
      InvMass12_2->Fill((lp1+ln1).M(), weight);
      InvMass34_2->Fill((lp2+ln2).M(), weight);

      // PolCosTheta
      //PolCosTheta_12->Fill(GetPolCosTheta(lp1, ln1), weight);
      //PolCosTheta_34->Fill(GetPolCosTheta(lp2, ln2), weight);
    }
  }
  std::cout << "End looping." << std::endl;

	// Formatting
  InvMass4l_1->SetMinimum(0);
  Z2Mass_1->SetMinimum(0);
  InvMass12_1->SetMinimum(0);
  InvMass34_1->SetMinimum(0);
  LepEnergy_1->SetMinimum(0);
  LepPt_1->SetMinimum(0);

  InvMass4l_2->SetMinimum(0);
  Z2Mass_2->SetMinimum(0);
  InvMass12_2->SetMinimum(0);
  InvMass34_2->SetMinimum(0);
  LepEnergy_2->SetMinimum(0);
  LepPt_2->SetMinimum(0);

  // Scaling
  if (!_norm){
    if (_isT1MC){
      _ntuple1->GetEntry(0);
      std::cout << "summedWeights1=" << summedWeights1 << std::endl;
      float histScaling1 = _xsec1 * _lumi / summedWeights1;
      InvMass4l_1->Scale(histScaling1);
      Z2Mass_1->Scale(histScaling1);
      InvMass12_1->Scale(histScaling1);
      InvMass34_1->Scale(histScaling1);
      LepEnergy_1->Scale(histScaling1);
      LepPt_1->Scale(histScaling1);
    }
    if (_isT2MC){
      _ntuple2->GetEntry(0);
      std::cout << "summedWeights2=" << summedWeights2 << std::endl;
      float histScaling2 = _xsec2 * _lumi / summedWeights2;
      InvMass4l_2->Scale(histScaling2);
      Z2Mass_2->Scale(histScaling2);
      InvMass12_2->Scale(histScaling2);
      InvMass34_2->Scale(histScaling2);
      LepEnergy_2->Scale(histScaling2);
      LepPt_2->Scale(histScaling2);
    }
  }
  else{
    InvMass4l_1->Scale(1.0/InvMass4l_1->Integral());
    Z2Mass_1->Scale(1.0/Z2Mass_1->Integral());
    InvMass12_1->Scale(1.0/InvMass12_1->Integral());
    InvMass34_1->Scale(1.0/InvMass34_1->Integral());
    LepEnergy_1->Scale(1.0/LepEnergy_1->Integral()); 
    LepPt_1->Scale(1.0/LepPt_1->Integral());

    InvMass4l_2->Scale(1.0/InvMass4l_2->Integral());
    Z2Mass_2->Scale(1.0/Z2Mass_2->Integral());
    InvMass12_2->Scale(1.0/InvMass12_2->Integral());
    InvMass34_2->Scale(1.0/InvMass34_2->Integral());
    LepEnergy_2->Scale(1.0/LepEnergy_2->Integral()); 
    LepPt_2->Scale(1.0/LepPt_2->Integral());
  }
  
  // Drawing
  if (plot || _makePlots){
    c->cd();
    std::vector<std::string> labels = {_label1, _label2};

    auto Draw = [&c, &dirname, _filetype=_filetype](TH1* hist){
      hist->Draw("hist");
      c->SaveAs((dirname + "/" + hist->GetName() + _filetype).c_str());
    };
    auto DrawSame = [&c, &dirname, _filetype=_filetype, &labels](TH1* hist1, TH1* hist2){
      Float_t max = (hist1->GetMaximum() > hist2->GetMaximum())?
        hist1->GetMaximum() : hist2->GetMaximum();
      hist1->SetMaximum(max*1.2);
      hist2->SetMaximum(max*1.2);

      hist1->Draw("hist");
      hist2->SetLineColor(kRed);
      hist2->Draw("hist same");

      UInt_t idx = 0;
      TLegend *legend = c->BuildLegend();
      for (const TObject *obj : *legend->GetListOfPrimitives())
        ((TLegendEntry*)obj)->SetLabel(labels[idx++].c_str());
      std::string histname = hist1->GetName();
      c->SaveAs((dirname + "/" + histname.substr(0, histname.size()-2) + _filetype).c_str());
    };

    Draw(InvMass4l_1);
    Draw(Z2Mass_1);
    Draw(InvMass12_1);
    Draw(InvMass34_1);
    Draw(LepEnergy_1);
    Draw(LepPt_1);

    Draw(InvMass4l_2);
    Draw(Z2Mass_2);
    Draw(InvMass12_2);
    Draw(InvMass34_2);
    Draw(LepEnergy_2);
    Draw(LepPt_2);

    TLegend *legend;

    DrawSame(InvMass4l_1, InvMass4l_2);
    DrawSame(Z2Mass_1, Z2Mass_2);
    DrawSame(InvMass12_1, InvMass12_2);
    DrawSame(InvMass34_1, InvMass34_2);
    DrawSame(LepEnergy_1, LepEnergy_2);
    DrawSame(LepPt_1, LepPt_2);
  }

  // Writing
  histout->cd();
  histout->rmdir(_channel.c_str());
  TDirectory *subdir = histout->mkdir(_channel.c_str());
  if (subdir != nullptr) subdir->cd();

  InvMass4l_1->Write();
  Z2Mass_1->Write();
  InvMass12_1->Write();
  InvMass34_1->Write();
  LepEnergy_1->Write();
  LepPt_1->Write();

  InvMass4l_2->Write();
  Z2Mass_2->Write();
  InvMass12_2->Write();
  InvMass34_2->Write();
  LepEnergy_2->Write();
  LepPt_2->Write();

  // Ratio Plots
  if (ratioplot || _makeRatios){
  }

  c->Close();
  histout->Close();
}

#ifndef __CLING__
#include "interface/argparse.h"
int main(int nargs, char *argv[]){
	auto parser = argparse::ArgumentParser(nargs, argv)
		.formatter_class(argparse::HelpFormatter::ArgumentDefaults);

  parser.add_argument<bool>("--mc1").def("false")
    .help("if true, consider first input to be MC");
  parser.add_argument<bool>("--mc2").def("false")
    .help("if true, consider second input to be MC");
  parser.add_argument<double>("-l", "--lumi").def("7.561502251") //2022CD
    .help("data luminosity to scale MC in fb-1");
  parser.add_argument<double>("-x", "--xsec1").def("1390")
    .help("cross-section of given MC process in fb for first input");
  parser.add_argument<double>("-X", "--xsec2").def("1390")
    .help("cross-section of given MC process in fb for second input");
  parser.add_argument("-c", "--channel").choices("eeee,eemm,mmmm")
    .help("process only the specified channel, otherwise all");
  parser.add_argument<bool>("-n", "--norm").def("false")
    .help("if true, scale histograms to 1");
  parser.add_argument("-m", "--mode").def("UPDATE")
    .help("option for creating the output file");
  parser.add_argument<bool>("--noplots").def("false")
    .help("if true, don't save plots to 'plots/' directory");
  parser.add_argument("-t", "--filetype").def(".png")
    .help("filetype for output plots");

	parser.add_argument("name")
		.help("name for output file (and optional plots subdirectory)");
  parser.add_argument("label1").metavar("first-label")
		.help("name for first input slimmed ntuple, as seen in 'slimmed/*.root");
  parser.add_argument("label2").metavar("second-label")
		.help("name for second input slimmed ntuple, as seen in 'slimmed/*.root");

  auto args = parser.parse_args();

	for (std::string channel : {"eeee", "eemm", "mmmm"}){
    if (!args["channel"].is_none() && args["channel"] != channel)
      continue;
		CompLooper l(args["name"].c_str(), channel.c_str(), args["label1"].c_str(), args["label2"].c_str());
    if (!args["noplots"].is_true()) l.SetMakePlots();
    if (args["norm"].is_true()) l.SetNorm();
    l.SetMode(args["mode"]);
    l.SetPlotFiletype(args["filetype"]);
		
		if (args["mc1"].is_true()){
      l.SetMC1();
      l.SetLumi(args["lumi"]);
      l.SetXsec1(args["xsec1"]);
		}
		if (args["mc2"].is_true()){
      l.SetMC2();
      l.SetLumi(args["lumi"]);
      l.SetXsec2(args["xsec2"]);
		}
		l.Loop();
	}
	return 0;
  CompLooper l("comp_4e", "eeee", "2022MC", "2022Data");
  l.SetMC1();
  l.SetMakePlots();
  l.SetLumi(7.561502251);//fb-1
  l.SetXsec1(1.390*1000);//fb
  l.Loop();
  return 0;
}
#endif
