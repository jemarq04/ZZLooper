#include "interface/ZZLooperBase.h"
#include <iostream>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TLegendEntry.h"
#include "ROOT/RDataFrame.hxx"

//NOTE: Only accepts 2022 skimmed ntuples
class ZZLooper : public ZZLooperBase {
  public:
    ZZLooper(const char *name, const char *channel, const char *filename);
    ~ZZLooper();

    //Debug
    void SetNorm(bool norm=true){_norm = norm;}
    void SetDeduplicate(bool deduplicate=true){_deduplicate=deduplicate;}

    void SetMode(std::string mode){_mode = mode;}
    void SetMakePlots(bool val=true);
    void SetPlotFiletype(std::string ft=".png");
    void Loop();
  private:
    bool _makePlots = false;
    bool _norm = false, _deduplicate = false;
    std::string _filetype = ".png", _mode = "UPDATE";
};

ZZLooper::ZZLooper(const char *name, const char *channel, const char *filename)
    : ZZLooperBase(name, channel, filename){};

ZZLooper::~ZZLooper(){}

void ZZLooper::SetMakePlots(bool val){_makePlots = val;}
void ZZLooper::SetPlotFiletype(std::string ft){
  if (ft.size() == 0) return;
  if (ft[0] != '.') ft = "." + ft;
  _filetype = ft;
}

void ZZLooper::Loop(){
  Init();

  // Options
  bool plot = false; //override for _makePlots

  // Setup
  TFile *histout = new TFile(("out/histout_" + _name + ".root").c_str(), "update");
  if (!histout->IsOpen()) return;
  histout->cd();
	TCanvas *c = new TCanvas("c", "Histogram Canvas");
  c->cd();
  bool makeDirs = (plot || _makePlots);
	
  // Optionally creating plot directories
  std::string dirname = "plots/" + _name;
  if (makeDirs){
    DIR *plotdir = opendir("plots");
    if (!plotdir) mkdir("plots", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);

    plotdir = opendir(dirname.c_str());
    if (!plotdir) mkdir(dirname.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);
    MakeIndexFile(dirname);

    dirname += "/" + _channel;
    plotdir = opendir(dirname.c_str());
    if (!plotdir) mkdir(dirname.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);
    MakeIndexFile(dirname);
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
  std::vector<Double_t> InvMass4l_binning   = {100.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0};
  std::vector<Double_t> InvMass_pair_binning = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120.};
  std::vector<Double_t> LepPt_binning        = {0., 30., 50., 100., 200.};

  // Declaring Histograms
  TH1F *InvMass4l = new TH1F("InvMass4l", "4-Lepton Invariant Mass", InvMass4l_binning.size()-1, &InvMass4l_binning[0]);
  TH1F *InvMass12 = new TH1F("InvMass12", "Primary Lepton Pair Invariant Mass", InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34 = new TH1F("InvMass34", "Secondary Lepton Pair Invariant Mass", InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy = new TH1F("LepEnergy", "Lepton Energy", 60, 0, 600);
  TH1F *LepPt = new TH1F("LepPt", "Lepton Transverse Momentum", LepPt_binning.size()-1, &LepPt_binning[0]);
  TH1F *PolCosTheta12 = new TH1F("PolCosTheta12", "cos#theta_{12}*", 50, -1, 1);
  TH1F *PolCosTheta34 = new TH1F("PolCosTheta34", "cos#theta_{34}*", 50, -1, 1);
  SetTitles(InvMass4l, "m_{4l} [GeV]");
  SetTitles(InvMass12, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy, "E_{l} [GeV]");
  SetTitles(LepPt, "p_{l,T} [GeV]");


  Long64_t nentries = _ntuple->GetEntries();
  ROOT::Math::PtEtaPhiEVector lp1, ln1, lp2, ln2;
	Float_t Z1mass, Z2mass;
  Float_t l1Pt, l1Eta, l1Phi, l1Energy;
  Float_t l2Pt, l2Eta, l2Phi, l2Energy;
  Float_t l3Pt, l3Eta, l3Phi, l3Energy;
  Float_t l4Pt, l4Eta, l4Phi, l4Energy;
  Int_t l1PdgId, l2PdgId, l3PdgId, l4PdgId;
  std::vector<ULong64_t> evts;
  std::cout << std::endl << "Begin looping over " << nentries << " entries..." << std::endl;
  for (unsigned int i=0; i<nentries; i++){
    _ntuple->GetEntry(i);
    if (!_isMC && _deduplicate){
      if (std::find(evts.begin(), evts.end(), evt) != evts.end())
        continue;
      else evts.push_back(evt);
    }
    if (_channel == "eeee"){
			Z1mass = e1_e2_Mass; Z2mass = e3_e4_Mass;
      l1Pt = e1Pt; l1Eta = e1Eta; l1Phi = e1Phi; l1Energy = e1Energy;
      l2Pt = e2Pt; l2Eta = e2Eta; l2Phi = e2Phi; l2Energy = e2Energy;
      l3Pt = e3Pt; l3Eta = e3Eta; l3Phi = e3Phi; l3Energy = e3Energy;
      l4Pt = e4Pt; l4Eta = e4Eta; l4Phi = e4Phi; l4Energy = e4Energy;
      l1PdgId = e1PdgId; l2PdgId = e2PdgId;
      l3PdgId = e3PdgId; l4PdgId = e4PdgId;
    }
    else if (_channel == "eemm"){
			if (std::fabs(e1_e2_Mass - Z_MASS) < std::fabs(m1_m2_Mass - Z_MASS)){
				Z1mass = e1_e2_Mass; Z2mass = m1_m2_Mass;
				l1Pt = e1Pt; l1Eta = e1Eta; l1Phi = e1Phi; l1Energy = e1Energy;
				l2Pt = e2Pt; l2Eta = e2Eta; l2Phi = e2Phi; l2Energy = e2Energy;
				l3Pt = m1Pt; l3Eta = m1Eta; l3Phi = m1Phi; l3Energy = m1Energy;
				l4Pt = m2Pt; l4Eta = m2Eta; l4Phi = m2Phi; l4Energy = m2Energy;
				l1PdgId = e1PdgId; l2PdgId = e2PdgId;
				l3PdgId = m1PdgId; l4PdgId = m2PdgId;
			}
			else{
				Z1mass = m1_m2_Mass; Z2mass = e1_e2_Mass;
				l1Pt = m1Pt; l1Eta = m1Eta; l1Phi = m1Phi; l1Energy = m1Energy;
				l2Pt = m2Pt; l2Eta = m2Eta; l2Phi = m2Phi; l2Energy = m2Energy;
				l3Pt = e1Pt; l3Eta = e1Eta; l3Phi = e1Phi; l3Energy = e1Energy;
				l4Pt = e2Pt; l4Eta = e2Eta; l4Phi = e2Phi; l4Energy = e2Energy;
				l1PdgId = m1PdgId; l2PdgId = m2PdgId;
				l3PdgId = e1PdgId; l4PdgId = e2PdgId;
			}
    }
    else if (_channel == "mmmm"){
			Z1mass = m1_m2_Mass; Z2mass = m3_m4_Mass;
      l1Pt = m1Pt; l1Eta = m1Eta; l1Phi = m1Phi; l1Energy = m1Energy;
      l2Pt = m2Pt; l2Eta = m2Eta; l2Phi = m2Phi; l2Energy = m2Energy;
      l3Pt = m3Pt; l3Eta = m3Eta; l3Phi = m3Phi; l3Energy = m3Energy;
      l4Pt = m4Pt; l4Eta = m4Eta; l4Phi = m4Phi; l4Energy = m4Energy;
      l1PdgId = m1PdgId; l2PdgId = m2PdgId;
      l3PdgId = m3PdgId; l4PdgId = m4PdgId;
    }
    // Set primary pair
    if (l1PdgId > 0){
      lp1 = ROOT::Math::PtEtaPhiEVector(l1Pt, l1Eta, l1Phi, l1Energy);
      ln1 = ROOT::Math::PtEtaPhiEVector(l2Pt, l2Eta, l2Phi, l2Energy);
    }
    else{
      ln1 = ROOT::Math::PtEtaPhiEVector(l1Pt, l1Eta, l1Phi, l1Energy);
      lp1 = ROOT::Math::PtEtaPhiEVector(l2Pt, l2Eta, l2Phi, l2Energy);
    }
    // Set secondary pair
    if (l3PdgId > 0){
      lp2 = ROOT::Math::PtEtaPhiEVector(l3Pt, l3Eta, l3Phi, l3Energy);
      ln2 = ROOT::Math::PtEtaPhiEVector(l4Pt, l4Eta, l4Phi, l4Energy);
    }
    else{
      ln2 = ROOT::Math::PtEtaPhiEVector(l3Pt, l3Eta, l3Phi, l3Energy);
      lp2 = ROOT::Math::PtEtaPhiEVector(l4Pt, l4Eta, l4Phi, l4Energy);
    }
    // Sort by Pt
		/*
    std::vector<ROOT::Math::PtEtaPhiEVector> leptons{lp1, ln1, lp2, ln2};
    std::sort(leptons.begin(), leptons.end(), 
      [](const ROOT::Math::PtEtaPhiEVector &v1, const ROOT::Math::PtEtaPhiEVector &v2){
        return v1.Pt() > v2.Pt();
      }
    );
		*/

		if (Z1mass > 60.0 && Z1mass < 120.0 && Z2mass > 60.0 && Z2mass < 120.0){ 
			Float_t weight = _isMC? genWeight : 1.;

			InvMass4l->Fill(Mass, weight);
			InvMass12->Fill(Z1mass, weight);
			InvMass34->Fill(Z2mass, weight);

			LepEnergy->Fill(l1Energy, weight);
			LepEnergy->Fill(l2Energy, weight);
			LepEnergy->Fill(l3Energy, weight);
			LepEnergy->Fill(l4Energy, weight);

			LepPt->Fill(l1Pt, weight);
			LepPt->Fill(l2Pt, weight);
			LepPt->Fill(l3Pt, weight);
			LepPt->Fill(l4Pt, weight);

			PolCosTheta12->Fill(GetPolCosTheta(lp1, ln1), weight);
			PolCosTheta34->Fill(GetPolCosTheta(lp2, ln2), weight);
		}
  }
  std::cout << "End looping." << std::endl;

  // Formatting
  InvMass4l->SetMinimum(0);
  InvMass12->SetMinimum(0);
  InvMass34->SetMinimum(0);
  PolCosTheta12->SetMinimum(0);
  PolCosTheta34->SetMinimum(0);
  LepEnergy->SetMinimum(0);
  LepPt->SetMinimum(0);

  // Scaling
  if (_norm){
    InvMass4l->Scale(1.0/InvMass4l->Integral());
    InvMass12->Scale(1.0/InvMass12->Integral());
    InvMass34->Scale(1.0/InvMass34->Integral());
    PolCosTheta12->Scale(1.0/PolCosTheta12->Integral());
    PolCosTheta34->Scale(1.0/PolCosTheta34->Integral());
    LepEnergy->Scale(1.0/LepEnergy->Integral());
    LepPt->Scale(1.0/LepPt->Integral());
  }
  else if (_isMC){
    _ntuple->GetEntry(0);
    float histScaling = _kfac * _xsec * _lumi / summedWeights;
    InvMass4l->Scale(histScaling);
    InvMass12->Scale(histScaling);
    InvMass34->Scale(histScaling);
    PolCosTheta12->Scale(histScaling);
    PolCosTheta34->Scale(histScaling);
    LepEnergy->Scale(histScaling);
    LepPt->Scale(histScaling);
  }

  // Drawing
  if (plot || _makePlots){
    c->cd();
    auto Draw = [&c, &dirname, _filetype=_filetype](TH1* hist){
      hist->Draw("hist");
      c->SaveAs((dirname + "/" + hist->GetName() + _filetype).c_str());
    };

    Draw(InvMass4l);
    Draw(InvMass12);
    Draw(InvMass34);
    Draw(PolCosTheta12);
    Draw(PolCosTheta34);
    Draw(LepEnergy);
    Draw(LepPt);
  }

  // Writing
  histout->cd();
  histout->rmdir(_channel.c_str());
  TDirectory *subdir = histout->mkdir(_channel.c_str());
  if (subdir != nullptr) subdir->cd();

  InvMass4l->Write();
  InvMass12->Write();
  InvMass34->Write();
  PolCosTheta12->Write();
  PolCosTheta34->Write();
  LepEnergy->Write();
  LepPt->Write();
  
  c->Close();
  histout->Close();
}

#ifndef __CLING__
#include "interface/argparse.h"
int main(int nargs, char *argv[]){
	auto parser = argparse::ArgumentParser(nargs, argv)
		.formatter_class(argparse::HelpFormatter::ArgumentDefaults);

  parser.add_argument<bool>("--mc").def("false")
    .help("if true, add additional MC branches to slimmed ntuple");
  parser.add_argument<double>("-l", "--lumi").def("7.561502251") //2022CD
    .help("data luminosity to scale MC in fb-1");
  parser.add_argument<double>("-x", "--xsec").def("1390") //qqZZ
    .help("cross-section of given MC process in fb");
	parser.add_argument<double>("-k", "--kfac").def("1.2") //qqZZ kfac ~1.1 plus missing ggZZ signal
		.help("k-factor of given MC process");
  parser.add_argument("-c", "--channels")
    .help("process only the specified comma-separated channels, otherwise all. options: eeee, eemm, mmmm");
  parser.add_argument<bool>("-n", "--norm").def("false")
    .help("if true, scale histograms to 1");
  parser.add_argument<bool>("--deduplicate").def("false")
    .help("if true, skip events whose event number has already been processed");
  parser.add_argument<bool>("-r", "--recreate").def("false")
    .help("if true, the output file will be recreated for the given channel(s)");
  parser.add_argument<bool>("--noplots").def("false")
    .help("if true, don't save plots to 'plots/' directory");
  parser.add_argument("-t", "--filetype").def(".png")
    .help("filetype for output plots");
  parser.add_argument("label").help("name for input dataset, as seen in 'slimmed/*.root");

  auto args = parser.parse_args();

  std::vector<std::string> channels;
  if (args["channels"].is_none()) channels = {"eeee", "eemm", "mmmm"};
  else{
    std::stringstream ss(args["channels"].str());
    for (std::string channel; std::getline(ss, channel, ',');){
      if (channel != "eeee" && channel != "eemm" && channel != "mmmm")
        parser.error("invalid channel " + channel);
      channels.push_back(channel);
    }
  }

  std::string filename = (std::string)"slimmed/" + args["label"].str() + ".root";
  int index=0;
  for (std::string channel : channels){
    ZZLooper l(args["label"].c_str(), channel.c_str(), filename.c_str());
    if ((index++)==0 && args["recreate"].is_true()) l.SetMode("recreate");
    l.SetMakePlots(!args["noplots"].is_true());
    l.SetNorm(args["norm"]);
    l.SetDeduplicate(args["deduplicate"]);
    l.SetPlotFiletype(args["filetype"]);

    if (args["mc"].is_true()){
      l.SetMC();
      l.SetLumi(args["lumi"]);
      l.SetXsec(args["xsec"]);
      l.SetKfac(args["kfac"]);
    }
    l.Loop();
  }
	
  return 0;
}
#endif
