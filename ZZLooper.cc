#include "interface/ZZLooperBase.h"
#include <iostream>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TLegendEntry.h"
#include "ROOT/RDataFrame.hxx"

#ifndef NOSF
#include "correction.h"
#endif

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

    double GetScaleFactor();

    void Loop(bool applyScaleFacs=false);
  private:
    bool _makePlots = false;
    bool _norm = false, _deduplicate = false;
    std::string _filetype = ".png", _mode = "UPDATE";
#ifndef NOSF
    std::unique_ptr<correction::CorrectionSet> _pileupSF, _eIdSF, _eRecoSF, _mIdSF;
#endif
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

double ZZLooper::GetScaleFactor(){
  double weight = 1.0;

#ifndef NOSF
  auto GetEleRecoSFName = [](float pt){
    if (pt < 20) return "RecoBelow20";
    else if (pt < 75) return "Reco20to75";
    else return "RecoAbove75";
  };
  if (_channel == "eeee"){
    if (_eIdSF != nullptr){}
    if (_eRecoSF != nullptr){
      const auto recoref = (*_eRecoSF->begin()).second;
      if (e1Pt > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt), e1Eta, e1Pt});
      if (e2Pt > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt), e2Eta, e2Pt});
      if (e3Pt > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e3Pt), e3Eta, e3Pt});
      if (e4Pt > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e4Pt), e4Eta, e4Pt});
    }
  }
  else if (_channel == "eemm"){
    if (_eIdSF != nullptr){}
    if (_eRecoSF != nullptr){
      const auto recoref = (*_eRecoSF->begin()).second;
      if (e1Pt > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt), e1Eta, e1Pt});
      if (e2Pt > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt), e2Eta, e2Pt});
    }
    if (_mIdSF != nullptr){
      if (m1Pt > 15) weight *= _mIdSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta), m1Pt, "nominal"});
      if (m2Pt > 15) weight *= _mIdSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta), m2Pt, "nominal"});
    }
  }
  else if (_channel == "mmmm"){
    if (_mIdSF != nullptr){
      if (m1Pt > 15) weight *= _mIdSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta), m1Pt, "nominal"});
      if (m2Pt > 15) weight *= _mIdSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta), m2Pt, "nominal"});
      if (m3Pt > 15) weight *= _mIdSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({std::fabs(m3Eta), m3Pt, "nominal"});
      if (m4Pt > 15) weight *= _mIdSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({std::fabs(m4Eta), m4Pt, "nominal"});
    }
  }
  if (_pileupSF != nullptr)
    weight *= (*_pileupSF->begin()).second->evaluate({nTruePU, "nominal"});
#endif
  return weight;
}

void ZZLooper::Loop(bool applyScaleFacs){
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

#ifndef NOSF
  // Setup correction sets
  if (applyScaleFacs && _isMC){
    std::string basename = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";
    _pileupSF = correction::CorrectionSet::from_file(basename + "/LUM/2022_Summer22/puWeights.json.gz");
    //_eIdSF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22/electron.json.gz");
    _eRecoSF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22/electron.json.gz");
    _mIdSF = correction::CorrectionSet::from_file(basename + "/MUO/2022_Summer22/muon_Z.json.gz");
  }
#endif

  // Binning
  std::vector<Double_t> InvMass4l_binning   = {100.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0};
  std::vector<Double_t> InvMass_pair_binning = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120.};
  std::vector<Double_t> LepPt_binning        = {0., 30., 50., 100., 200.};
  
  // Declaring Histograms
  gROOT->cd();
  std::vector<std::string> histnames = {
    "InvMass4l", "InvMass12", "InvMass34", "LepEnergy", "LepPt", "LepEta",
    "PolCosTheta12", "PolCosTheta34",
    "LepZZIsoVal", "LepZZIsoPass", "LepSIP3D"
  };

  TH1F *InvMass4l = new TH1F("InvMass4l", "4-Lepton Invariant Mass", InvMass4l_binning.size()-1, &InvMass4l_binning[0]);
  TH1F *InvMass12 = new TH1F("InvMass12", "Primary Lepton Pair Invariant Mass", InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34 = new TH1F("InvMass34", "Secondary Lepton Pair Invariant Mass", InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy = new TH1F("LepEnergy", "Lepton Energy", 60, 0, 600);
  TH1F *LepPt = new TH1F("LepPt", "Lepton Transverse Momentum", LepPt_binning.size()-1, &LepPt_binning[0]);
  TH1F *LepEta = new TH1F("LepEta", "Lepton Pseudorapidity", 20, -5, 5);
  TH1F *PolCosTheta12 = new TH1F("PolCosTheta12", "cos#theta_{12}*", 50, -1, 1);
  TH1F *PolCosTheta34 = new TH1F("PolCosTheta34", "cos#theta_{34}*", 50, -1, 1);
  TH1F *LepZZIsoVal = new TH1F("LepZZIsoVal", "Isolation Value", 10, 0, 0.5);
  TH1F *LepZZIsoPass = new TH1F("LepZZIsoPass", "Isolation Cut", 2, 0, 2);
  TH1F *LepSIP3D = new TH1F("LepSIP3D", "Sip3D", 10, 0, 5);
  SetTitles(InvMass4l, "m_{4l} [GeV]");
  SetTitles(InvMass12, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy, "E_{l} [GeV]");
  SetTitles(LepPt, "p_{l,T} [GeV]");
  SetTitles(LepEta, "#eta");
  SetTitles(PolCosTheta12, "cos#theta_{12}*");
  SetTitles(PolCosTheta34, "cos#theta_{34}*");
  SetTitles(LepZZIsoVal, "Isolation");
  SetTitles(LepZZIsoPass, "Pass");
  SetTitles(LepSIP3D, "SIP3D");


  Long64_t nentries = _ntuple->GetEntries();
  ROOT::Math::PtEtaPhiEVector lp1, ln1, lp2, ln2;
  Float_t Z1mass, Z2mass;
  Float_t l1Pt, l1Eta, l1Phi, l1Energy;
  Float_t l2Pt, l2Eta, l2Phi, l2Energy;
  Float_t l3Pt, l3Eta, l3Phi, l3Energy;
  Float_t l4Pt, l4Eta, l4Phi, l4Energy;
  Int_t l1PdgId, l2PdgId, l3PdgId, l4PdgId;
  Float_t l1ZZIso, l2ZZIso, l3ZZIso, l4ZZIso;
  Bool_t l1ZZIsoPass, l2ZZIsoPass, l3ZZIsoPass, l4ZZIsoPass;
  Float_t l1SIP3D, l2SIP3D, l3SIP3D, l4SIP3D;
  std::vector<ULong64_t> evts;
  std::cout << std::endl << "Begin looping over " << nentries << " entries..." << std::endl;
  Int_t num_eemm = 0, num_mmee = 0;
  for (unsigned int i=0; i<nentries; i++){
    Bool_t mmee = false;
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
      l1ZZIso = e1ZZIso; l2ZZIso = e2ZZIso;
      l3ZZIso = e3ZZIso; l4ZZIso = e4ZZIso;
      l1ZZIsoPass = e1ZZIsoPass; l2ZZIsoPass = e2ZZIsoPass;
      l3ZZIsoPass = e3ZZIsoPass; l4ZZIsoPass = e4ZZIsoPass;
      l1SIP3D = e1SIP3D; l2SIP3D = e2SIP3D;
      l3SIP3D = e3SIP3D; l4SIP3D = e4SIP3D;
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
        l1ZZIso = e1ZZIso; l2ZZIso = e2ZZIso;
        l3ZZIso = m1ZZIso; l4ZZIso = m2ZZIso;
        l1ZZIsoPass = e1ZZIsoPass; l2ZZIsoPass = e2ZZIsoPass;
        l3ZZIsoPass = m1ZZIsoPass; l4ZZIsoPass = m2ZZIsoPass;
        l1SIP3D = e1SIP3D; l2SIP3D = e2SIP3D;
        l3SIP3D = m1SIP3D; l4SIP3D = m2SIP3D;
      }
      else{
        mmee = true;
        Z1mass = m1_m2_Mass; Z2mass = e1_e2_Mass;
        l1Pt = m1Pt; l1Eta = m1Eta; l1Phi = m1Phi; l1Energy = m1Energy;
        l2Pt = m2Pt; l2Eta = m2Eta; l2Phi = m2Phi; l2Energy = m2Energy;
        l3Pt = e1Pt; l3Eta = e1Eta; l3Phi = e1Phi; l3Energy = e1Energy;
        l4Pt = e2Pt; l4Eta = e2Eta; l4Phi = e2Phi; l4Energy = e2Energy;
        l1PdgId = m1PdgId; l2PdgId = m2PdgId;
        l3PdgId = e1PdgId; l4PdgId = e2PdgId;
        l1ZZIso = m1ZZIso; l2ZZIso = m2ZZIso;
        l3ZZIso = e1ZZIso; l4ZZIso = e2ZZIso;
        l1ZZIsoPass = m1ZZIsoPass; l2ZZIsoPass = m2ZZIsoPass;
        l3ZZIsoPass = e1ZZIsoPass; l4ZZIsoPass = e2ZZIsoPass;
        l1SIP3D = m1SIP3D; l2SIP3D = m2SIP3D;
        l3SIP3D = e1SIP3D; l4SIP3D = e2SIP3D;
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
      l1ZZIso = m1ZZIso; l2ZZIso = m2ZZIso;
      l3ZZIso = m3ZZIso; l4ZZIso = m4ZZIso;
      l1ZZIsoPass = m1ZZIsoPass; l2ZZIsoPass = m2ZZIsoPass;
      l3ZZIsoPass = m3ZZIsoPass; l4ZZIsoPass = m4ZZIsoPass;
      l1SIP3D = m1SIP3D; l2SIP3D = m2SIP3D;
      l3SIP3D = m3SIP3D; l4SIP3D = m4SIP3D;
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
    std::sort(leptons.begin(), leptons.end(), 
    std::vector<ROOT::Math::PtEtaPhiEVector> leptons{lp1, ln1, lp2, ln2};
      [](const ROOT::Math::PtEtaPhiEVector &v1, const ROOT::Math::PtEtaPhiEVector &v2){
        return v1.Pt() > v2.Pt();
      }
    );
    */

    bool ZZIsoPass = l1ZZIsoPass && l2ZZIsoPass && l3ZZIsoPass && l4ZZIsoPass;
    if (Z1mass > 60.0 && Z1mass < 120.0 && Z2mass > 60.0 && Z2mass < 120.0 && ZZIsoPass){ 
      Float_t weight = _isMC? genWeight : 1.;
      if (applyScaleFacs && _isMC) weight *= GetScaleFactor();

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

      LepEta->Fill(l1Eta, weight);
      LepEta->Fill(l2Eta, weight);
      LepEta->Fill(l3Eta, weight);
      LepEta->Fill(l4Eta, weight);

      PolCosTheta12->Fill(GetPolCosTheta(lp1, ln1), weight);
      PolCosTheta34->Fill(GetPolCosTheta(lp2, ln2), weight);

      LepZZIsoVal->Fill(l1ZZIso, weight);
      LepZZIsoVal->Fill(l2ZZIso, weight);
      LepZZIsoVal->Fill(l3ZZIso, weight);
      LepZZIsoVal->Fill(l4ZZIso, weight);

      LepZZIsoPass->Fill(l1ZZIsoPass, weight);
      LepZZIsoPass->Fill(l2ZZIsoPass, weight);
      LepZZIsoPass->Fill(l3ZZIsoPass, weight);
      LepZZIsoPass->Fill(l4ZZIsoPass, weight);

      LepSIP3D->Fill(l1SIP3D, weight);
      LepSIP3D->Fill(l2SIP3D, weight);
      LepSIP3D->Fill(l3SIP3D, weight);
      LepSIP3D->Fill(l4SIP3D, weight);

      if (_channel == "eemm"){
        if (mmee) num_mmee++;
        else num_eemm++;
      }
    }
  }
  std::cout << "End looping." << std::endl;
  
  // Writing README
  if (plot || _makePlots){
    std::ofstream readme((dirname + "/README").c_str());
    if (readme.is_open()){
      readme << "Events: " << InvMass4l->GetEntries() << std::endl;
      if (_channel == "eemm"){
        readme << "  - eemm: " << num_eemm << std::endl;
        readme << "  - mmee: " << num_mmee << std::endl;
      }
      if (_isMC && !_norm){
        float histScaling = _kfac * _xsec * _lumi / summedWeights;
        readme << "Scaling: " << histScaling << std::endl;
        readme << "  - Cross-section: " << _xsec << " fb" << std::endl;
        readme << "  - Luminosity: " << _lumi << " fb-1" << std::endl;
        readme << "  - Summed Weights: " << summedWeights << std::endl;
        readme << "  - k-factor: " << _kfac << std::endl;
      }
      if (applyScaleFacs) readme << "Efficiency SFs applied." << std::endl;
      readme.close();
    }
  }
  
  // Creating subdirs in histfile
  histout->cd();
  histout->rmdir(_channel.c_str());
  TDirectory *subdir = histout->mkdir(_channel.c_str());
  if (subdir != nullptr) subdir->cd();
  
  // Creating MC scale factors
  _ntuple->GetEntry(0);
  float histScaling = _kfac * _xsec * _lumi / summedWeights;
  
  // Defining draw functions
  auto Draw = [&c, &dirname, _filetype=_filetype](TH1* hist){
    hist->Draw("hist");
    c->SaveAs((dirname + "/" + hist->GetName() + _filetype).c_str());
  };
  
  // Finalizing
  c->cd();
  for (auto histname : histnames){
    TH1F *hist = (TH1F*)gROOT->Get(histname.c_str());
    if (hist == nullptr){
      std::cout << "Couldn't find hist. Skipping " << histname << std::endl;
      break;
    }
    // Formatting
    hist->SetMinimum(0);

    // Scaling
    if (_norm)
      hist->Scale(1.0/hist->Integral());
    else if (_isMC)
      hist->Scale(histScaling);

    // Drawing
    if (plot || _makePlots)
      Draw(hist);

    // Writing
    hist->Write();

    delete hist;
  }

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
  parser.add_argument<double>("-l", "--lumi").def("37.846170084") //2022
    .help("data luminosity to scale MC in fb-1");
  parser.add_argument<double>("-x", "--xsec").def("1390") //qqZZ
    .help("cross-section of given MC process in fb");
  parser.add_argument<double>("-k", "--kfac").def("1.2") //qqZZ kfac ~1.1 plus missing ggZZ signal
    .help("k-factor of given MC process");
#ifndef NOSF
  parser.add_argument<bool>("--sf").def("false")
    .help("if true, apply correctionlib scale factors");
#endif
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
    if ((index++)==0 && args["recreate"]) l.SetMode("recreate");
    l.SetMakePlots(!args["noplots"]);
    l.SetNorm(args["norm"]);
    l.SetDeduplicate(args["deduplicate"]);
    l.SetPlotFiletype(args["filetype"]);

    if (args["mc"].is_true()){
      l.SetMC();
      l.SetLumi(args["lumi"]);
      l.SetXsec(args["xsec"]);
      l.SetKfac(args["kfac"]);
    }
#ifndef NOSF
    l.Loop(args["sf"]);
#else
		l.Loop();
#endif
  }
  
  return 0;
}
#endif
