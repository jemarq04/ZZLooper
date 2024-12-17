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

#ifndef NOSF
#include "correction.h"
#endif

//NOTE: Only accepts 2022 skimmed ntuples
class CompLooper : public CompLooperBase {
  private:
    enum Ntuple {First, Second, EE};
  public:
    CompLooper(const char *name, const char *channel, const char *label1, const char *label2);
    ~CompLooper();

    //Debug
    void SetDebug(bool debug=true){_debug = debug;}
    void SetNorm(bool norm=true){_norm = norm;}
    void SetDeduplicate(bool deduplicate=true){_deduplicate=deduplicate;}

    void SetMode(std::string mode){_mode = mode;}
    void SetMakePlots(bool val=true){_makePlots = val;}
    void SetMakeRatios(bool val=true){_makeRatios = val;}
    void SetPlotFiletype(std::string ft=".png");

    double GetScaleFactor(Ntuple ntuple);

    void Loop(bool applyScaleFacs=false);
  private:
    std::string FindFile(const char *label);
    bool _debug = false;
    bool _makePlots = false, _makeRatios = false;
    bool _norm = false, _deduplicate = false;
    std::string _label1, _label2;
    std::string _filetype = ".png", _mode = "UPDATE";
#ifndef NOSF
    std::unique_ptr<correction::CorrectionSet> _pileupSF, _eIdSF, _eRecoSF, _mIdSF;
    std::unique_ptr<correction::CorrectionSet> _pileupEESF, _eIdEESF, _eRecoEESF, _mIdEESF;
#endif
};

CompLooper::CompLooper(const char *name, const char *channel, const char *label1, const char *label2)
    : _label1(label1), _label2(label2),
      CompLooperBase(name, channel, 
          ((std::string)"slimmed/" + label1 + ".root").c_str(), 
          ((std::string)"slimmed/" + label2 + ".root").c_str()
      ){
};

CompLooper::~CompLooper(){}

void CompLooper::SetPlotFiletype(std::string ft){
  if (ft.size() == 0) return;
  if (ft[0] != '.') ft = "." + ft;
  _filetype = ft;
}

double CompLooper::GetScaleFactor(Ntuple ntuple){
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
      const auto recorefEE = (*_eRecoEESF->begin()).second;
      if (ntuple == Ntuple::First){
        if (e1Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt1), e1Eta1, e1Pt1});
        if (e2Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt1), e2Eta1, e2Pt1});
        if (e3Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e3Pt1), e3Eta1, e3Pt1});
        if (e4Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e4Pt1), e4Eta1, e4Pt1});
      }
      else if (ntuple == Ntuple::Second){
        if (e1Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt2), e1Eta2, e1Pt2});
        if (e2Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt2), e2Eta2, e2Pt2});
        if (e3Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e3Pt2), e3Eta2, e3Pt2});
        if (e4Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e4Pt2), e4Eta2, e4Pt2});
      }
      else{
        if (e1PtEE > 10) weight *= recorefEE->evaluate({"2022Re-recoE+PromptFG", "sf", GetEleRecoSFName(e1PtEE), e1EtaEE, e1PtEE});
        if (e2PtEE > 10) weight *= recorefEE->evaluate({"2022Re-recoE+PromptFG", "sf", GetEleRecoSFName(e2PtEE), e2EtaEE, e2PtEE});
        if (e3PtEE > 10) weight *= recorefEE->evaluate({"2022Re-recoE+PromptFG", "sf", GetEleRecoSFName(e3PtEE), e3EtaEE, e3PtEE});
        if (e4PtEE > 10) weight *= recorefEE->evaluate({"2022Re-recoE+PromptFG", "sf", GetEleRecoSFName(e4PtEE), e4EtaEE, e4PtEE});
      }
    }
  }
  else if (_channel == "eemm"){
    if (_eIdSF != nullptr){}
    if (_eRecoSF != nullptr){
      const auto recoref = (*_eRecoSF->begin()).second;
      const auto recorefEE = (*_eRecoEESF->begin()).second;
      if (ntuple == Ntuple::First){
        if (e1Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt1), e1Eta1, e1Pt1});
        if (e2Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt1), e2Eta1, e2Pt1});
      }
      else if (ntuple == Ntuple::Second){
        if (e1Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt2), e1Eta2, e1Pt2});
        if (e2Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt2), e2Eta2, e2Pt2});
      }
      else{
        if (e1PtEE > 10) weight *= recorefEE->evaluate({"2022Re-recoE+PromptFG", "sf", GetEleRecoSFName(e1PtEE), e1EtaEE, e1PtEE});
        if (e2PtEE > 10) weight *= recorefEE->evaluate({"2022Re-recoE+PromptFG", "sf", GetEleRecoSFName(e2PtEE), e2EtaEE, e2PtEE});
      }
    }
    if (_mIdSF != nullptr){
      if (ntuple == Ntuple::First){
        if (m1Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta1), m1Pt1, "nominal"});
        if (m2Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta1), m2Pt1, "nominal"});
      }
      else if (ntuple == Ntuple::Second){
        if (m1Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta2), m1Pt2, "nominal"});
        if (m2Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta2), m2Pt2, "nominal"});
      }
      else{
        if (m1PtEE > 15) weight *= _mIdEESF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1EtaEE), m1PtEE, "nominal"});
        if (m2PtEE > 15) weight *= _mIdEESF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2EtaEE), m2PtEE, "nominal"});
      }
    }
  }
  else if (_channel == "mmmm"){
    if (_mIdSF != nullptr){
      if (ntuple == Ntuple::First){
        if (m1Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta1), m1Pt1, "nominal"});
        if (m2Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta1), m2Pt1, "nominal"});
        if (m3Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m3Eta1), m3Pt1, "nominal"});
        if (m4Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m4Eta1), m4Pt1, "nominal"});
      }
      else if (ntuple == Ntuple::Second){
        if (m1Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta2), m1Pt2, "nominal"});
        if (m2Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta2), m2Pt2, "nominal"});
        if (m3Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m3Eta2), m3Pt2, "nominal"});
        if (m4Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m4Eta2), m4Pt2, "nominal"});
      }
      else{
        if (m1PtEE > 15) weight *= _mIdEESF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1EtaEE), m1PtEE, "nominal"});
        if (m2PtEE > 15) weight *= _mIdEESF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2EtaEE), m2PtEE, "nominal"});
        if (m3PtEE > 15) weight *= _mIdEESF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m3EtaEE), m3PtEE, "nominal"});
        if (m4PtEE > 15) weight *= _mIdEESF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m4EtaEE), m4PtEE, "nominal"});
      }
    }
  }
  if (_pileupSF != nullptr){
    if (ntuple == Ntuple::First)
      weight *= (*_pileupSF->begin()).second->evaluate({nTruePU1, "nominal"});
    else if (ntuple == Ntuple::Second)
      weight *= (*_pileupSF->begin()).second->evaluate({nTruePU2, "nominal"});
    else
      weight *= (*_pileupEESF->begin()).second->evaluate({nTruePUEE, "nominal"});
  }
#endif
  return weight;
}

void CompLooper::Loop(bool applyScaleFacs){
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

    dirname += "/" + _channel;
    plotdir = opendir(dirname.c_str());
    if (!plotdir) mkdir(dirname.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);
    MakeIndexFile(dirname);

    plotdir = opendir((dirname + "/" + _label1).c_str());
    if (!plotdir) mkdir((dirname + "/" + _label1).c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);
    MakeIndexFile(dirname + "/" + _label1);

    plotdir = opendir((dirname + "/" + _label2).c_str());
    if (!plotdir) mkdir((dirname + "/" + _label2).c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    else closedir(plotdir);
    MakeIndexFile(dirname + "/" + _label2);
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
  if (applyScaleFacs && (_isT1MC || _isT2MC)){
    std::string basename = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";
    _pileupSF = correction::CorrectionSet::from_file(basename + "/LUM/2022_Summer22/puWeights.json.gz");
    _pileupEESF = correction::CorrectionSet::from_file(basename + "/LUM/2022_Summer22EE/puWeights.json.gz");
    //_eIdSF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22/electron.json.gz");
    //_eIdEESF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22EE/electron.json.gz");
    _eRecoSF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22/electron.json.gz");
    _eRecoEESF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22EE/electron.json.gz");
    _mIdSF = correction::CorrectionSet::from_file(basename + "/MUO/2022_Summer22/muon_Z.json.gz");
    _mIdEESF = correction::CorrectionSet::from_file(basename + "/MUO/2022_Summer22EE/muon_Z.json.gz");
    //if (_pileupSF != nullptr) std::cout << "Applying pileupSF..." << std::endl;
    //if (_eIdSF != nullptr) std::cout << "Applying eIdSF..." << std::endl;
    //if (_eRecoSF != nullptr) std::cout << "Applying eRecoSF..." << std::endl;
    //if (_mIdSF != nullptr) std::cout << "Applying mIdSF..." << std::endl;
  }
#endif

  // Binning
  std::vector<Double_t> InvMass_4l_binning   = {100.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0};
  std::vector<Double_t> InvMass_pair_binning = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120.};
  std::vector<Double_t> LepPt_binning        = {0., 5., 30., 50., 100., 200.};

  // Declaring Histograms
  std::vector<std::string> histnames = {
    "InvMass4l", "InvMass12", "InvMass34", "LepEnergy", "LepPt", "LepEta",
    "PolCosTheta12", "PolCosTheta34",
    "LepZZIsoVal", "LepZZIsoPass", "LepSIP3D"
  };

  gROOT->cd();
  TH1F *InvMass4l_1 = new TH1F("InvMass4l_1", TString::Format("4-Lepton Invariant Mass (%s)", _label1.c_str()), InvMass_4l_binning.size()-1, &InvMass_4l_binning[0]);
  TH1F *InvMass12_1 = new TH1F("InvMass12_1", TString::Format("Primary Lepton Pair Invariant Mass (%s)", _label1.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34_1 = new TH1F("InvMass34_1", TString::Format("Secondary Lepton Pair Invariant Mass (%s)", _label1.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy_1 = new TH1F("LepEnergy_1", TString::Format("Lepton Energy (%s)", _label1.c_str()), 60, 0, 600);
  TH1F *LepPt_1 = new TH1F("LepPt_1", TString::Format("Lepton Transverse Momentum (%s)", _label1.c_str()), LepPt_binning.size()-1, &LepPt_binning[0]);
  TH1F *LepEta_1 = new TH1F("LepEta_1", TString::Format("Lepton Psuedorapidity (%s)", _label1.c_str()), 20, -5, 5);
  TH1F *PolCosTheta12_1 = new TH1F("PolCosTheta12_1", TString::Format("cos#theta_{12}* (%s)", _label1.c_str()), 50, -1, 1);
  TH1F *PolCosTheta34_1 = new TH1F("PolCosTheta34_1", TString::Format("cos#theta_{34}* (%s)", _label1.c_str()), 50, -1, 1);
  TH1F *LepZZIsoVal_1 = new TH1F("LepZZIsoVal_1", TString::Format("Isolation Value (%s)", _label1.c_str()), 10, 0, 0.5);
  TH1F *LepZZIsoPass_1 = new TH1F("LepZZIsoPass_1", TString::Format("Isolation Cut (%s)", _label1.c_str()), 2, 0, 2);
  TH1F *LepSIP3D_1 = new TH1F("LepSIP3D_1", TString::Format("Sip3D (%s)", _label1.c_str()), 10, 0, 5);
  SetTitles(InvMass4l_1, "m_{4l} [GeV]");
  SetTitles(InvMass12_1, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34_1, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy_1, "E_{l} [GeV]");
  SetTitles(LepPt_1, "p_{l,T} [GeV]");
  SetTitles(LepEta_1, "#eta");
  SetTitles(PolCosTheta12_1, "cos#theta_{12}*");
  SetTitles(PolCosTheta34_1, "cos#theta_{34}*");
  SetTitles(LepZZIsoVal_1, "Isolation");
  SetTitles(LepZZIsoPass_1, "Pass");
  SetTitles(LepSIP3D_1, "SIP3D");

  TH1F *InvMass4l_2 = new TH1F("InvMass4l_2", TString::Format("4-Lepton Invariant Mass (%s)", _label2.c_str()), InvMass_4l_binning.size()-1, &InvMass_4l_binning[0]);
  TH1F *InvMass12_2 = new TH1F("InvMass12_2", TString::Format("Primary Lepton Pair Invariant Mass (%s)", _label2.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34_2 = new TH1F("InvMass34_2", TString::Format("Secondary Lepton Pair Invariant Mass (%s)", _label2.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy_2 = new TH1F("LepEnergy_2", TString::Format("Lepton Energy (%s)", _label2.c_str()), 60, 0, 600);
  TH1F *LepPt_2 = new TH1F("LepPt_2", TString::Format("Lepton Transverse Momentum (%s)", _label2.c_str()), LepPt_binning.size()-1, &LepPt_binning[0]);
  TH1F *LepEta_2 = new TH1F("LepEta_2", TString::Format("Lepton Psuedorapidity (%s)", _label2.c_str()), 20, -5, 5);
  TH1F *PolCosTheta12_2 = new TH1F("PolCosTheta12_2", TString::Format("cos#theta_{12}* (%s)", _label2.c_str()), 50, -1, 1);
  TH1F *PolCosTheta34_2 = new TH1F("PolCosTheta34_2", TString::Format("cos#theta_{34}* (%s)", _label2.c_str()), 50, -1, 1);
  TH1F *LepZZIsoVal_2 = new TH1F("LepZZIsoVal_2", TString::Format("Isolation Value (%s)", _label2.c_str()), 10, 0, 0.5);
  TH1F *LepZZIsoPass_2 = new TH1F("LepZZIsoPass_2", TString::Format("Isolation Cut (%s)", _label2.c_str()), 2, 0, 2);
  TH1F *LepSIP3D_2 = new TH1F("LepSIP3D_2", TString::Format("Sip3D (%s)", _label2.c_str()), 10, 0, 5);
  SetTitles(InvMass4l_2, "m_{4l} [GeV]");
  SetTitles(InvMass12_2, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34_2, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy_2, "E_{l} [GeV]");
  SetTitles(LepPt_2, "p_{l,T} [GeV]");
  SetTitles(LepEta_2, "#eta");
  SetTitles(PolCosTheta12_2, "cos#theta_{12}*");
  SetTitles(PolCosTheta34_2, "cos#theta_{34}*");
  SetTitles(LepZZIsoVal_2, "Isolation");
  SetTitles(LepZZIsoPass_2, "Pass");
  SetTitles(LepSIP3D_2, "SIP3D");

  Long64_t nentries1 = _ntuple1->GetEntries();
  Long64_t nentries2 = _ntuple2->GetEntries();
  if (_doEE){
    Long64_t nentriesEE= _ntupleEE->GetEntries();
    nentries1 = nentries1 > nentriesEE? nentries1 : nentriesEE;
  }
  Long64_t nentries = nentries1 > nentries2? nentries1 : nentries2;

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
  std::vector<ULong64_t> evts1={}, evts2={};
  std::cout << std::endl << "Begin looping over " << nentries << " entries..." << std::endl;
  for (unsigned int i=0; i<nentries; i++){
    if (_ntuple1->GetEntry(i) && std::find(evts1.begin(), evts1.end(), evt1) == evts1.end()){
      if (!_isT1MC && _deduplicate) evts1.push_back(evt1);
      if (_channel == "eeee"){
        Z1mass = e1_e2_Mass1; Z2mass = e3_e4_Mass1;
        l1Pt = e1Pt1; l1Eta = e1Eta1; l1Phi = e1Phi1; l1Energy = e1Energy1;
        l2Pt = e2Pt1; l2Eta = e2Eta1; l2Phi = e2Phi1; l2Energy = e2Energy1;
        l3Pt = e3Pt1; l3Eta = e3Eta1; l3Phi = e3Phi1; l3Energy = e3Energy1;
        l4Pt = e4Pt1; l4Eta = e4Eta1; l4Phi = e4Phi1; l4Energy = e4Energy1;
        l1PdgId = e1PdgId1; l2PdgId = e2PdgId1;
        l3PdgId = e3PdgId1; l4PdgId = e4PdgId1;
        l1ZZIso = e1ZZIso1; l2ZZIso = e2ZZIso1;
        l3ZZIso = e3ZZIso1; l4ZZIso = e4ZZIso1;
        l1ZZIsoPass = e1ZZIsoPass1; l2ZZIsoPass = e2ZZIsoPass1;
        l3ZZIsoPass = e3ZZIsoPass1; l4ZZIsoPass = e4ZZIsoPass1;
        l1SIP3D = e1SIP3D1; l2SIP3D = e2SIP3D1;
        l3SIP3D = e3SIP3D1; l4SIP3D = e4SIP3D1;
      }
      else if (_channel == "eemm"){
        if (std::fabs(e1_e2_Mass1 - Z_MASS) < std::fabs(m1_m2_Mass1 - Z_MASS)){
          Z1mass = e1_e2_Mass1; Z2mass = m1_m2_Mass1;
          l1Pt = e1Pt1; l1Eta = e1Eta1; l1Phi = e1Phi1; l1Energy = e1Energy1;
          l2Pt = e2Pt1; l2Eta = e2Eta1; l2Phi = e2Phi1; l2Energy = e2Energy1;
          l3Pt = m1Pt1; l3Eta = m1Eta1; l3Phi = m1Phi1; l3Energy = m1Energy1;
          l4Pt = m2Pt1; l4Eta = m2Eta1; l4Phi = m2Phi1; l4Energy = m2Energy1;
          l1PdgId = e1PdgId1; l2PdgId = e2PdgId1;
          l3PdgId = m1PdgId1; l4PdgId = m2PdgId1;
          l1ZZIso = e1ZZIso1; l2ZZIso = e2ZZIso1;
          l3ZZIso = m3ZZIso1; l4ZZIso = m4ZZIso1;
          l1ZZIsoPass = e1ZZIsoPass1; l2ZZIsoPass = e2ZZIsoPass1;
          l3ZZIsoPass = m3ZZIsoPass1; l4ZZIsoPass = m4ZZIsoPass1;
          l1SIP3D = e1SIP3D1; l2SIP3D = e2SIP3D1;
          l3SIP3D = m3SIP3D1; l4SIP3D = m4SIP3D1;
        }
        else{
          Z1mass = m1_m2_Mass1; Z2mass = e1_e2_Mass1;
          l1Pt = m1Pt1; l1Eta = m1Eta1; l1Phi = m1Phi1; l1Energy = m1Energy1;
          l2Pt = m2Pt1; l2Eta = m2Eta1; l2Phi = m2Phi1; l2Energy = m2Energy1;
          l3Pt = e1Pt1; l3Eta = e1Eta1; l3Phi = e1Phi1; l3Energy = e1Energy1;
          l4Pt = e2Pt1; l4Eta = e2Eta1; l4Phi = e2Phi1; l4Energy = e2Energy1;
          l1PdgId = m1PdgId1; l2PdgId = m2PdgId1;
          l3PdgId = e1PdgId1; l4PdgId = e2PdgId1;
          l1ZZIso = m1ZZIso1; l2ZZIso = m2ZZIso1;
          l3ZZIso = e3ZZIso1; l4ZZIso = e4ZZIso1;
          l1ZZIsoPass = m1ZZIsoPass1; l2ZZIsoPass = m2ZZIsoPass1;
          l3ZZIsoPass = e3ZZIsoPass1; l4ZZIsoPass = e4ZZIsoPass1;
          l1SIP3D = m1SIP3D1; l2SIP3D = m2SIP3D1;
          l3SIP3D = e3SIP3D1; l4SIP3D = e4SIP3D1;
        }
      }
      else if (_channel == "mmmm"){
        Z1mass = m1_m2_Mass1; Z2mass = m3_m4_Mass1;
        l1Pt = m1Pt1; l1Eta = m1Eta1; l1Phi = m1Phi1; l1Energy = m1Energy1;
        l2Pt = m2Pt1; l2Eta = m2Eta1; l2Phi = m2Phi1; l2Energy = m2Energy1;
        l3Pt = m3Pt1; l3Eta = m3Eta1; l3Phi = m3Phi1; l3Energy = m3Energy1;
        l4Pt = m4Pt1; l4Eta = m4Eta1; l4Phi = m4Phi1; l4Energy = m4Energy1;
        l1PdgId = m1PdgId1; l2PdgId = m2PdgId1;
        l3PdgId = m3PdgId1; l4PdgId = m4PdgId1;
        l1ZZIso = m1ZZIso1; l2ZZIso = m2ZZIso1;
        l3ZZIso = m3ZZIso1; l4ZZIso = m4ZZIso1;
        l1ZZIsoPass = m1ZZIsoPass1; l2ZZIsoPass = m2ZZIsoPass1;
        l3ZZIsoPass = m3ZZIsoPass1; l4ZZIsoPass = m4ZZIsoPass1;
        l1SIP3D = m1SIP3D1; l2SIP3D = m2SIP3D1;
        l3SIP3D = m3SIP3D1; l4SIP3D = m4SIP3D1;
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

      bool ZZIsoPass = l1ZZIsoPass && l2ZZIsoPass && l3ZZIsoPass && l4ZZIsoPass;
      if (Z1mass > 60.0 && Z1mass < 120.0 && Z2mass > 60.0 && Z2mass < 120.0 && ZZIsoPass){ 
        Float_t weight = 1.;
        if (_isT1MC){
          weight *= genWeight1;
          if (applyScaleFacs) weight *= GetScaleFactor(Ntuple::First);
          if (_doEE) weight *= 0.2;
          if (!_norm) weight /= summedWeights1;
        }

        InvMass4l_1->Fill(Mass1, weight);
        InvMass12_1->Fill(Z1mass, weight);
        InvMass34_1->Fill(Z2mass, weight);

        LepEnergy_1->Fill(l1Energy, weight);
        LepEnergy_1->Fill(l2Energy, weight);
        LepEnergy_1->Fill(l3Energy, weight);
        LepEnergy_1->Fill(l4Energy, weight);

        LepPt_1->Fill(l1Pt, weight);
        LepPt_1->Fill(l2Pt, weight);
        LepPt_1->Fill(l3Pt, weight);
        LepPt_1->Fill(l4Pt, weight);

        LepEta_1->Fill(l1Eta, weight);
        LepEta_1->Fill(l2Eta, weight);
        LepEta_1->Fill(l3Eta, weight);
        LepEta_1->Fill(l4Eta, weight);

        PolCosTheta12_1->Fill(GetPolCosTheta(lp1, ln1), weight);
        PolCosTheta34_1->Fill(GetPolCosTheta(lp2, ln2), weight);

        LepZZIsoVal_1->Fill(l1ZZIso, weight);
        LepZZIsoVal_1->Fill(l2ZZIso, weight);
        LepZZIsoVal_1->Fill(l3ZZIso, weight);
        LepZZIsoVal_1->Fill(l4ZZIso, weight);

        LepZZIsoPass_1->Fill(l1ZZIsoPass, weight);
        LepZZIsoPass_1->Fill(l2ZZIsoPass, weight);
        LepZZIsoPass_1->Fill(l3ZZIsoPass, weight);
        LepZZIsoPass_1->Fill(l4ZZIsoPass, weight);

        LepSIP3D_1->Fill(l1SIP3D, weight);
        LepSIP3D_1->Fill(l2SIP3D, weight);
        LepSIP3D_1->Fill(l3SIP3D, weight);
        LepSIP3D_1->Fill(l4SIP3D, weight);
      }
    }
    if (_ntuple2->GetEntry(i) && std::find(evts2.begin(), evts2.end(), evt2) == evts2.end()){
      if (!_isT2MC && _deduplicate) evts2.push_back(evt2);
      if (_channel == "eeee"){
        Z1mass = e1_e2_Mass2; Z2mass = e3_e4_Mass2;
        l1Pt = e1Pt2; l1Eta = e1Eta2; l1Phi = e1Phi2; l1Energy = e1Energy2;
        l2Pt = e2Pt2; l2Eta = e2Eta2; l2Phi = e2Phi2; l2Energy = e2Energy2;
        l3Pt = e3Pt2; l3Eta = e3Eta2; l3Phi = e3Phi2; l3Energy = e3Energy2;
        l4Pt = e4Pt2; l4Eta = e4Eta2; l4Phi = e4Phi2; l4Energy = e4Energy2;
        l1PdgId = e1PdgId2; l2PdgId = e2PdgId2;
        l3PdgId = e3PdgId2; l4PdgId = e4PdgId2;
        l1ZZIso = e1ZZIso2; l2ZZIso = e2ZZIso2;
        l3ZZIso = e3ZZIso2; l4ZZIso = e4ZZIso2;
        l1ZZIsoPass = e1ZZIsoPass2; l2ZZIsoPass = e2ZZIsoPass2;
        l3ZZIsoPass = e3ZZIsoPass2; l4ZZIsoPass = e4ZZIsoPass2;
        l1SIP3D = e1SIP3D2; l2SIP3D = e2SIP3D2;
        l3SIP3D = e3SIP3D2; l4SIP3D = e4SIP3D2;
      }
      else if (_channel == "eemm"){
        if (std::fabs(e1_e2_Mass2 - Z_MASS) < std::fabs(m1_m2_Mass2 - Z_MASS)){
          Z1mass = e1_e2_Mass2; Z2mass = m1_m2_Mass2;
          l1Pt = e1Pt2; l1Eta = e1Eta2; l1Phi = e1Phi2; l1Energy = e1Energy2;
          l2Pt = e2Pt2; l2Eta = e2Eta2; l2Phi = e2Phi2; l2Energy = e2Energy2;
          l3Pt = m1Pt2; l3Eta = m1Eta2; l3Phi = m1Phi2; l3Energy = m1Energy2;
          l4Pt = m2Pt2; l4Eta = m2Eta2; l4Phi = m2Phi2; l4Energy = m2Energy2;
          l1PdgId = e1PdgId2; l2PdgId = e2PdgId2;
          l3PdgId = m1PdgId2; l4PdgId = m2PdgId2;
          l1ZZIso = e1ZZIso2; l2ZZIso = e2ZZIso2;
          l3ZZIso = m3ZZIso2; l4ZZIso = m4ZZIso2;
          l1ZZIsoPass = e1ZZIsoPass2; l2ZZIsoPass = e2ZZIsoPass2;
          l3ZZIsoPass = m3ZZIsoPass2; l4ZZIsoPass = m4ZZIsoPass2;
          l1SIP3D = e1SIP3D2; l2SIP3D = e2SIP3D2;
          l3SIP3D = m3SIP3D2; l4SIP3D = m4SIP3D2;
        }
        else{
          Z1mass = m1_m2_Mass2; Z2mass = e1_e2_Mass2;
          l1Pt = m1Pt2; l1Eta = m1Eta2; l1Phi = m1Phi2; l1Energy = m1Energy2;
          l2Pt = m2Pt2; l2Eta = m2Eta2; l2Phi = m2Phi2; l2Energy = m2Energy2;
          l3Pt = e1Pt2; l3Eta = e1Eta2; l3Phi = e1Phi2; l3Energy = e1Energy2;
          l4Pt = e2Pt2; l4Eta = e2Eta2; l4Phi = e2Phi2; l4Energy = e2Energy2;
          l1PdgId = m1PdgId2; l2PdgId = m2PdgId2;
          l3PdgId = e1PdgId2; l4PdgId = e2PdgId2;
          l1ZZIso = m1ZZIso2; l2ZZIso = m2ZZIso2;
          l3ZZIso = e3ZZIso2; l4ZZIso = e4ZZIso2;
          l1ZZIsoPass = m1ZZIsoPass2; l2ZZIsoPass = m2ZZIsoPass2;
          l3ZZIsoPass = e3ZZIsoPass2; l4ZZIsoPass = e4ZZIsoPass2;
          l1SIP3D = m1SIP3D2; l2SIP3D = m2SIP3D2;
          l3SIP3D = e3SIP3D2; l4SIP3D = e4SIP3D2;
        }
      }
      else if (_channel == "mmmm"){
        Z1mass = m1_m2_Mass2; Z2mass = m3_m4_Mass2;
        l1Pt = m1Pt2; l1Eta = m1Eta2; l1Phi = m1Phi2; l1Energy = m1Energy2;
        l2Pt = m2Pt2; l2Eta = m2Eta2; l2Phi = m2Phi2; l2Energy = m2Energy2;
        l3Pt = m3Pt2; l3Eta = m3Eta2; l3Phi = m3Phi2; l3Energy = m3Energy2;
        l4Pt = m4Pt2; l4Eta = m4Eta2; l4Phi = m4Phi2; l4Energy = m4Energy2;
        l1PdgId = m1PdgId2; l2PdgId = m2PdgId2;
        l3PdgId = m3PdgId2; l4PdgId = m4PdgId2;
        l1ZZIso = m1ZZIso2; l2ZZIso = m2ZZIso2;
        l3ZZIso = m3ZZIso2; l4ZZIso = m4ZZIso2;
        l1ZZIsoPass = m1ZZIsoPass2; l2ZZIsoPass = m2ZZIsoPass2;
        l3ZZIsoPass = m3ZZIsoPass2; l4ZZIsoPass = m4ZZIsoPass2;
        l1SIP3D = m1SIP3D2; l2SIP3D = m2SIP3D2;
        l3SIP3D = m3SIP3D2; l4SIP3D = m4SIP3D2;
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

      bool ZZIsoPass = l1ZZIsoPass && l2ZZIsoPass && l3ZZIsoPass && l4ZZIsoPass;
      if (Z1mass > 60.0 && Z1mass < 120.0 && Z2mass > 60.0 && Z2mass < 120.0 && ZZIsoPass){ 
        Float_t weight = 1.;
        if (_isT2MC){
          weight *= genWeight2;
          if (applyScaleFacs) weight *= GetScaleFactor(Ntuple::Second);
          if (_doEE) weight *= 0.2;
          if (!_norm) weight /= summedWeights2;
        }

        InvMass4l_2->Fill(Mass2, weight);
        InvMass12_2->Fill(Z1mass, weight);
        InvMass34_2->Fill(Z2mass, weight);

        LepEnergy_2->Fill(l1Energy, weight);
        LepEnergy_2->Fill(l2Energy, weight);
        LepEnergy_2->Fill(l3Energy, weight);
        LepEnergy_2->Fill(l4Energy, weight);

        LepPt_2->Fill(l1Pt, weight);
        LepPt_2->Fill(l2Pt, weight);
        LepPt_2->Fill(l3Pt, weight);
        LepPt_2->Fill(l4Pt, weight);

        LepEta_2->Fill(l1Eta, weight);
        LepEta_2->Fill(l2Eta, weight);
        LepEta_2->Fill(l3Eta, weight);
        LepEta_2->Fill(l4Eta, weight);

        PolCosTheta12_2->Fill(GetPolCosTheta(lp1, ln1), weight);
        PolCosTheta34_2->Fill(GetPolCosTheta(lp2, ln2), weight);

        LepZZIsoVal_2->Fill(l1ZZIso, weight);
        LepZZIsoVal_2->Fill(l2ZZIso, weight);
        LepZZIsoVal_2->Fill(l3ZZIso, weight);
        LepZZIsoVal_2->Fill(l4ZZIso, weight);

        LepZZIsoPass_2->Fill(l1ZZIsoPass, weight);
        LepZZIsoPass_2->Fill(l2ZZIsoPass, weight);
        LepZZIsoPass_2->Fill(l3ZZIsoPass, weight);
        LepZZIsoPass_2->Fill(l4ZZIsoPass, weight);

        LepSIP3D_2->Fill(l1SIP3D, weight);
        LepSIP3D_2->Fill(l2SIP3D, weight);
        LepSIP3D_2->Fill(l3SIP3D, weight);
        LepSIP3D_2->Fill(l4SIP3D, weight);
      }
    }
    if (_doEE && _ntupleEE->GetEntry(i)){
      if (_channel == "eeee"){
        Z1mass = e1_e2_MassEE; Z2mass = e3_e4_MassEE;
        l1Pt = e1PtEE; l1Eta = e1EtaEE; l1Phi = e1PhiEE; l1Energy = e1EnergyEE;
        l2Pt = e2PtEE; l2Eta = e2EtaEE; l2Phi = e2PhiEE; l2Energy = e2EnergyEE;
        l3Pt = e3PtEE; l3Eta = e3EtaEE; l3Phi = e3PhiEE; l3Energy = e3EnergyEE;
        l4Pt = e4PtEE; l4Eta = e4EtaEE; l4Phi = e4PhiEE; l4Energy = e4EnergyEE;
        l1PdgId = e1PdgIdEE; l2PdgId = e2PdgIdEE;
        l3PdgId = e3PdgIdEE; l4PdgId = e4PdgIdEE;
        l1ZZIso = e1ZZIsoEE; l2ZZIso = e2ZZIsoEE;
        l3ZZIso = e3ZZIsoEE; l4ZZIso = e4ZZIsoEE;
        l1ZZIsoPass = e1ZZIsoPassEE; l2ZZIsoPass = e2ZZIsoPassEE;
        l3ZZIsoPass = e3ZZIsoPassEE; l4ZZIsoPass = e4ZZIsoPassEE;
        l1SIP3D = e1SIP3DEE; l2SIP3D = e2SIP3DEE;
        l3SIP3D = e3SIP3DEE; l4SIP3D = e4SIP3DEE;
      }
      else if (_channel == "eemm"){
        if (std::fabs(e1_e2_MassEE - Z_MASS) < std::fabs(m1_m2_MassEE - Z_MASS)){
          Z1mass = e1_e2_MassEE; Z2mass = m1_m2_MassEE;
          l1Pt = e1PtEE; l1Eta = e1EtaEE; l1Phi = e1PhiEE; l1Energy = e1EnergyEE;
          l2Pt = e2PtEE; l2Eta = e2EtaEE; l2Phi = e2PhiEE; l2Energy = e2EnergyEE;
          l3Pt = m1PtEE; l3Eta = m1EtaEE; l3Phi = m1PhiEE; l3Energy = m1EnergyEE;
          l4Pt = m2PtEE; l4Eta = m2EtaEE; l4Phi = m2PhiEE; l4Energy = m2EnergyEE;
          l1PdgId = e1PdgIdEE; l2PdgId = e2PdgIdEE;
          l3PdgId = m1PdgIdEE; l4PdgId = m2PdgIdEE;
          l1ZZIso = e1ZZIsoEE; l2ZZIso = e2ZZIsoEE;
          l3ZZIso = m3ZZIsoEE; l4ZZIso = m4ZZIsoEE;
          l1ZZIsoPass = e1ZZIsoPassEE; l2ZZIsoPass = e2ZZIsoPassEE;
          l3ZZIsoPass = m3ZZIsoPassEE; l4ZZIsoPass = m4ZZIsoPassEE;
          l1SIP3D = e1SIP3DEE; l2SIP3D = e2SIP3DEE;
          l3SIP3D = m3SIP3DEE; l4SIP3D = m4SIP3DEE;
        }
        else{
          Z1mass = m1_m2_MassEE; Z2mass = e1_e2_MassEE;
          l1Pt = m1PtEE; l1Eta = m1EtaEE; l1Phi = m1PhiEE; l1Energy = m1EnergyEE;
          l2Pt = m2PtEE; l2Eta = m2EtaEE; l2Phi = m2PhiEE; l2Energy = m2EnergyEE;
          l3Pt = e1PtEE; l3Eta = e1EtaEE; l3Phi = e1PhiEE; l3Energy = e1EnergyEE;
          l4Pt = e2PtEE; l4Eta = e2EtaEE; l4Phi = e2PhiEE; l4Energy = e2EnergyEE;
          l1PdgId = m1PdgIdEE; l2PdgId = m2PdgIdEE;
          l3PdgId = e1PdgIdEE; l4PdgId = e2PdgIdEE;
          l1ZZIso = m1ZZIsoEE; l2ZZIso = m2ZZIsoEE;
          l3ZZIso = e3ZZIsoEE; l4ZZIso = e4ZZIsoEE;
          l1ZZIsoPass = m1ZZIsoPassEE; l2ZZIsoPass = m2ZZIsoPassEE;
          l3ZZIsoPass = e3ZZIsoPassEE; l4ZZIsoPass = e4ZZIsoPassEE;
          l1SIP3D = m1SIP3DEE; l2SIP3D = m2SIP3DEE;
          l3SIP3D = e3SIP3DEE; l4SIP3D = e4SIP3DEE;
        }
      }
      else if (_channel == "mmmm"){
        Z1mass = m1_m2_MassEE; Z2mass = m3_m4_MassEE;
        l1Pt = m1PtEE; l1Eta = m1EtaEE; l1Phi = m1PhiEE; l1Energy = m1EnergyEE;
        l2Pt = m2PtEE; l2Eta = m2EtaEE; l2Phi = m2PhiEE; l2Energy = m2EnergyEE;
        l3Pt = m3PtEE; l3Eta = m3EtaEE; l3Phi = m3PhiEE; l3Energy = m3EnergyEE;
        l4Pt = m4PtEE; l4Eta = m4EtaEE; l4Phi = m4PhiEE; l4Energy = m4EnergyEE;
        l1PdgId = m1PdgIdEE; l2PdgId = m2PdgIdEE;
        l3PdgId = m3PdgIdEE; l4PdgId = m4PdgIdEE;
        l1ZZIso = m1ZZIsoEE; l2ZZIso = m2ZZIsoEE;
        l3ZZIso = m3ZZIsoEE; l4ZZIso = m4ZZIsoEE;
        l1ZZIsoPass = m1ZZIsoPassEE; l2ZZIsoPass = m2ZZIsoPassEE;
        l3ZZIsoPass = m3ZZIsoPassEE; l4ZZIsoPass = m4ZZIsoPassEE;
        l1SIP3D = m1SIP3DEE; l2SIP3D = m2SIP3DEE;
        l3SIP3D = m3SIP3DEE; l4SIP3D = m4SIP3DEE;
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

      bool ZZIsoPass = l1ZZIsoPass && l2ZZIsoPass && l3ZZIsoPass && l4ZZIsoPass;
      if (Z1mass > 60.0 && Z1mass < 120.0 && Z2mass > 60.0 && Z2mass < 120.0 && ZZIsoPass){ 
        Float_t weight = genWeightEE;
        if (applyScaleFacs) weight *= GetScaleFactor(Ntuple::EE);
        weight *= 0.8;
        if (!_norm) weight /= summedWeightsEE;

        if (_isT1MC){
          InvMass4l_1->Fill(MassEE, weight);
          InvMass12_1->Fill(Z1mass, weight);
          InvMass34_1->Fill(Z2mass, weight);

          LepEnergy_1->Fill(l1Energy, weight);
          LepEnergy_1->Fill(l2Energy, weight);
          LepEnergy_1->Fill(l3Energy, weight);
          LepEnergy_1->Fill(l4Energy, weight);

          LepPt_1->Fill(l1Pt, weight);
          LepPt_1->Fill(l2Pt, weight);
          LepPt_1->Fill(l3Pt, weight);
          LepPt_1->Fill(l4Pt, weight);

          LepEta_1->Fill(l1Eta, weight);
          LepEta_1->Fill(l2Eta, weight);
          LepEta_1->Fill(l3Eta, weight);
          LepEta_1->Fill(l4Eta, weight);

          PolCosTheta12_1->Fill(GetPolCosTheta(lp1, ln1), weight);
          PolCosTheta34_1->Fill(GetPolCosTheta(lp2, ln2), weight);

          LepZZIsoVal_1->Fill(l1ZZIso, weight);
          LepZZIsoVal_1->Fill(l2ZZIso, weight);
          LepZZIsoVal_1->Fill(l3ZZIso, weight);
          LepZZIsoVal_1->Fill(l4ZZIso, weight);

          LepZZIsoPass_1->Fill(l1ZZIsoPass, weight);
          LepZZIsoPass_1->Fill(l2ZZIsoPass, weight);
          LepZZIsoPass_1->Fill(l3ZZIsoPass, weight);
          LepZZIsoPass_1->Fill(l4ZZIsoPass, weight);

          LepSIP3D_1->Fill(l1SIP3D, weight);
          LepSIP3D_1->Fill(l2SIP3D, weight);
          LepSIP3D_1->Fill(l3SIP3D, weight);
          LepSIP3D_1->Fill(l4SIP3D, weight);
        }
        else{
          InvMass4l_2->Fill(MassEE, weight);
          InvMass12_2->Fill(Z1mass, weight);
          InvMass34_2->Fill(Z2mass, weight);

          LepEnergy_2->Fill(l1Energy, weight);
          LepEnergy_2->Fill(l2Energy, weight);
          LepEnergy_2->Fill(l3Energy, weight);
          LepEnergy_2->Fill(l4Energy, weight);

          LepPt_2->Fill(l1Pt, weight);
          LepPt_2->Fill(l2Pt, weight);
          LepPt_2->Fill(l3Pt, weight);
          LepPt_2->Fill(l4Pt, weight);

          LepEta_2->Fill(l1Eta, weight);
          LepEta_2->Fill(l2Eta, weight);
          LepEta_2->Fill(l3Eta, weight);
          LepEta_2->Fill(l4Eta, weight);

          PolCosTheta12_2->Fill(GetPolCosTheta(lp1, ln1), weight);
          PolCosTheta34_2->Fill(GetPolCosTheta(lp2, ln2), weight);

          LepZZIsoVal_2->Fill(l1ZZIso, weight);
          LepZZIsoVal_2->Fill(l2ZZIso, weight);
          LepZZIsoVal_2->Fill(l3ZZIso, weight);
          LepZZIsoVal_2->Fill(l4ZZIso, weight);

          LepZZIsoPass_2->Fill(l1ZZIsoPass, weight);
          LepZZIsoPass_2->Fill(l2ZZIsoPass, weight);
          LepZZIsoPass_2->Fill(l3ZZIsoPass, weight);
          LepZZIsoPass_2->Fill(l4ZZIsoPass, weight);

          LepSIP3D_2->Fill(l1SIP3D, weight);
          LepSIP3D_2->Fill(l2SIP3D, weight);
          LepSIP3D_2->Fill(l3SIP3D, weight);
          LepSIP3D_2->Fill(l4SIP3D, weight);
        }
      }
    }
  }
  std::cout << "End looping." << std::endl;

  // Writing README
  if (plot || _makePlots){
    std::ofstream readme((dirname + "/README").c_str());
    if (readme.is_open()){
      readme << _label1 << " Events: " << InvMass4l_1->GetEntries() << std::endl;
      readme << _label2 << " Events: " << InvMass4l_2->GetEntries() << std::endl;
      if (_isT1MC && !_norm){
        _ntuple1->GetEntry(0);
        float histScaling1 = _kfac1 * _xsec1 * _lumi / summedWeights1;
        readme << _label1 << " Scaling: " << histScaling1 << std::endl;
        readme << "  - Cross-section: " << _xsec1 << " fb" << std::endl;
        readme << "  - Luminosity: " << _lumi << " fb-1" << std::endl;
        readme << "  - Summed Weights: " << summedWeights1 << std::endl;
        readme << "  - k-factor: " << _kfac1 << std::endl;
      }
      if (_isT2MC && !_norm){
        _ntuple2->GetEntry(0);
        float histScaling2 = _kfac2 * _xsec2 * _lumi / summedWeights2;
        readme << _label2 << " Scaling: " << histScaling2 << std::endl;
        readme << "  - Cross-section: " << _xsec2 << " fb" << std::endl;
        readme << "  - Luminosity: " << _lumi << " fb-1" << std::endl;
        readme << "  - Summed Weights: " << summedWeights2 << std::endl;
        readme << "  - k-factor: " << _kfac2 << std::endl;
      }
      if (applyScaleFacs) readme << "Efficiency SFs applied." << std::endl;
      if (_doEE) readme << "MC pre- and postEE samples weighted." << std::endl;
      readme.close();
    }
  }

  // Creating subdirs in histfile
  histout->cd();
  histout->rmdir(_channel.c_str());
  TDirectory *subdir = histout->mkdir(_channel.c_str());
  TDirectory *subsubdir1 = nullptr, *subsubdir2 = nullptr;
  if (subdir != nullptr){
    subdir->cd();
    subsubdir1 = subdir->mkdir(_label1.c_str());
    subsubdir2 = subdir->mkdir(_label2.c_str());
  }

  // Creating MC scale factors
  // Division by summedWeights is done when filling histos bc of EE analysis
  float histScaling1 = _kfac1 * _xsec1 * _lumi;
  float histScaling2 = _kfac2 * _xsec2 * _lumi;

  // Defining draw functions
  std::vector<std::string> labels = {_label1, _label2};
  auto Draw = [&c, &dirname, _filetype=_filetype, &labels](TH1* hist, Ntuple ntuple){
    int idx = ntuple == Ntuple::First? 0 : 1;
    if (idx < 0) idx = 0;
    else if (idx >= labels.size()) idx = labels.size()-1;
    hist->Draw("hist");
    c->SaveAs((dirname + "/" + labels[idx] + "/" + hist->GetName() + _filetype).c_str());
  };
	auto DrawSame = [&c, &dirname, _filetype=_filetype, &labels](TH1* hist1, TH1* hist2){
		std::string histname = hist1->GetName();
		histname = histname.substr(0, histname.size()-2);
		TH1 *hist1_clone = (TH1*)hist1->Clone(histname.c_str());
		std::string title = hist1->GetTitle();
		if (title.find(labels[0]) != std::string::npos){
			title = title.substr(0, title.find(labels[0])-2);
			hist1_clone->SetTitle(title.c_str());
		}
		Float_t max = (hist1->GetMaximum() > hist2->GetMaximum())?
			hist1->GetMaximum() : hist2->GetMaximum();
		hist1_clone->SetMaximum(max*1.2);
		hist2->SetMaximum(max*1.2);

		hist1_clone->Draw("hist");
		hist2->SetLineColor(kRed);
		hist2->Draw("hist same");

		UInt_t idx = 0;
		TLegend *legend = c->BuildLegend();
		for (const TObject *obj : *legend->GetListOfPrimitives())
			((TLegendEntry*)obj)->SetLabel(labels[idx++].c_str());
		c->SaveAs((dirname + "/" + histname + _filetype).c_str());
	};

  // Finalizing
  for (auto histname : histnames){
    TH1F *hist1 = (TH1F*)gROOT->Get((histname + "_1").c_str());
    TH1F *hist2 = (TH1F*)gROOT->Get((histname + "_2").c_str());
    if (hist1 == nullptr || hist2 == nullptr){
      std::cout << "Couldn't find hist(s). Skipping " << histname << std::endl;
      break;
    }
    // Formatting
    hist1->SetMinimum(0);
    hist2->SetMinimum(0);

    // Scaling
    if (!_norm){
      if (_isT1MC) hist1->Scale(histScaling1);
      if (_isT2MC) hist2->Scale(histScaling2);
    }
    else{
      hist1->Scale(1.0/hist1->Integral());
      hist2->Scale(1.0/hist2->Integral());
    }

    // Drawing
    if (plot || _makePlots){
      c->cd();
      Draw(hist1, Ntuple::First);
      Draw(hist2, Ntuple::Second);
      DrawSame(hist1, hist2);
    }

    // Writing
    histout->cd();
    if (subsubdir1 != nullptr) subsubdir1->cd();
    hist1->Write();
    if (subsubdir2 != nullptr) subsubdir2->cd();
    hist2->Write();

    delete hist1;
    delete hist2;
  }
  
  // Ratio Plots
  if (ratioplot || _makeRatios){
  }

  c->Close();
  histout->Close();
}

#ifndef __CLING__
#include "interface/argparse.h"
int main(int nargs, char *argv[]){
  float kfac;
  kfac = 1.2; //qqZZ kfac ~1.1 plus esimated missing ggZZ signal
  //kfac = 1.0835 * 1.01587 * 1.0004874;//qqZZ kfac * expected ggZZ signal * expected EWK signal
  kfac = 1.0835;
  std::stringstream ss_kfac; ss_kfac << kfac;

  auto parser = argparse::ArgumentParser(nargs, argv)
    .formatter_class(argparse::HelpFormatter::ArgumentDefaults);

  parser.add_argument<bool>("--debug").def("false")
    .help("if true, perform extra code for testing purposes");
  parser.add_argument<bool>("--mc1").def("false")
    .help("if true, consider first input to be MC");
  parser.add_argument<bool>("--mc2").def("false")
    .help("if true, consider second input to be MC");
  parser.add_argument<double>("-l", "--lumi").def("37.846170084") //2022
    .help("data luminosity to scale MC in fb-1");
  parser.add_argument<double>("-x", "--xsec1").def("1256") //qqZZ
    .help("cross-section of given MC process in fb for first input");
  parser.add_argument<double>("-X", "--xsec2").def("1256") //qqZZ
    .help("cross-section of given MC process in fb for second input");
  parser.add_argument<double>("-k", "--kfac1").def(ss_kfac.str()) //qqZZ kfac ~1.1 plus missing ggZZ signal
    .help("k-factor of given MC process for first input");
  parser.add_argument<double>("-K", "--kfac2").def(ss_kfac.str()) //qqZZ kfac ~1.1 plus missing ggZZ signal
    .help("k-factor of given MC process for second input");
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
  parser.add_argument<bool>("--EE").def("false")
    .help("if true, the looper will analyze MC samples with weights for pre- and postEE comparison.");
  parser.add_argument<bool>("-r", "--recreate").def("false")
    .help("if true, the output file will be recreated for the given channel(s)");
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
  int index=0;
  for (std::string channel : channels){
    CompLooper l(args["name"].c_str(), channel.c_str(), args["label1"].c_str(), args["label2"].c_str());
    if ((index++)==0 && args["recreate"]) l.SetMode("recreate");
    l.SetMakePlots(!args["noplots"]);
    l.SetDebug(args["debug"]);
    l.SetNorm(args["norm"]);
    l.SetDeduplicate(args["deduplicate"]);
    l.SetPlotFiletype(args["filetype"]);
    l.SetEE(args["EE"]);
    
    if (args["mc1"].is_true()){
      l.SetMC1();
      l.SetLumi(args["lumi"]);
      l.SetXsec1(args["xsec1"]);
      l.SetKfac1(args["kfac1"]);
    }
    if (args["mc2"].is_true()){
      l.SetMC2();
      l.SetLumi(args["lumi"]);
      l.SetXsec2(args["xsec2"]);
      l.SetKfac2(args["kfac2"]);
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
