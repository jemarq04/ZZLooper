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
    enum Ntuple {First, Second};
  public:
    CompLooper(const char *name, const char *channel, const char *label1, const char *label2);
    ~CompLooper();

    //Debug
    void SetNorm(bool norm=true){_norm = norm;}
    void SetDeduplicate(bool deduplicate=true){_deduplicate=deduplicate;}

    void SetMode(std::string mode){_mode = mode;}
    void SetMakePlots(bool val=true);
    void SetMakeRatios(bool val=true);
    void SetPlotFiletype(std::string ft=".png");

    double GetScaleFactor(Ntuple ntuple);

    void Loop(bool applyScaleFacs=false);
  private:
    std::string FindFile(const char *label);
    static constexpr const char *CHANNEL = "eemm";
    bool _makePlots = false, _makeRatios = false;
    bool _norm = false, _deduplicate = false;
    std::string _label1, _label2;
    std::string _filetype = ".png", _mode = "UPDATE";
#ifndef NOSF
    std::unique_ptr<correction::CorrectionSet> _pileupSF, _eIdSF, _eRecoSF, _mIdSF;
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

void CompLooper::SetMakePlots(bool val){_makePlots = val;}
void CompLooper::SetMakeRatios(bool val){_makeRatios = val;}
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
      if (ntuple == Ntuple::First){
        if (e1Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt1), e1Eta1, e1Pt1});
        if (e2Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt1), e2Eta1, e2Pt1});
        if (e3Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e3Pt1), e3Eta1, e3Pt1});
        if (e4Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e4Pt1), e4Eta1, e4Pt1});
      }
      else{
        if (e1Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt2), e1Eta2, e1Pt2});
        if (e2Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt2), e2Eta2, e2Pt2});
        if (e3Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e3Pt2), e3Eta2, e3Pt2});
        if (e4Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e4Pt2), e4Eta2, e4Pt2});
      }
    }
  }
  else if (_channel == "eemm"){
    if (_eIdSF != nullptr){}
    if (_eRecoSF != nullptr){
      const auto recoref = (*_eRecoSF->begin()).second;
      if (ntuple == Ntuple::First){
        if (e1Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt1), e1Eta1, e1Pt1});
        if (e2Pt1 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt1), e2Eta1, e2Pt1});
      }
      else{
        if (e1Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e1Pt2), e1Eta2, e1Pt2});
        if (e2Pt2 > 10) weight *= recoref->evaluate({"2022Re-recoBCD", "sf", GetEleRecoSFName(e2Pt2), e2Eta2, e2Pt2});
      }
    }
    if (_mIdSF != nullptr){
      if (ntuple == Ntuple::First){
        if (m1Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta1), m1Pt1, "nominal"});
        if (m2Pt1 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta1), m2Pt1, "nominal"});
      }
      else{
        if (m1Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta2), m1Pt2, "nominal"});
        if (m2Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta2), m2Pt2, "nominal"});
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
      else{
        if (m1Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m1Eta2), m1Pt2, "nominal"});
        if (m2Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m2Eta2), m2Pt2, "nominal"});
        if (m3Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m3Eta2), m3Pt2, "nominal"});
        if (m4Pt2 > 15) weight *= _mIdSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({std::fabs(m4Eta2), m4Pt2, "nominal"});
      }
    }
  }
  if (_pileupSF != nullptr){
    if (ntuple == Ntuple::First)
      weight *= (*_pileupSF->begin()).second->evaluate({nTruePU1, "nominal"});
    else
      weight *= (*_pileupSF->begin()).second->evaluate({nTruePU2, "nominal"});
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
    //_eIdSF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22/electron.json.gz");
    _eRecoSF = correction::CorrectionSet::from_file(basename + "/EGM/2022_Summer22/electron.json.gz");
    _mIdSF = correction::CorrectionSet::from_file(basename + "/MUO/2022_Summer22/muon_Z.json.gz");
    //if (_pileupSF != nullptr) std::cout << "Applying pileupSF..." << std::endl;
    //if (_eIdSF != nullptr) std::cout << "Applying eIdSF..." << std::endl;
    //if (_eRecoSF != nullptr) std::cout << "Applying eRecoSF..." << std::endl;
    //if (_mIdSF != nullptr) std::cout << "Applying mIdSF..." << std::endl;
  }
#endif

  // Binning
  std::vector<Double_t> InvMass_4l_binning   = {100.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0};
  std::vector<Double_t> InvMass_pair_binning = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120.};
  std::vector<Double_t> LepPt_binning        = {0., 30., 50., 100., 200.};

  // Declaring Histograms
  TH1F *InvMass4l_1 = new TH1F("InvMass4l_1", TString::Format("4-Lepton Invariant Mass (%s)", _label1.c_str()), InvMass_4l_binning.size()-1, &InvMass_4l_binning[0]);
  TH1F *InvMass12_1 = new TH1F("InvMass12_1", TString::Format("Primary Lepton Pair Invariant Mass (%s)", _label1.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34_1 = new TH1F("InvMass34_1", TString::Format("Secondary Lepton Pair Invariant Mass (%s)", _label1.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy_1 = new TH1F("LepEnergy_1", TString::Format("Lepton Energy (%s)", _label1.c_str()), 60, 0, 600);
  TH1F *LepPt_1 = new TH1F("LepPt_1", TString::Format("Lepton Transverse Momentum (%s)", _label1.c_str()), LepPt_binning.size()-1, &LepPt_binning[0]);
  TH1F *PolCosTheta12_1 = new TH1F("PolCosTheta12_1", TString::Format("cos#theta_{12}* (%s)", _label1.c_str()), 50, -1, 1);
  TH1F *PolCosTheta34_1 = new TH1F("PolCosTheta34_1", TString::Format("cos#theta_{34}* (%s)", _label1.c_str()), 50, -1, 1);
  SetTitles(InvMass4l_1, "m_{4l} [GeV]");
  SetTitles(InvMass12_1, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34_1, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy_1, "E_{l} [GeV]");
  SetTitles(LepPt_1, "p_{l,T} [GeV]");
  SetTitles(PolCosTheta12_1, "cos#theta_{12}*");
  SetTitles(PolCosTheta34_1, "cos#theta_{34}*");

  TH1F *InvMass4l_2 = new TH1F("InvMass4l_2", TString::Format("4-Lepton Invariant Mass (%s)", _label2.c_str()), InvMass_4l_binning.size()-1, &InvMass_4l_binning[0]);
  TH1F *InvMass12_2 = new TH1F("InvMass12_2", TString::Format("Primary Lepton Pair Invariant Mass (%s)", _label2.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *InvMass34_2 = new TH1F("InvMass34_2", TString::Format("Secondary Lepton Pair Invariant Mass (%s)", _label2.c_str()), InvMass_pair_binning.size()-1, &InvMass_pair_binning[0]);
  TH1F *LepEnergy_2 = new TH1F("LepEnergy_2", TString::Format("Lepton Energy (%s)", _label2.c_str()), 60, 0, 600);
  TH1F *LepPt_2 = new TH1F("LepPt_2", TString::Format("Lepton Transverse Momentum (%s)", _label2.c_str()), LepPt_binning.size()-1, &LepPt_binning[0]);
  TH1F *PolCosTheta12_2 = new TH1F("PolCosTheta12_2", TString::Format("cos#theta_{12}* (%s)", _label2.c_str()), 50, -1, 1);
  TH1F *PolCosTheta34_2 = new TH1F("PolCosTheta34_2", TString::Format("cos#theta_{34}* (%s)", _label2.c_str()), 50, -1, 1);
  SetTitles(InvMass4l_2, "m_{4l} [GeV]");
  SetTitles(InvMass12_2, "m_{Z_{1}} [GeV]");
  SetTitles(InvMass34_2, "m_{Z_{2}} [GeV]");
  SetTitles(LepEnergy_2, "E_{l} [GeV]");
  SetTitles(LepPt_2, "p_{l,T} [GeV]");
  SetTitles(PolCosTheta12_2, "cos#theta_{12}*");
  SetTitles(PolCosTheta34_2, "cos#theta_{34}*");

  Long64_t nentries1 = _ntuple1->GetEntries();
  Long64_t nentries2 = _ntuple2->GetEntries();
  Long64_t nentries = nentries1 > nentries2? nentries1 : nentries2;

  ROOT::Math::PtEtaPhiEVector lp1, ln1, lp2, ln2;
  Float_t Z1mass, Z2mass;
  Float_t l1Pt, l1Eta, l1Phi, l1Energy;
  Float_t l2Pt, l2Eta, l2Phi, l2Energy;
  Float_t l3Pt, l3Eta, l3Phi, l3Energy;
  Float_t l4Pt, l4Eta, l4Phi, l4Energy;
  Int_t l1PdgId, l2PdgId, l3PdgId, l4PdgId;
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
        }
        else{
          Z1mass = m1_m2_Mass1; Z2mass = e1_e2_Mass1;
          l1Pt = m1Pt1; l1Eta = m1Eta1; l1Phi = m1Phi1; l1Energy = m1Energy1;
          l2Pt = m2Pt1; l2Eta = m2Eta1; l2Phi = m2Phi1; l2Energy = m2Energy1;
          l3Pt = e1Pt1; l3Eta = e1Eta1; l3Phi = e1Phi1; l3Energy = e1Energy1;
          l4Pt = e2Pt1; l4Eta = e2Eta1; l4Phi = e2Phi1; l4Energy = e2Energy1;
          l1PdgId = m1PdgId1; l2PdgId = m2PdgId1;
          l3PdgId = e1PdgId1; l4PdgId = e2PdgId1;
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
        Float_t weight = _isT1MC? genWeight1 : 1.;
        if (applyScaleFacs && _isT1MC) weight *= GetScaleFactor(Ntuple::First);

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

        PolCosTheta12_1->Fill(GetPolCosTheta(lp1, ln1), weight);
        PolCosTheta34_1->Fill(GetPolCosTheta(lp2, ln2), weight);
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
        }
        else{
          Z1mass = m1_m2_Mass2; Z2mass = e1_e2_Mass2;
          l1Pt = m1Pt2; l1Eta = m1Eta2; l1Phi = m1Phi2; l1Energy = m1Energy2;
          l2Pt = m2Pt2; l2Eta = m2Eta2; l2Phi = m2Phi2; l2Energy = m2Energy2;
          l3Pt = e1Pt2; l3Eta = e1Eta2; l3Phi = e1Phi2; l3Energy = e1Energy2;
          l4Pt = e2Pt2; l4Eta = e2Eta2; l4Phi = e2Phi2; l4Energy = e2Energy2;
          l1PdgId = m1PdgId2; l2PdgId = m2PdgId2;
          l3PdgId = e1PdgId2; l4PdgId = e2PdgId2;
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
        Float_t weight = _isT2MC? genWeight2 : 1.;
        if (applyScaleFacs && _isT2MC) weight *= GetScaleFactor(Ntuple::Second);

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

        PolCosTheta12_2->Fill(GetPolCosTheta(lp1, ln1), weight);
        PolCosTheta34_2->Fill(GetPolCosTheta(lp2, ln2), weight);
      }
    }
  }
  std::cout << "End looping." << std::endl;

  // Formatting
  InvMass4l_1->SetMinimum(0);
  InvMass12_1->SetMinimum(0);
  InvMass34_1->SetMinimum(0);
  LepEnergy_1->SetMinimum(0);
  LepPt_1->SetMinimum(0);
  PolCosTheta12_1->SetMinimum(0);
  PolCosTheta34_1->SetMinimum(0);

  InvMass4l_2->SetMinimum(0);
  InvMass12_2->SetMinimum(0);
  InvMass34_2->SetMinimum(0);
  LepEnergy_2->SetMinimum(0);
  LepPt_2->SetMinimum(0);
  PolCosTheta12_2->SetMinimum(0);
  PolCosTheta34_2->SetMinimum(0);

  // Scaling
  if (!_norm){
    if (_isT1MC){
      _ntuple1->GetEntry(0);
      //std::cout << "summedWeights1=" << summedWeights1 << std::endl;
      float histScaling1 = _kfac1 * _xsec1 * _lumi / summedWeights1;
      InvMass4l_1->Scale(histScaling1);
      InvMass12_1->Scale(histScaling1);
      InvMass34_1->Scale(histScaling1);
      LepEnergy_1->Scale(histScaling1);
      LepPt_1->Scale(histScaling1);
      PolCosTheta12_1->Scale(histScaling1);
      PolCosTheta34_1->Scale(histScaling1);
    }
    if (_isT2MC){
      _ntuple2->GetEntry(0);
      //std::cout << "summedWeights2=" << summedWeights2 << std::endl;
      float histScaling2 = _kfac2 * _xsec2 * _lumi / summedWeights2;
      InvMass4l_2->Scale(histScaling2);
      InvMass12_2->Scale(histScaling2);
      InvMass34_2->Scale(histScaling2);
      LepEnergy_2->Scale(histScaling2);
      LepPt_2->Scale(histScaling2);
      PolCosTheta12_2->Scale(histScaling2);
      PolCosTheta34_2->Scale(histScaling2);
    }
  }
  else{
    InvMass4l_1->Scale(1.0/InvMass4l_1->Integral());
    InvMass12_1->Scale(1.0/InvMass12_1->Integral());
    InvMass34_1->Scale(1.0/InvMass34_1->Integral());
    LepEnergy_1->Scale(1.0/LepEnergy_1->Integral()); 
    LepPt_1->Scale(1.0/LepPt_1->Integral());
    PolCosTheta12_1->Scale(1.0/PolCosTheta12_1->Integral());
    PolCosTheta34_1->Scale(1.0/PolCosTheta34_1->Integral());

    InvMass4l_2->Scale(1.0/InvMass4l_2->Integral());
    InvMass12_2->Scale(1.0/InvMass12_2->Integral());
    InvMass34_2->Scale(1.0/InvMass34_2->Integral());
    LepEnergy_2->Scale(1.0/LepEnergy_2->Integral()); 
    LepPt_2->Scale(1.0/LepPt_2->Integral());
    PolCosTheta12_2->Scale(1.0/PolCosTheta12_2->Integral());
    PolCosTheta34_2->Scale(1.0/PolCosTheta34_2->Integral());
  }
  
  // Drawing
  if (plot || _makePlots){
    c->cd();
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

    Draw(InvMass4l_1, Ntuple::First);
    Draw(InvMass12_1, Ntuple::First);
    Draw(InvMass34_1, Ntuple::First);
    Draw(LepEnergy_1, Ntuple::First);
    Draw(LepPt_1, Ntuple::First);
    Draw(PolCosTheta12_1, Ntuple::First);
    Draw(PolCosTheta34_1, Ntuple::First);

    Draw(InvMass4l_2, Ntuple::Second);
    Draw(InvMass12_2, Ntuple::Second);
    Draw(InvMass34_2, Ntuple::Second);
    Draw(LepEnergy_2, Ntuple::Second);
    Draw(LepPt_2, Ntuple::Second);
    Draw(PolCosTheta12_2, Ntuple::Second);
    Draw(PolCosTheta34_2, Ntuple::Second);

    TLegend *legend;

    DrawSame(InvMass4l_1, InvMass4l_2);
    DrawSame(InvMass12_1, InvMass12_2);
    DrawSame(InvMass34_1, InvMass34_2);
    DrawSame(LepEnergy_1, LepEnergy_2);
    DrawSame(LepPt_1, LepPt_2);
    DrawSame(PolCosTheta12_1, PolCosTheta12_2);
    DrawSame(PolCosTheta34_1, PolCosTheta34_2);
  }

  // Writing
  histout->cd();
  histout->rmdir(_channel.c_str());
  TDirectory *subdir = histout->mkdir(_channel.c_str());
  if (subdir != nullptr){
    subdir->cd();
    subdir->rmdir(_label1.c_str());
    TDirectory *subsubdir = subdir->mkdir(_label1.c_str());
    if (subsubdir != nullptr) subsubdir->cd();
  }

  InvMass4l_1->Write();
  InvMass12_1->Write();
  InvMass34_1->Write();
  LepEnergy_1->Write();
  LepPt_1->Write();
  PolCosTheta12_1->Write();
  PolCosTheta34_1->Write();

  if (subdir != nullptr){
    subdir->cd();
    subdir->rmdir(_label2.c_str());
    TDirectory *subsubdir = subdir->mkdir(_label2.c_str());
    if (subsubdir != nullptr) subsubdir->cd();
  }

  InvMass4l_2->Write();
  InvMass12_2->Write();
  InvMass34_2->Write();
  LepEnergy_2->Write();
  LepPt_2->Write();
  PolCosTheta12_2->Write();
  PolCosTheta34_2->Write();

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
  parser.add_argument<double>("-x", "--xsec1").def("1390") //qqZZ
    .help("cross-section of given MC process in fb for first input");
  parser.add_argument<double>("-X", "--xsec2").def("1390") //qqZZ
    .help("cross-section of given MC process in fb for second input");
  parser.add_argument<double>("-k", "--kfac1").def("1.2") //qqZZ kfac ~1.1 plus missing ggZZ signal
    .help("k-factor of given MC process for first input");
  parser.add_argument<double>("-K", "--kfac2").def("1.2") //qqZZ kfac ~1.1 plus missing ggZZ signal
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
    if ((index++)==0 && args["recreate"].is_true()) l.SetMode("recreate");
    l.SetMakePlots(!args["noplots"].is_true());
    l.SetNorm(args["norm"]);
    l.SetDeduplicate(args["deduplicate"]);
    l.SetPlotFiletype(args["filetype"]);
    
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
