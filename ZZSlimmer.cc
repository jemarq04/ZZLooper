#include "interface/ZZSlimmerBase.h"
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

//Slims and applies selections to two ntuples and saves to output for CompLooper
class ZZSlimmer : public ZZSlimmerBase {
  public:
    ZZSlimmer(const char *name, const char *channel, const char *filenames="<NONE>");
    ~ZZSlimmer();

    void SetLabels(std::string title1, std::string title2){
      _label1 = title1;
      _label2 = title2;
    }
    void SetMode(std::string mode){_mode = mode;}
    void SetNoCuts(bool nocuts=true){_nocuts = nocuts;}
    void Slim();
  private:
    std::string _label1, _label2;
    std::string _mode = "update";
    bool _nocuts = false;
};

ZZSlimmer::ZZSlimmer(const char *name, const char *channel, const char *filenames)
    : ZZSlimmerBase(name, channel, filenames){
};

ZZSlimmer::~ZZSlimmer(){}

void ZZSlimmer::Slim(){
  // Setup
  TFile *outfile = new TFile(("slimmed/" + _name + ".root").c_str(), _mode.c_str());
  if (!outfile->IsOpen()) return;
  outfile->cd();

  if (_verbose) std::cout << "Slimming " << _channel << " channel..." << std::endl;
  
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
  
  // Activate certain branches
  std::vector<std::string> branchnames = {"evt", "run", "nvtx", "lumi", "Mass"};
  for (std::string var : {"Pt", "Eta", "Phi", "Energy", "PdgId"}){
    branchnames.push_back(l1+var);
    branchnames.push_back(l2+var);
    branchnames.push_back(l3+var);
    branchnames.push_back(l4+var);
  }
  _ntuple->SetBranchStatus("*", 0);
  for (std::string &branchname : branchnames){
    _ntuple->SetBranchStatus(branchname.c_str(), 1);
  }
  if (_isMC){
    _ntuple->SetBranchStatus("genWeight", 1);
    _ntuple->SetBranchStatus("nTruePU", 1);
  }
  _ntuple->SetBranchStatus((l1 + "_" + l2 + "_Mass").c_str(), 1);
  _ntuple->SetBranchStatus((l3 + "_" + l4 + "_Mass").c_str(), 1);

  // Define cutstring
  std::string cutstring = "";
  for (std::string lep : {l1, l2, l3, l4}){
    cutstring += lep + "ZZTightIDNoVtx && ";
    _ntuple->SetBranchStatus((lep + "ZZTightIDNoVtx").c_str(), 1);
    if (lep[0] == 'e'){
      cutstring += lep + "SIP3D < 4.0 && "; 
    }
    _ntuple->SetBranchStatus((lep + "SIP3D").c_str(), 1);
    _ntuple->SetBranchStatus((lep + "ZZIsoPass").c_str(), 1);
    _ntuple->SetBranchStatus((lep + "ZZIso").c_str(), 1);
  }
  cutstring = cutstring.substr(0, cutstring.size()-4);
  if (_nocuts) cutstring = "";
  else if (_verbose) std::cout << "- Applying cuts..." << std::endl;
  
  // Create trees
  if (_verbose) std::cout << "- Copying tree..." << std::endl;
  TTree *outtree = _ntuple->CopyTree(cutstring.c_str());

  // Add summedWeights
  if (_isMC){
    if (_verbose) std::cout << "- Adding summedWeights branch..." << std::endl;
    Float_t summedWeights = *ROOT::RDataFrame(*_meta).Sum("summedWeights");
    TBranch *b_summedWeights = outtree->Branch("summedWeights", &summedWeights, "summedWeights/F");
    for (Long64_t i=0; i<outtree->GetEntries(); i++)
      b_summedWeights->Fill();
  }
  
  // Writing
  if (_verbose) std::cout << "- Writing..." << std::endl;
  outfile->cd();
  outfile->rmdir(_channel.c_str());
  TDirectory *subdir = outfile->mkdir(_channel.c_str());
  if (subdir != nullptr) subdir->cd();
  outtree->Write("ntuple");

  outfile->Close();
  if (_verbose) std::cout << "- Done." << std::endl;
}

#ifndef __CLING__
#include "interface/argparse.h"

int main(int nargs, char *argv[]){
  auto parser = argparse::ArgumentParser(nargs, argv)
    .formatter_class(argparse::HelpFormatter::ArgumentDefaults);

  parser.add_argument<bool>("-v", "--verbose").def("false")
    .help("if true, script will be more verbose");
  parser.add_argument<bool>("--nocuts").def("false")
    .help("if true, don't apply any cuts to data");
  parser.add_argument<bool>("--mc").def("false")
    .help("if true, add additional MC branches to slimmed ntuple");
  parser.add_argument("-c", "--channels")
    .help("process only the specified comma-separated channels, otherwise all. options: eeee, eemm, mmmm");
  parser.add_argument<bool>("-r", "--recreate").def("false")
    .help("if true, the output file will be recreated for the given channel(s)");
  parser.add_argument("-f", "--filename")
    .help("use the given file for the list of datasets instead of searching by label");
  parser.add_argument("label").help("name for input dataset, as seen in 'inputs/*.dat");

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

  std::string filename = args["filename"].is_none()?
    (std::string)"inputs/" + args["label"].str() + ".dat" : args["filename"];

  int index=0;
  for (std::string channel : channels){
    ZZSlimmer l(args["label"].c_str(), channel.c_str());
    l.SetVerbose(args["verbose"]);
    l.AddFromFile(filename.c_str());
    l.SetMC(args["mc"]);
    l.SetNoCuts(args["nocuts"]);
    if ((index++)==0 && args["recreate"].is_true()) l.SetMode("recreate");
    l.Slim();
  }
  
  return 0;
}
#endif
