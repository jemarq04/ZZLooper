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
    void Loop();
  private:
    std::string _label1, _label2;
    std::string _mode = "update";
};

ZZSlimmer::ZZSlimmer(const char *name, const char *channel, const char *filenames)
    : ZZSlimmerBase(name, channel, filenames){
};

ZZSlimmer::~ZZSlimmer(){}

void ZZSlimmer::Loop(){
  // Setup
  TFile *outfile = new TFile(("slimmed/" + _name + ".root").c_str(), _mode.c_str());
  if (!outfile->IsOpen()) return;
  outfile->cd();
	
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
  std::vector<std::string> branchnames = {"Mass"};
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
  if (_isMC) _ntuple->SetBranchStatus("genWeight", 1);
  _ntuple->SetBranchStatus((l1 + "_" + l2 + "_Mass").c_str(), 1);
  _ntuple->SetBranchStatus((l3 + "_" + l4 + "_Mass").c_str(), 1);

  // Define cutstring
  std::string cutstring = "";
  for (std::string lep : {l1, l2, l3, l4}){
    cutstring += lep + "ZZTightIDNoVtx && ";
    _ntuple->SetBranchStatus((lep + "ZZTightIDNoVtx").c_str(), 1);
    if (lep[0] == 'e'){
      cutstring += lep + "SIP3D < 4.0 && "; 
      _ntuple->SetBranchStatus((lep + "SIP3D").c_str(), 1);
    }
  }
  cutstring = cutstring.substr(0, cutstring.size()-4);
  
  // Create trees
  TTree *outtree = _ntuple->CopyTree(cutstring.c_str());

  // Add summedWeights
  if (_isMC){
    Float_t summedWeights = *ROOT::RDataFrame(*_meta).Sum("summedWeights");
    TBranch *b_summedWeights = outtree->Branch("summedWeights", &summedWeights, "summedWeights/F");
    for (Long64_t i=0; i<outtree->GetEntries(); i++)
      b_summedWeights->Fill();
  }
  
  // Writing
  outfile->cd();
  outfile->rmdir(_channel.c_str());
  TDirectory *subdir = outfile->mkdir(_channel.c_str());
  if (subdir != nullptr) subdir->cd();
  outtree->Write("ntuple");

  outfile->Close();
}

#ifndef __CLING__
#include "interface/argparse.h"

int main(int nargs, char *argv[]){
	auto parser = argparse::ArgumentParser(nargs, argv)
		.formatter_class(argparse::HelpFormatter::ArgumentDefaults);

  parser.add_argument<bool>("--mc").def("false")
    .help("if true, add additional MC branches to slimmed ntuple");
  parser.add_argument("-c", "--channel").choices("eeee,eemm,mmmm")
    .help("process only the specified channel, otherwise all");
  parser.add_argument("-m", "--mode").def("UPDATE")
    .help("option for creating the output file");
  parser.add_argument("-f", "--filename")
    .help("use the given file for the list of datasets instead of searching by label");
  parser.add_argument("label").help("name for input dataset, as seen in 'inputs/*.dat");

  auto args = parser.parse_args();

  std::string filename = args["filename"].is_none()?
    (std::string)"inputs/" + args["label"].str() + ".dat" : args["filename"];
  std::cout << "Getting list of datasets from " << filename << std::endl;
  for (std::string channel : {"eeee", "eemm", "mmmm"}){
    if (!args["channel"].is_none() && args["channel"] != channel)
      continue;
    std::cout << "Slimming " << channel << " channel..." << std::endl;
    ZZSlimmer l(args["label"].c_str(), channel.c_str());
    l.AddFromFile(filename.c_str());
    if (args["mc"].is_true()) l.SetMC();
    l.SetMode(args["mode"]);
    l.Loop();
  }
	
  return 0;
}
#endif
