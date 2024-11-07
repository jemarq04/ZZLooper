#ifndef ZZSLIMMERBASE_H
#define ZZSLIMMERBASE_H

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

#include "TH1.h"
#include "TChain.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"

//NOTE: Only accepts 2022 skimmed ntuples
class ZZSlimmerBase{
  public:
    ZZSlimmerBase(const char *name, const char *channel, const char *filenames="<NONE>");
    ~ZZSlimmerBase();

    void SetMC(bool isMC=true){_isMC = isMC;}

    bool AddFromFile(const char *filenames);
    virtual void Slim() = 0;
  protected:
    const double Z_MASS = 91.187; // GeV, given by PDG 2020

    TChain *_ntuple, *_meta;
    std::string _name;
    std::string _channel;

    bool _isMC = false;
};

ZZSlimmerBase::ZZSlimmerBase(const char *name, const char *channel, const char *filenames)
  : _name(name), _channel(channel){
  _ntuple = new TChain((_channel + "/ntuple").c_str());
  _meta = new TChain("metaInfo/metaInfo");
  if (strcmp(filenames, "<NONE>")){
    _ntuple->Add(filenames);
    _meta->Add(filenames);
    std::cout << "Added " << filenames << std::endl;
  }
}

ZZSlimmerBase::~ZZSlimmerBase(){
  if (_ntuple) delete _ntuple;
  if (_meta) delete _meta;
}

bool ZZSlimmerBase::AddFromFile(const char *filename){
  std::ifstream infile(filename);
  if (!infile.is_open()) return false;
  bool added = true;
  for (std::string line; getline(infile, line);){
    if (line.size() && line[0] == '#') continue;
    if (!_ntuple->AddFile(line.c_str()) || !_meta->AddFile(line.c_str())){
      added = false;
      break;
    }
  }
  infile.close();
  return added;
}

#endif//ZZSLIMMERBASE_H
