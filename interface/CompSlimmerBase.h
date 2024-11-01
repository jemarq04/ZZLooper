#ifndef COMPSLIMMERBASE_H
#define COMPSLIMMERBASE_H

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

#include "TH1.h"
#include "TChain.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"

//NOTE: Only accepts 2022 skimmed ntuples
class CompSlimmerBase{
  public:
    CompSlimmerBase(const char *name, const char *channel, const char *filenames1="<NONE>", const char *filenames2="<NONE>");
    ~CompSlimmerBase();

    bool AddFromFile1(const char *filename);
    bool AddFromFile2(const char *filename);

    void SetMC1(bool isMC=true){_isT1MC = isMC;}
    void SetMC2(bool isMC=true){_isT2MC = isMC;}
    void SetMC(bool isMC1=true, bool isMC2=true){SetMC1(isMC1); SetMC2(isMC2);}

    virtual void Loop() = 0;
  protected:
    void SetTitles(TH1 *hist, std::string x, std::string y="Number of Entries") const;
    double GetDeltaPhi(double phi1, double phi2) const;
    double GetPolCosTheta(const ROOT::Math::PtEtaPhiEVector &lep, const ROOT::Math::PtEtaPhiEVector &z) const;

    const double Z_MASS = 91.187; // GeV, given by PDG 2020

    TChain *_ntuple1, *_ntuple2;
    TChain *_meta1, *_meta2;
    std::string _name;
    std::string _channel;

    bool _isT1MC = false, _isT2MC = false;
};

CompSlimmerBase::CompSlimmerBase(const char *name, const char *channel, const char *filenames1, const char *filenames2) 
  : _name(name), _channel(channel){
  _ntuple1 = new TChain((_channel + "/ntuple").c_str());
  _ntuple2 = new TChain((_channel + "/ntuple").c_str());
  _meta1 = new TChain("metaInfo/metaInfo");
  _meta2 = new TChain("metaInfo/metaInfo");

  if (strcmp(filenames1, "<NONE>")){
    _ntuple1->Add(filenames1);
    _meta1->Add(filenames1);
    std::cout << "Added to first chain " << filenames1 << std::endl;
  }
  if (strcmp(filenames2, "<NONE>")){
    _ntuple2->Add(filenames2);
    _meta2->Add(filenames2);
    std::cout << "Added to second chain " << filenames2 << std::endl;
  }
}

CompSlimmerBase::~CompSlimmerBase(){
  if (_ntuple1) delete _ntuple1;
  if (_ntuple2) delete _ntuple2;
  if (_meta1) delete _meta1;
  if (_meta2) delete _meta2;
}

bool CompSlimmerBase::AddFromFile1(const char *filename){
  std::ifstream infile(filename);
  if (!infile.is_open()) return false;
  bool added = true;
  for(std::string line; getline(infile, line);){
    if (line.size() && line[0] == '#') continue;
    if (!_ntuple1->AddFile(line.c_str()) || !_meta1->AddFile(line.c_str())){
      added = false;
      break;
    }
  }
  infile.close();
  if (added) std::cout << "Added to first chain from file " << filename << std::endl;
  return added;
}

bool CompSlimmerBase::AddFromFile2(const char *filename){
  std::ifstream infile(filename);
  if (!infile.is_open()) return false;
  bool added = true;
  for(std::string line; getline(infile, line);){
    if (line.size() && line[0] == '#') continue;
    if (!_ntuple2->AddFile(line.c_str()) || !_meta2->AddFile(line.c_str())){
      added = false;
      break;
    }
  }
  infile.close();
  if (added) std::cout << "Added to second chain from file " << filename << std::endl;
  return added;
}

void CompSlimmerBase::SetTitles(TH1 *hist, std::string x, std::string y) const{
  if (hist->InheritsFrom("TH3")) return;
  hist->GetXaxis()->SetTitle(x.c_str());
  hist->GetYaxis()->SetTitle(y.c_str());
  //hist->GetYaxis()->SetTitleOffset(1.4);
}

double CompSlimmerBase::GetDeltaPhi(double phi1, double phi2) const{
  double deltaphi = phi1 - phi2;
  while (deltaphi < 0) deltaphi += (2 * M_PI);
  while (deltaphi > M_PI) deltaphi -= (2*M_PI);
  return deltaphi;
}

double CompSlimmerBase::GetPolCosTheta(const ROOT::Math::PtEtaPhiEVector &lp, const ROOT::Math::PtEtaPhiEVector &ln) const{
  ROOT::Math::PtEtaPhiEVector z = lp+ln;
  ROOT::Math::Boost boost(z.BoostToCM());
  ROOT::Math::PtEtaPhiEVector lep_boost = boost(ln);
  return lep_boost.Vect().Dot(z.Vect()) / std::sqrt( lep_boost.Vect().Mag2() * z.Vect().Mag2() );
}

#endif//COMPSLIMMERBASE_H
