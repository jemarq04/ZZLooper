#include <iostream>
using namespace std;

void DumpNtupleBranches_comp(string filename){
  TFile *infile = TFile::Open(filename.c_str());
  vector<string> channels = {"eeee", "eemm", "mmmm"};
  vector<string> comps = {"1", "2"};
  for (const auto& channel : channels){
    TTree *tree = (TTree*)infile->Get((channel + "/ntuple").c_str());

    TObjArray *branches = tree->GetListOfBranches();
    TKey *key;
    TIter next(branches);
    ofstream fout(("comp/" + channel + "_vars.txt").c_str());
    while ((key = (TKey*)next())){
      TLeaf *leaf = tree->GetLeaf(key->GetName());
      if (!leaf){
        cout << "leaf not found: " << key->GetName() << endl;
        continue;
      }
      //cout << "Branch: " << key->GetName() << ", Type: " << leaf->GetTypeName() << endl;
      string tname = leaf->GetTypeName();
      bool isBaseType = (tname.size() >= 2)? 0==tname.compare(tname.size()-2, 2, "_t") : false;
      for (int i=0; i<2; i++){
        if (isBaseType)
          fout << "//" << tname << " " << key->GetName() << comps[i] << ";" << endl;
        else
          fout << "//std::" << tname << " *" << key->GetName() << comps[i] << " = 0;" << endl;
      }
      if (isBaseType) fout << "//" << tname << " " << key->GetName() << "EE;" << endl;
      else fout << "//std::" << tname << " *" << key->GetName() << "EE = 0;" << endl;
    }
    fout.close();
    next.Reset();
    fout.open(("comp/" + channel + "_addr.txt").c_str());
    while ((key = (TKey*)next())){
      for (int i=0; i<2; i++)
        fout << "//_ntuple" << comps[i] << "->SetBranchAddress(\"" << key->GetName() << "\", &" << key->GetName() << comps[i] << ");" << endl;
      fout << "//_ntupleEE->SetBranchAddress(\"" << key->GetName() << "\", &" << key->GetName() << "EE);" << endl;
    }
    fout.close();
  }

  infile->Close();
}
