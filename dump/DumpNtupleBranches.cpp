#include <iostream>
using namespace std;

void DumpNtupleBranches(string filename){
  TFile *infile = TFile::Open(filename.c_str());
  vector<string> channels = {"eeee", "eemm", "mmmm"};
  for (const auto& channel : channels){
    TTree *tree = (TTree*)infile->Get((channel + "/ntuple").c_str());

    TObjArray *branches = tree->GetListOfBranches();
    TKey *key;
    TIter next(branches);
    ofstream fout(("out/" + channel + "_vars.txt").c_str());
    while ((key = (TKey*)next())){
      TLeaf *leaf = tree->GetLeaf(key->GetName());
      if (!leaf){
        cout << "leaf not found: " << key->GetName() << endl;
        continue;
      }
      //cout << "Branch: " << key->GetName() << ", Type: " << leaf->GetTypeName() << endl;
      string tname = leaf->GetTypeName();
      bool isBaseType = (tname.size() >= 2)? 0==tname.compare(tname.size()-2, 2, "_t") : false;
      if (isBaseType)
        fout << "//" << tname << " " << key->GetName() << ";" << endl;
      else
        fout << "//std::" << tname << " *" << key->GetName() << " = 0;" << endl;

    }
    fout.close();
    next.Reset();
    fout.open(("out/" + channel + "_addr.txt").c_str());
    while ((key = (TKey*)next()))
      fout << "//_ntuple" << "->SetBranchAddress(\"" << key->GetName() << "\", &" << key->GetName() << ");" << endl;
    fout.close();
  }

  infile->Close();
}
