#include <assert.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <unistd.h>
#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TObjArray.h>
#include <TCollection.h>
#include <TBranch.h>
#include <TDataType.h>

#include <sxmc/ttree_io.h>
 
namespace sxmc {
  namespace io {

int read_float_vector_ttree(const std::string &filename,
                            std::vector<float> &data, 
                            std::vector<unsigned int> &rank,
                            std::vector<std::string> &ttree_fields) {
  TFile* f;
  if (filename.compare(0, 1, ".") == 0) {
    char filePath[256];
    assert(getcwd(filePath, sizeof(filePath)) != NULL);
    printf("FILENAME: %s\n",
           (std::string(filePath)+filename.substr(1)).c_str());
    f = TFile::Open((std::string(filePath) + filename.substr(1)).c_str());
  }
  else {
    f = TFile::Open(filename.c_str());
  }

  if (!f) {
    std::cout << "read_float_vector_ttree: Problem opening "
              << filename << std::endl;
    return -1;
  }
    
  // get whatever ttree is in there
  TKey* k = dynamic_cast<TKey*>(f->GetListOfKeys()->First());
  TTree* t = dynamic_cast<TTree*>(f->Get(k->GetName()));
  int num_entries = t->GetEntries();

  // we will now find all the branches in the ttree
  // and read out all those that can be stored as floats
  TObjArray* oa = t->GetListOfBranches();
  TIter next(oa);
  TBranch* temp;
  std::vector<TBranch*> temp_branches, branches;
  std::vector<EDataType> temp_types, types;

  while ((temp = dynamic_cast<TBranch*>(next()))) {
    temp_branches.push_back(temp);
    temp_types.push_back(kOther_t);
  }

  bool first_file = (ttree_fields.size() == 0);

  // Make a list of all the branches that we can convert to floats
  TClass* cls;
  for (size_t i=0; i<temp_branches.size(); i++) {
    temp_branches[i]->GetExpectedType(cls, temp_types[i]);

    if (temp_types[i] == 3 || temp_types[i] == 5 ||
        temp_types[i] == 8 || temp_types[i] == 18) {
      types.push_back(temp_types[i]);
      branches.push_back(temp_branches[i]);
      
      if (first_file) {
        ttree_fields.push_back(temp_branches[i]->GetName());
      }
      else {
        if (branches.back()->GetName() != ttree_fields[branches.size()-1]) {
          // our second root file is not organized the same!
          return -1;
        }
      }
    }
  }

  int* int_values = new int[branches.size()];
  float* float_values = new float[branches.size()];
  double* double_values = new double[branches.size()];
  bool* bool_values = new bool[branches.size()];

  for (size_t i=0;i<branches.size();i++) {
    switch (types[i]) {
      case 3:
        t->SetBranchAddress(branches[i]->GetName(), &int_values[i],
                            &branches[i]);
        break;
      case 5:
        t->SetBranchAddress(branches[i]->GetName(), &float_values[i],
                            &branches[i]);
        break;
      case 8:
        t->SetBranchAddress(branches[i]->GetName(), &double_values[i],
                            &branches[i]);
        break;
      case 18:
        t->SetBranchAddress(branches[i]->GetName(), &(bool_values[i]),
                            &branches[i]);
        break;
      default:
        return -1;
    }
  }

  if (first_file) {
    rank.push_back(num_entries); 
    rank.push_back(branches.size());
  }
  else {
    rank[0] += num_entries;
  }

  int oldsize = data.size();
  data.resize(oldsize + num_entries * branches.size());

  // read in the values
  for (int i=0; i<num_entries; i++) {
    t->GetEntry(i);
    for (size_t j=0; j<branches.size(); j++) {
      if (types[j] == 3) {
        data[oldsize + i*branches.size()+j] = \
          static_cast<float>(int_values[j]);
      }
      else if (types[j] == 5) {
        data[oldsize + i*branches.size()+j] = \
          static_cast<float>(float_values[j]);
      }
      else if (types[j] == 8) {
        data[oldsize + i*branches.size()+j] = \
          static_cast<float>(double_values[j]);
      }
      else if (types[j] == 18) {
        data[oldsize + i*branches.size()+j] = \
          static_cast<float>(bool_values[j]);
      }
    }
  }

  f->Close();

  delete[] int_values;
  delete[] float_values;
  delete[] double_values;
  delete[] bool_values;

  return 0;
}

  }  // namespace io
}  // namespace sxmc

