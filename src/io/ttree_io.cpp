#include <unistd.h>
#include <cassert>
#include <cstdio>
#include <vector>
#include <string>
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

int read_float_vector_ttree(const std::string& filename,
                            std::vector<float>& data, 
                            std::vector<unsigned int>& rank,
                            std::vector<std::string>& ttree_fields) {
  TFile* f = NULL;
  if (filename.compare(0, 1, ".") == 0) {
    char filePath[512];
    assert(getcwd(filePath, sizeof(filePath)) != NULL);
    printf("FILENAME: %s\n",
           (std::string(filePath) + filename.substr(1)).c_str());
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
    
  // Get the first object in the file
  TKey* k = dynamic_cast<TKey*>(f->GetListOfKeys()->First());
  assert(k);
  TTree* t = dynamic_cast<TTree*>(f->Get(k->GetName()));
  assert(t);
  int num_entries = t->GetEntries();

  // Find all the branches in the TTree, read out all those that
  // can be stored as floats
  TObjArray* oa = t->GetListOfBranches();
  assert(oa);
  TIter branch_iter(oa);
  TBranch* temp;
  std::vector<TBranch*> temp_branches, branches;
  std::vector<EDataType> temp_types, types;

  while ((temp = dynamic_cast<TBranch*>(branch_iter()))) {
    temp_branches.push_back(temp);
    temp_types.push_back(kOther_t);
  }

  bool first_file = ttree_fields.empty();

  // Make a list of all the branches that we can convert to floats
  TClass* cls;
  for (size_t i=0; i<temp_branches.size(); i++) {
    temp_branches[i]->GetExpectedType(cls, temp_types[i]);

    if (temp_types[i] == kInt_t    || temp_types[i] == kFloat_t ||
        temp_types[i] == kDouble_t || temp_types[i] == kBool_t) {
      types.push_back(temp_types[i]);
      branches.push_back(temp_branches[i]);
      
      if (first_file) {
        ttree_fields.push_back(temp_branches[i]->GetName());
      }
      else {
        if (branches.back()->GetName() != ttree_fields[branches.size()-1]) {
          // Our second root file is not organized the same!
          std::cerr << "read_float_vector_ttree: Format mismatch between "
                       "input files!" << std::endl;
          return -1;
        }
      }
    }
  }

  int* int_values = new int[branches.size()];
  float* float_values = new float[branches.size()];
  double* double_values = new double[branches.size()];
  bool* bool_values = new bool[branches.size()];

  for (size_t i=0; i<branches.size(); i++) {
    switch (types[i]) {
      case kInt_t:
        t->SetBranchAddress(branches[i]->GetName(), &int_values[i],
                            &branches[i]);
        break;
      case kFloat_t:
        t->SetBranchAddress(branches[i]->GetName(), &float_values[i],
                            &branches[i]);
        break;
      case kDouble_t:
        t->SetBranchAddress(branches[i]->GetName(), &double_values[i],
                            &branches[i]);
        break;
      case kBool_t:
        t->SetBranchAddress(branches[i]->GetName(), &bool_values[i],
                            &branches[i]);
        break;
      default:
        std::cerr << "read_float_vector_ttree: Attempted to read a non-"
                     "convertible field." << std::endl;
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

  size_t oldsize = data.size();
  data.resize(oldsize + num_entries * branches.size());

  // Read the values into the data output vector
  for (int i=0; i<num_entries; i++) {
    t->GetEntry(i);
    for (size_t j=0; j<branches.size(); j++) {
      size_t index = oldsize + i * branches.size() + j;
      if (types[j] == kInt_t) {
        data[index] = static_cast<float>(int_values[j]);
      }
      else if (types[j] == kFloat_t) {
        data[index] = static_cast<float>(float_values[j]);
      }
      else if (types[j] == kDouble_t) {
        data[index] = static_cast<float>(double_values[j]);
      }
      else if (types[j] == kBool_t) {
        data[index] = static_cast<float>(bool_values[j]);
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

