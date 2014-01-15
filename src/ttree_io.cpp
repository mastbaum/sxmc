#include <TTree.h>

#include <sxmc/ttree_io.h>

#include <stdio.h>
#include <vector>
#include <string>
#include <unistd.h>
#include <iostream>
 
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TObjArray.h>
#include <TCollection.h>
#include <TBranch.h>
#include <TDataType.h>
 
int read_float_vector_ttree(const std::string &filename,
                            std::vector<float> &data, 
                            std::vector<unsigned int> &rank,
                            std::vector<std::string> &ttree_fields) {
  TFile *f;
  if (filename.compare(0,1,".") == 0){
    char filePath[256];
    getcwd(filePath, sizeof(filePath));
    printf("FILENAME: %s\n",(std::string(filePath)+filename.substr(1)).c_str());
    f = TFile::Open((std::string(filePath)+filename.substr(1)).c_str());
  }else{
    f = TFile::Open(filename.c_str());
  }

  if (!f){
    std::cout << "Problem opening " << filename << std::endl;
    return -1;
  }
    

  // get whatever ttree is in there
  TKey *k = (TKey*) f->GetListOfKeys()->First();
  TTree *t = (TTree*) f->Get(k->GetName());
  int num_entries = t->GetEntries();

  // we will now find all the branches in the ttree
  // and read out all those that can be stored as floats
  TObjArray *oa = t->GetListOfBranches();
  TIter next(oa);
  TBranch* temp;
  std::vector<TBranch*> temp_branches,branches;
  std::vector<EDataType> temp_types,types;
  std::vector<int> int_values;
  std::vector<float> float_values;
  std::vector<double> double_values;
  std::vector<my_bool> bool_values;
  while ( (temp = (TBranch *) next()) ) {
    temp_branches.push_back(temp);
    temp_types.push_back(kOther_t);
  }

  int first_file = 1;
  if (ttree_fields.size() > 0)
    first_file = 0;

  // Make a list of all the branches that we can convert to floats
  TClass *my_class;
  for (size_t i=0;i<temp_branches.size();i++){
    temp_branches[i]->GetExpectedType(my_class,temp_types[i]);

    if (temp_types[i] == 3 || temp_types[i] == 5 ||
        temp_types[i] == 8 || temp_types[i] == 18){
      types.push_back(temp_types[i]);
      branches.push_back(temp_branches[i]);
      
      if (first_file){
        ttree_fields.push_back(temp_branches[i]->GetName());
      }else{
        if (branches.back()->GetName() != ttree_fields[branches.size()-1]){
          // our second root file is not organized the same!
          return -1;
        }
      }

      int_values.push_back(-1);
      float_values.push_back(-1);
      double_values.push_back(-1);
      //bool_values.push_back(false);
      my_bool boolfix;
      boolfix.the_bool = false;
      bool_values.push_back(boolfix);
    }
  }

  for (size_t i=0;i<branches.size();i++){
    if (types[i] == 3){
      t->SetBranchAddress(branches[i]->GetName(),&int_values[i],&branches[i]);
    }else if (types[i] == 5){
      t->SetBranchAddress(branches[i]->GetName(),&float_values[i],&branches[i]);
    }else if (types[i] == 8){
      t->SetBranchAddress(branches[i]->GetName(),&double_values[i],&branches[i]);
    }else if (types[i] == 18){
      t->SetBranchAddress(branches[i]->GetName(),&(bool_values[i].the_bool),&branches[i]);
    }else{
      // something is wrong!
      return -1;
    }
  }

  if (first_file){
    rank.push_back(num_entries); 
    rank.push_back(branches.size());
  }else{
    rank[0] += num_entries;
  }
  int oldsize = data.size();
  data.resize(oldsize + num_entries*branches.size());

  // read in the values
  for (int i=0;i<num_entries;i++){
    t->GetEntry(i);
    for (size_t j=0;j<branches.size();j++){
      if (types[j] == 3){
        data[oldsize + i*branches.size()+j] = (float) int_values[j];
      }else if (types[j] == 5){
        data[oldsize + i*branches.size()+j] = (float) float_values[j];
      }else if (types[j] == 8){
        data[oldsize + i*branches.size()+j] = (float) double_values[j];
      }else if (types[j] == 18){
        data[oldsize + i*branches.size()+j] = (float) (bool_values[j].the_bool);
      }
    }
  }

  f->Close();
  return 0;
}
