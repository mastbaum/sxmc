#include <iostream>
#include <assert.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TNtuple.h>
#include "generator.h"
#include "utils.h"

#define VERBOSE true

TNtuple* FakeDataGenerator::operator()(float* norms, bool poisson) {
  if (VERBOSE) {
    std::cout << "FakeDataGenerator::make_dataset: Generating dataset..."
              << std::endl;
  }
  TNtuple* nt = new TNtuple("tev", "events", "r:e");
  double r = 0;
  double e = 0;
  for (size_t i=0; i<this->pdfs.size(); i++) {
    float nexpected = norms[i];

    // shortcut: generate only events passing cuts
    if (this->pdfs[i]->IsA() == TH2F::Class()) {
      TH2F* ht = dynamic_cast<TH2F*>(this->pdfs[i]);
      int nobserved = nexpected;
      if (poisson && nexpected > 0) {
        nobserved = gRandom->Poisson(nexpected);
      }
      for (int i=0; i<nobserved; i++) {
        do {
          ht->GetRandom2(r, e);
        } while(e > e_range.max || e < e_range.min ||
                r > r_range.max || r < r_range.min);
        nt->Fill(r, e);
      }
      if (VERBOSE) {
        std::cout << "FakeDataGenerator::make_dataset: " << i << ": "
                  << nobserved << " events (" << nexpected << " expected)"
                  << std::endl;
      }
    }
    else if (this->pdfs[i]->IsA() == TH1D::Class()) {
      TH1D* ht = dynamic_cast<TH1D*>(this->pdfs[i]);
      int nobserved = nexpected;
      if (poisson && nexpected > 0) {
        nobserved = gRandom->Poisson(nexpected);
      }
      for (int i=0; i<nobserved; i++) {
        do {
          e = ht->GetRandom();
        } while(e > e_range.max || e < e_range.min);
        nt->Fill(r, e);
      }
      if (VERBOSE) {
        std::cout << "FakeDataGenerator::make_dataset: " << i << ": "
                  << nobserved << " events (" << nexpected << " expected)"
                  << std::endl;
      }
    }
    else {
      std::cerr << "FakeDataGenerator::make_dataset: Unknown histogram class "
                << pdfs[i]->ClassName() << std::endl;
      assert(false);
    }
  }

  return nt;
}

