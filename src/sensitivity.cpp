#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>
#include <TFile.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom.h>
#include "config.h"
#include "generator.h"
#include "nll.h"
#include "mcmc.h"

// find the x value that a fraction 1-CL/2 of the distribution falls above
double upper_limit(TH1F* h, double cl=0.682) {
  double tail = (1.0 - cl) / 2;
  double integral = 0;
  int thisbin = h->GetNbinsX();
  while (integral < tail * h->Integral()) {
    integral = h->Integral(thisbin--, h->GetNbinsX());
  }
  return h->GetBinLowEdge(thisbin);
}

// run an ensemble of independent fake experiments, fit with mcmc, and
// tabulate background-fluctuation sensitivity
TH1F* ensemble(std::vector<Signal>& signals, Range<float>& e_range,
               Range<float>& r_range, unsigned steps,
               float burnin_fraction, std::string signal_name,
               float confidence, unsigned nexperiments) {
  TH1F* limits = new TH1F("limits", "Background-fluctuation sensitivity",
                          100, 0, 50);

  for (unsigned i=0; i<nexperiments; i++) {
    std::cout << "Experiment " << i << " / " << nexperiments << std::endl;
    // make fake data
    FakeDataGenerator gen(signals, e_range, r_range);
    float* norms = new float[signals.size()];
    for (size_t i=0; i<signals.size(); i++) {
      norms[i] = signals[i].nexpected;
    }
    TNtuple* data = gen(norms);

    // run mcmc
    TNtuple* lspace = mcmc(signals, data, steps, burnin_fraction);

    // calculate signal sensitivity
    TH1F hproj("hproj", "hproj", 1000, 0, 100);
    lspace->Draw((signal_name + ">>hproj").c_str(), "", "goff");
    double limit = upper_limit(&hproj, confidence);
    limits->Fill(limit);

    TCanvas c1;
    hproj.Rebin(10);
    hproj.Draw();
    c1.SaveAs("lproj.pdf");
    TFile f("lspace.root", "recreate");
    lspace->Write();
    f.Close();

    delete data;
    delete lspace;
    delete[] norms;
  }

  return limits;
}


int main(int argc, char* argv[]) {
  assert(argc == 2);
  gRandom->SetSeed(0);

  // load configuration from json file
  std::string config_filename = std::string(argv[1]);
  FitConfig fc(config_filename);
  fc.print();

  // run ensemble
  TH1F* limits = ensemble(fc.signals, fc.e_range, fc.r_range, fc.steps,
                         fc.burnin_fraction, fc.signal_name,
                         fc.confidence, fc.experiments);

  TCanvas c1;
  limits->Draw();
  c1.SaveAs((fc.output_file + "_limits.pdf").c_str());
  c1.SaveAs((fc.output_file + "_limits.C").c_str());

  double xq[1] = {0.5};
  double yq[1];
  limits->GetQuantiles(1, yq, xq);

  std::cout << "Median " << fc.confidence * 100 << "\% limit: " << yq[0]
            << std::endl;

  delete limits;

  return 0;
}

