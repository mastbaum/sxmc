#include <iostream>
#include <vector>

#include <sxmc/signal.h>
#include <sxmc/observable.h>
#include <sxmc/systematic.h>
#include <sxmc/generator.h>

std::vector<float>
make_fake_dataset(std::vector<Signal>& signals,
                  std::vector<Systematic>& systematics,
                  std::vector<Observable>& observables,
                  bool poisson,
                  std::string signal) {
  std::cout << "make_fake_dataset: Generating dataset..." << std::endl;

  std::vector<double> syst_vals;
  for (size_t i=0; i<systematics.size(); i++) {
    for (size_t j=0; j<systematics[i].npars; j++) {
      syst_vals.push_back(systematics[i].means[j]);
    }
  }

  std::vector<float> upper;
  std::vector<float> lower;
  for (size_t i=0; i<observables.size(); i++) {
    upper.push_back(observables[i].upper);
    lower.push_back(observables[i].lower);
  }

  std::vector<float> events;
  std::vector<unsigned> observed(signals.size());
  for (size_t i=0; i<signals.size(); i++) {
    // Skip the named signal
    if (signals[i].name.find(signal) != std::string::npos) {
      std::cout << "make_fake_dataset: Skipping signal "
                << signal << std::endl;
      continue;
    }

    double eff = signals[i].get_efficiency(systematics);
    double nevents = signals[i].nexpected * eff;

    observed[i] = \
      dynamic_cast<pdfz::EvalHist*>(signals[i].histogram)->RandomSample(
        events, nevents, syst_vals, upper, lower, poisson,
        signals[i].dataset);

    std::cout << "make_fake_dataset: " << signals[i].name << ": "
              << observed[i] << " events "
              << "(" << nevents << " expected, efficiency = "
              << eff << ")" << std::endl;
  }

  return events;
}

