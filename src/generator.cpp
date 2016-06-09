#include <iostream>
#include <vector>

#include <sxmc/generator.h>
#include <sxmc/signals.h>

std::vector<float>
make_fake_dataset(std::vector<Signal>& signals,
                  std::vector<Systematic>& systematics,
                  std::vector<Observable>& observables,
                  std::vector<double> params, bool poisson) {
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
    observed[i] = \
      dynamic_cast<pdfz::EvalHist*>(signals[i].histogram)->RandomSample(
        events, params[i], syst_vals, upper, lower, poisson,
        signals[i].dataset);

    std::cout << "make_fake_dataset: " << signals[i].name << ": "
              << observed[i] << " events "
              << "(" << params[i] << " expected)" << std::endl;
  }

  return events;
}

