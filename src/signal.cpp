#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <sxmc/signal.h>
#include <sxmc/pdfz.h>
#include <sxmc/ttree_io.h>

Signal::Signal(std::string _name, std::string _title,
               std::string _filename, unsigned _dataset,
               Source _source, double _nexpected,
               std::vector<std::string>& sample_fields,
               std::vector<Observable>& observables,
               std::vector<Observable>& cuts,
               std::vector<Systematic>& systematics)
    : name(_name), title(_title), filename(_filename), dataset(_dataset),
      source(_source), nexpected(_nexpected) {
  std::vector<float> data;
  std::vector<unsigned int> rank;
  std::vector<std::string> ttree_fields;

  int rc = \
    sxmc::io::read_float_vector_ttree(filename, data, rank, ttree_fields);
  assert(rc >= 0);

  this->n_mc = rank[0];
  std::vector<float> samples(this->n_mc * sample_fields.size());

  // If user provided a scale factor for MC generation rather than a rate,
  // nexpected is set negative.
  if (this->nexpected < 0) {
    this->nexpected *= -1.0 * n_mc;
  }

  read_dataset_to_samples(samples, data, this->dataset,
                          sample_fields, ttree_fields, cuts);

  // Build the histogram evaluator
  build_pdfz(samples, sample_fields.size(), observables, systematics);

  // Evaluate histogram at mean of systematics to see how many
  // of our samples fall within our observable min and max limits
  //set_efficiency(systematics);

  for (size_t i=0; i<systematics.size(); i++) {
    this->systematic_names.push_back(systematics[i].name);
  }
}


void Signal::read_dataset_to_samples(std::vector<float>& samples,
                                     std::vector<float>& dataset,
                                     unsigned dataset_id,
                                     std::vector<std::string>& sample_fields,
                                     std::vector<std::string>& dataset_fields,
                                     std::vector<Observable>& cuts) {
  // Build a lookup table for cuts
  std::vector<bool> field_has_cut(dataset_fields.size(), false);
  std::vector<double> field_cut_lower(dataset_fields.size());
  std::vector<double> field_cut_upper(dataset_fields.size());
  for (size_t i=0; i<dataset_fields.size(); i++) {
    for (size_t j=0; j<cuts.size(); j++) {
      if (cuts[j].field == dataset_fields[i]) {
        field_has_cut[i] = true;
        field_cut_lower[i] = cuts[j].lower;
        field_cut_upper[i] = cuts[j].upper;
      }
    }
  }

  // Build a map from sample array index to dataset array index
  // The last sample field is reserved for the dataset id tag
  std::vector<size_t> sample_to_dataset_map;
  for (size_t i=0; i<sample_fields.size()-1; i++) {
    size_t index = (std::find(dataset_fields.begin(), dataset_fields.end(),
                              sample_fields[i]) -
                    dataset_fields.begin());
    sample_to_dataset_map.push_back(index);
  }

  // Loop over events and add them to the sample array
  size_t ndata = static_cast<size_t>(dataset.size() / dataset_fields.size());
  size_t sample_index = 0;
  for (size_t i=0; i<ndata; i++) {
    bool event_valid = true;

    // Apply cuts
    for (size_t j=0; j<dataset_fields.size(); j++){
      float data = dataset[i * dataset_fields.size() + j];
      if (field_has_cut[j] &&
          (data < field_cut_lower[j] || data > field_cut_upper[j])) {
        event_valid = false;
        break;
      }
    }

    if (!event_valid) {
      continue;
    }

    // Copy event data from dataset to sample array
    for (size_t j=0; j<sample_fields.size()-1; j++) {
      float data = dataset[i * dataset_fields.size() + sample_to_dataset_map[j]];
      samples[sample_index * sample_fields.size() + j] = data;
    }
    samples[sample_index * sample_fields.size() + sample_fields.size() - 1] = dataset_id;
    sample_index++;
  }
  samples.resize(sample_index * sample_fields.size());
}


void Signal::build_pdfz(std::vector<float> &samples, int nfields,
                        std::vector<Observable>& observables,
                        std::vector<Systematic>& systematics) {
  // Build bin and limit arrays
  std::vector<double> lower(observables.size());
  std::vector<double> upper(observables.size());
  std::vector<int> nbins(observables.size());
  for (size_t i=0; i<(size_t) nfields; i++) {
    for (size_t j=0; j<observables.size(); j++) {
      if (observables[j].field_index == i) {
        lower[i] = observables[j].lower;
        upper[i] = observables[j].upper;
        nbins[i] = observables[j].bins;
        break;
      }
    }
  }

  // Build the histogram evaluator
  this->histogram = \
    new pdfz::EvalHist(samples, nfields, observables.size(),
                       lower, upper, nbins, this->dataset);

  for (size_t i=0; i<systematics.size(); i++) {
    Systematic* syst = &systematics[i];

    // Indices for these systematic parameters
    hemi::Array<short>* pars = new hemi::Array<short>(syst->npars, true);
    for (unsigned i=0; i<syst->pidx.size(); i++) {
      pars->writeOnlyHostPtr()[i] = syst->pidx[i];
    }

    size_t o_field = syst->observable_field_index;
    size_t t_field = syst->truth_field_index;

    if (syst->type == pdfz::Systematic::SHIFT) {
      this->histogram->AddSystematic(
        pdfz::ShiftSystematic(o_field, pars));
    }
    else if (syst->type == pdfz::Systematic::SCALE) {
      this->histogram->AddSystematic(
        pdfz::ScaleSystematic(o_field, pars));
    }
    else if (syst->type == pdfz::Systematic::CTSCALE) {
      this->histogram->AddSystematic(
        pdfz::CosThetaScaleSystematic(o_field, pars));
    }
    else if (syst->type == pdfz::Systematic::RESOLUTION_SCALE) {
      this->histogram->AddSystematic(
        pdfz::ResolutionScaleSystematic(o_field, t_field, pars));
    }
    else {
      std::cerr << "Signal::build_pdfz: Unknown systematic ID "
                << (int)syst->type << std::endl;
      assert(false);
    }
  }
}


void Signal::print() const {
  std::cout << "  " << this->name << std::endl
    << "    Title: \"" << this->title << "\"" << std::endl
    << "    Filename: " << this->filename << std::endl
    << "    Dataset: " << this->dataset << std::endl
    << "    Source: " << this->source.name << std::endl
    << "    Nexpected: " << this->nexpected << std::endl
    << "    Nmc: " << this->n_mc << std::endl
    << "    Systematics: ";
  for (size_t i=0; i<this->systematic_names.size(); i++) {
    std::cout << this->systematic_names[i] << " ";
  }
  std::cout << std::endl;
}

