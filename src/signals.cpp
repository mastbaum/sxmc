#include <algorithm>
#include <iostream>
#include <string>
#include <TMath.h>

#include <sxmc/signals.h>
#include <sxmc/pdfz.h>
#include <sxmc/ttree_io.h>

void Signal::read_dataset_to_samples(std::vector<float>& samples,
                                     std::vector<float>& dataset,
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
  std::vector<size_t> sample_to_dataset_map;
  for (size_t i=0; i<sample_fields.size(); i++) {
    size_t index = (std::find(dataset_fields.begin(), dataset_fields.end(),
                              sample_fields[i]) -
                    dataset_fields.begin());
    sample_to_dataset_map.push_back(index);
  }

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
    for (size_t j=0; j<sample_fields.size(); j++) {
      float data = dataset[i * dataset_fields.size() +
                           sample_to_dataset_map[j]];
      samples[sample_index * sample_fields.size() + j] = data;
    }
    sample_index++;
  }
  samples.resize(sample_index * sample_fields.size());
}


void Signal::apply_exclusions(std::vector<float>& samples,
                              std::vector<std::string>& sample_fields,
                              std::vector<int>& weights,
                              std::vector<Observable>& observables) {
  // Build a lookup table for excluded regions in the observables
  std::vector<bool> field_has_exclude(sample_fields.size(), false);
  std::vector<double> field_exclude_lower(sample_fields.size());
  std::vector<double> field_exclude_upper(sample_fields.size());
  for (size_t i=0; i<sample_fields.size(); i++) {
    for (size_t j=0; j<observables.size(); j++) {
      if (observables[j].field == sample_fields[i] &&
          observables[j].exclude) {
        field_has_exclude[i] = true;
        field_exclude_lower[i] = observables[j].exclude_min;
        field_exclude_upper[i] = observables[j].exclude_max;
      }
    }
  }

  size_t nsamples = static_cast<size_t>(samples.size() / sample_fields.size());
  size_t sample_index = 0;
  for (size_t i=0; i<nsamples; i++) {
    // Apply cuts
    unsigned excludes_total = 0;
    unsigned excludes_cut = 0;

    for (size_t j=0; j<sample_fields.size(); j++) {
      float v = samples[i * sample_fields.size() + j];

      // Use the union of excluded regions in observable ranges, i.e. only
      // cut an event if it is excluded in all observables
      if (field_has_exclude[j]) {
        excludes_total++;
        if (v >= field_exclude_lower[j] && v <= field_exclude_upper[j]) {
          excludes_cut++;
        }
      }
    }

    if (excludes_total > 0 && excludes_cut == excludes_total) {
      continue;
    }

    // Fill the samples array with passing events in place
    for (size_t j=0; j<sample_fields.size(); j++) {
      float data = samples[i * sample_fields.size() + j];
      samples[sample_index * sample_fields.size() + j] = data;
    }

    if (weights.size() > 0) {
      weights[sample_index] = weights[i];
    }

    sample_index++;
  }

  samples.resize(sample_index * sample_fields.size());
  weights.resize(sample_index);
}


// TODO: Tidy up
void Signal::set_efficiency(std::vector<Systematic> &systematics) {
  hemi::Array<double> param_buffer(systematics.size(), true);
  param_buffer.writeOnlyHostPtr();
  for (size_t i=0; i<systematics.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = systematics[i].mean;                             
  }
  hemi::Array<unsigned> norms_buffer(1, true);
  norms_buffer.writeOnlyHostPtr();
  this->histogram->SetNormalizationBuffer(&norms_buffer);
  this->histogram->SetParameterBuffer(&param_buffer);
  dynamic_cast<pdfz::EvalHist*>(this->histogram)->EvalAsync(false);
  dynamic_cast<pdfz::EvalHist*>(this->histogram)->EvalFinished();

  // efficiency is the number of events that make it into the histogram over the number of physical events input
  // note that this is dependent on the systematics, and for now it is calculated with all systematics at means 
  this->nevents = norms_buffer.readOnlyHostPtr()[0];
  this->efficiency = this->nevents / (double) (this->n_mc);
  // nexpected = physical events expected * efficiency
  this->nexpected *= this->efficiency;
  // uncertainty scales by efficiency as well
  this->sigma *= this->efficiency;

  std::cout << "Signal::Signal: " << this->nevents << " events remaining. " << this->n_mc-this->nevents << " events cut out of " << this->n_mc << " events. Total efficiency: " << this->efficiency << std::endl;
}


void Signal::build_pdfz(std::vector<float> &samples,
                        std::vector<int> &weights, int nfields,
                        std::vector<Observable> &observables,
                        std::vector<Systematic> &systematics) {
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
    new pdfz::EvalHist(samples, weights, nfields, observables.size(),
                       lower, upper, nbins);

  for (size_t i=0; i<systematics.size(); i++) {
    Systematic* syst = &systematics[i];

    size_t o_field = syst->observable_field_index;
    size_t t_field = syst->truth_field_index;

    if (syst->type == pdfz::Systematic::SHIFT) {
      this->histogram->AddSystematic(
        pdfz::ShiftSystematic(o_field, i));
    }
    else if (syst->type == pdfz::Systematic::SCALE) {
      this->histogram->AddSystematic(
        pdfz::ScaleSystematic(o_field, i));
    }
    else if (syst->type == pdfz::Systematic::RESOLUTION_SCALE) {
      this->histogram->AddSystematic(
        pdfz::ResolutionScaleSystematic(o_field, t_field, i));
    }
    else {
      std::cerr << "Signal::build_pdfz: Unknown systematic ID "
                << (int)syst->type << std::endl;
      assert(false);
    }
  }
}


Signal::Signal(std::string _name, std::string _title, float _nexpected,
               float _sigma, std::string _category,
               std::vector<std::string>& sample_fields,
               std::vector<Observable>& observables,
               std::vector<Observable>& cuts,
               std::vector<Systematic>& systematics,
               std::vector<std::string>& filenames)
    : name(_name), title(_title), category(_category), nexpected(_nexpected),
      sigma(_sigma), efficiency(1) {
  std::vector<float> dataset;
  std::vector<unsigned int> rank;
  std::vector<std::string> ttree_fields;

  for (size_t i=0; i<filenames.size(); i++) {
    int code = sxmc::io::read_float_vector_ttree(filenames[i], dataset,
                                                 rank, ttree_fields);
    assert(code >= 0);
  }

  this->n_mc = rank[0];
  std::vector<float> samples(this->n_mc * sample_fields.size());

  read_dataset_to_samples(samples, dataset, sample_fields, ttree_fields, cuts);
  apply_exclusions(samples, sample_fields, observables);

  // Create default weights
  std::vector<int> weights(samples.size() / sample_fields.size(), 1);

  // Build the histogram evaluator
  build_pdfz(samples, weights, sample_fields.size(), observables, systematics);

  // Evaluate histogram at mean of systematics to see how many
  // of our samples fall within our observable min and max limits
  set_efficiency(systematics);
}


Signal::Signal(std::string _name, std::string _title, float _nexpected,
               float _sigma, std::string _category,
               std::vector<Observable>& observables,
               std::vector<Observable>& cuts,
               std::vector<Systematic>& systematics,
               std::vector<float>& samples,
               std::vector<std::string>& sample_fields,
               std::vector<int> &weights)
    : name(_name), title(_title), category(_category), nexpected(_nexpected),
      sigma(_sigma), efficiency(1) {
  this->n_mc = 0;
  for (size_t i=0; i<weights.size(); i++)
    this->n_mc += weights[i];

  apply_exclusions(samples, sample_fields, weights, observables);

  // build the histogram evaluator
  build_pdfz(samples, weights, sample_fields.size(), observables, systematics);

  // Evaluate histogram at mean of systematics to see how many
  // of our samples fall within our observable min and max limits
  set_efficiency(systematics);
}

