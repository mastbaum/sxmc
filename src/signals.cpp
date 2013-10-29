#include <algorithm>
#include <iostream>
#include <TMath.h>

#include <sxmc/signals.h>
#include <sxmc/pdfz.h>
#include <sxmc/hdf5_io.h>

//FIXME do we need to include non observed fields so they can be used for systematic effects?
Signal::Signal(std::string _name, std::string _title, float _nexpected, float _sigma, std::string _category,
    std::vector<Observable>& observables,
    std::vector<Systematic>& systematics,
    std::vector<float>& _samples, std::vector<int> &_weights) : name(_name), title(_title), category(_category), nexpected(_nexpected),
  sigma(_sigma), efficiency(1)
{

  assert(_samples.size() % observables.size() == 0);
  this->nevents = 0;
  for (size_t i=0;i<_weights.size();i++)
    this->nevents += _weights[i];

  std::vector<float> samples(_samples.size());
  std::vector<int> weights(_weights.size());

  std::vector<bool> field_has_exclude(observables.size(), false);
  std::vector<double> field_exclude_lower(observables.size());
  std::vector<double> field_exclude_upper(observables.size());
  for (size_t i=0; i<observables.size(); i++) {
    if (observables[i].exclude) {
      field_has_exclude[i] = true;
      field_exclude_lower[i] = observables[i].exclude_min;
      field_exclude_upper[i] = observables[i].exclude_max;
    }
  }



  size_t sample_index = 0;
  for (size_t i=0; i<_weights.size(); i++) {
    // apply cuts
    size_t excludes_total = 0;
    size_t excludes_cut = 0;

    // exclude any samples in our exclude region
    for (size_t j=0; j<observables.size(); j++) {
      float v = _samples[i * observables.size() + j];
      // union of excluded regions in observable ranges
      if (field_has_exclude[j]) {
        excludes_total++;
        if (v >= field_exclude_lower[j] && v <= field_exclude_upper[j]) {
          excludes_cut++;
        }
      }
    }
    if (excludes_total > 0 && excludes_cut == excludes_total)
      continue;

    for (size_t j=0; j<observables.size(); j++){
      samples[sample_index * observables.size() + j] = _samples[i * observables.size() + j];
    }
    weights[sample_index] = _weights[i];
    sample_index++;
  }

  // perhaps we skipped some events due to cuts
  if (sample_index != _weights.size()) {
    int weighted_nevents = 0;
    for (size_t i=0;i<sample_index;i++)
      weighted_nevents += weights[i];
    double rescale = 1.0 * weighted_nevents/ this->nevents;
    std::cout << "Signal::Signal: " << this->nevents - weighted_nevents <<
      " of " << this->nevents << " events cut (" << rescale << ")" << std::endl;
    samples.resize(sample_index * observables.size());
    weights.resize(sample_index);
    this->efficiency *= rescale;
    this->nexpected *= rescale;
    this->sigma *= rescale;
    this->nevents = weighted_nevents;
  }

  // build bin and limit arrays
  std::vector<double> lower(observables.size());
  std::vector<double> upper(observables.size());
  std::vector<int> nbins(observables.size());
  for (size_t i=0; i<observables.size(); i++) {
    lower[i] = observables[i].lower;
    upper[i] = observables[i].upper;
    nbins[i] = observables[i].bins;
  }

  // build the histogram evaluator
  this->histogram = new pdfz::EvalHist(samples, weights, observables.size(),
      observables.size(),
      lower, upper, nbins);

  for (size_t i=0; i<systematics.size(); i++) {
    Systematic* syst = &systematics[i];

    size_t o_field = syst->observable_field_index;
    size_t t_field = syst->truth_field_index;

    if (syst->type == pdfz::Systematic::SHIFT) {
      this->histogram->AddSystematic(pdfz::ShiftSystematic(o_field, i));
    }
    else if (syst->type == pdfz::Systematic::SCALE) {
      this->histogram->AddSystematic(pdfz::ScaleSystematic(o_field, i));
    }
    else if (syst->type == pdfz::Systematic::RESOLUTION_SCALE) {
      this->histogram->AddSystematic(pdfz::ResolutionScaleSystematic(o_field,
            t_field,
            i));
    }
    else {
      std::cerr << "Signal::Signal: Unknown systematic ID "
        << (int)syst->type << std::endl;
      assert(false);
    }
  }

  // rescale expectations using fractions falling inside the pdfs
  hemi::Array<double> param_buffer(systematics.size(), true);
  param_buffer.writeOnlyHostPtr();
  for (size_t i=0; i<systematics.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = systematics[i].mean;                             
  }

  hemi::Array<unsigned> norms_buffer(1, true);
  norms_buffer.writeOnlyHostPtr();

  this->histogram->SetNormalizationBuffer(&norms_buffer);
  this->histogram->SetParameterBuffer(&param_buffer);
  dynamic_cast<pdfz::EvalHist*>(this->histogram)->CreateHistogram();
  if (this->nevents > 0) {
    float rescale = (1.0 * norms_buffer.readOnlyHostPtr()[0] / this->nevents);
    this->efficiency *= rescale;
    this->nexpected *= rescale;
    this->sigma *= rescale;
    std::cout << "Signal::Signal: " << this->nevents - norms_buffer.readOnlyHostPtr()[0] <<
      " of " << this->nevents << " events outside limits (" << rescale << "). Total efficiency: " << this->efficiency << std::endl;
    /*
       if (this->signals[i].name == this->signal_name) {
       this->signal_eff *= rescale;
       std::cout << "FitConfig::FitConfig: Signal efficiency: "
       << this->signal_eff << std::endl;
       }
       */
  }
}



Signal::Signal(std::string _name, std::string _title, float _nexpected, float _sigma, std::string _category,
    std::vector<std::string>& hdf5_fields,
    std::vector<size_t>& sample_fields,
    std::vector<Observable>& observables,
    std::vector<Observable>& cuts,
    std::vector<Systematic>& systematics,
    std::vector<std::string>& filenames) : name(_name), title(_title), category(_category), nexpected(_nexpected),
  sigma(_sigma), efficiency(1)
{
  // FIXME add offset and rank check to to read_float_vector_hdf5, to enable
  // chaining of multiple files
  std::vector<float> dataset;
  std::vector<unsigned> rank;
  for (size_t i=0; i<filenames.size(); i++) {
    int code = read_float_vector_hdf5(filenames[i], this->name,
        dataset, rank);
    assert(code >= 0);
  }

  // 2. copy over the relevant data into sample array
  this->nevents = rank[0];
  std::vector<float> samples(this->nevents * sample_fields.size());

  std::vector<size_t> cuts_hdf5_index;
  for (size_t i=0; i<cuts.size(); i++) {
    std::string field_name = cuts[i].field;
    size_t index = (std::find(hdf5_fields.begin(), hdf5_fields.end(),
          field_name) -
        hdf5_fields.begin());

    cuts_hdf5_index.push_back(index);  // index in hdf5 file
  }

  std::vector<bool> field_has_cut(hdf5_fields.size(), false);
  std::vector<double> field_cut_lower(hdf5_fields.size());
  std::vector<double> field_cut_upper(hdf5_fields.size());
  for (size_t i=0; i<hdf5_fields.size(); i++) {
    for (size_t j=0; j<cuts.size(); j++) {
      if (cuts_hdf5_index[j] == i) {
        field_has_cut[i] = true;
        field_cut_lower[i] = cuts[j].lower;
        field_cut_upper[i] = cuts[j].upper;
      }
    }
  }

  std::vector<bool> field_has_exclude(hdf5_fields.size(), false);
  std::vector<double> field_exclude_lower(hdf5_fields.size());
  std::vector<double> field_exclude_upper(hdf5_fields.size());
  for (size_t i=0; i<hdf5_fields.size(); i++) {
    for (size_t j=0; j<observables.size(); j++) {
      std::string field_name = observables[j].field;
      size_t index = (std::find(hdf5_fields.begin(), hdf5_fields.end(),
            field_name) -
          hdf5_fields.begin());

      if (observables[j].field == hdf5_fields[i] &&
          observables[j].exclude) {
        field_has_exclude[i] = true;
        field_exclude_lower[i] = observables[j].exclude_min;
        field_exclude_upper[i] = observables[j].exclude_max;
      }
    }
  }

  // hack to do (r/r_av)^3 since it's not in the hdf5 files yet
  std::vector<bool> do_r3_transform(sample_fields.size(), false);
  for (size_t i=0; i<observables.size(); i++) {
    if (observables[i].name != "R_CUBED") {
      continue;
    }

    std::string field_name = observables[i].field;
    size_t field = (std::find(hdf5_fields.begin(), hdf5_fields.end(),
          field_name) -
        hdf5_fields.begin());
    assert(field < hdf5_fields.size());

    for (size_t j=0; j<sample_fields.size(); j++) {
      if (field == sample_fields[j]) {
        do_r3_transform[j] = true;
        std::cout << "Signal::Signal: Doing R^3 transformation on "
          << sample_fields[j] << std::endl;
      }
    }
  }

  assert(rank[1] == hdf5_fields.size());
  size_t sample_index = 0;
  for (size_t i=0; i<this->nevents; i++) {
    // apply cuts
    bool event_valid = true;
    size_t excludes_total = 0;
    size_t excludes_cut = 0;

    for (size_t j=0; j<hdf5_fields.size(); j++) {
      float v = dataset[i * rank[1] + j];

      // cut on unobserved observables
      if (field_has_cut[j] &&
          (v < field_cut_lower[j] || v > field_cut_upper[j])) {
        event_valid = false;
        break;
      }

      // union of excluded regions in observable ranges
      if (field_has_exclude[j]) {
        excludes_total++;
        if (v >= field_exclude_lower[j] && v <= field_exclude_upper[j]) {
          excludes_cut++;
        }
      }
    }

    if (!event_valid ||
        (excludes_total > 0 && excludes_cut == excludes_total)) {
      continue;
    }

    for (size_t j=0; j<sample_fields.size(); j++) {
      float data = dataset[i * rank[1] + sample_fields[j]];
      if (do_r3_transform[j]) {
        data = TMath::Power(data / 6005.0, 3);
      }
      samples[sample_index * sample_fields.size() + j] = data;

    }
    sample_index++;
  }

  // perhaps we skipped some events due to cuts
  if (sample_index != this->nevents) {
    std::cout << "Signal::Signal: " << this->nevents - sample_index <<
      " of " << this->nevents << " events cut (" << 1.0*sample_index/this->nevents << ")" << std::endl;
    samples.resize(sample_index * sample_fields.size());
    this->efficiency *= (1.0 * sample_index / this->nevents);
    this->nexpected *= (1.0 * sample_index / this->nevents);
    this->sigma *= (1.0 * sample_index / this->nevents);
    this->nevents = sample_index;
  }

  // build bin and limit arrays
  std::vector<double> lower(observables.size());
  std::vector<double> upper(observables.size());
  std::vector<int> nbins(observables.size());
  for (size_t i=0; i<sample_fields.size(); i++) {
    for (size_t j=0; j<observables.size(); j++) {
      if (observables[j].field_index == i) {
        lower[i] = observables[j].lower;
        upper[i] = observables[j].upper;
        nbins[i] = observables[j].bins;
        break;
      }
    }
  }

  // create default weights
  std::vector<int> weights(samples.size(),1);

  // build the histogram evaluator
  this->histogram = new pdfz::EvalHist(samples, weights, sample_fields.size(),
      observables.size(),
      lower, upper, nbins);

  for (size_t i=0; i<systematics.size(); i++) {
    Systematic* syst = &systematics[i];

    size_t o_field = syst->observable_field_index;
    size_t t_field = syst->truth_field_index;

    if (syst->type == pdfz::Systematic::SHIFT) {
      this->histogram->AddSystematic(pdfz::ShiftSystematic(o_field, i));
    }
    else if (syst->type == pdfz::Systematic::SCALE) {
      this->histogram->AddSystematic(pdfz::ScaleSystematic(o_field, i));
    }
    else if (syst->type == pdfz::Systematic::RESOLUTION_SCALE) {
      this->histogram->AddSystematic(pdfz::ResolutionScaleSystematic(o_field,
            t_field,
            i));
    }
    else {
      std::cerr << "Signal::Signal: Unknown systematic ID "
        << (int)syst->type << std::endl;
      assert(false);
    }
  }

  // rescale expectations using fractions falling inside the pdfs
  hemi::Array<double> param_buffer(systematics.size(), true);
  param_buffer.writeOnlyHostPtr();
  for (size_t i=0; i<systematics.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = systematics[i].mean;                             
  }

  hemi::Array<unsigned> norms_buffer(1, true);
  norms_buffer.writeOnlyHostPtr();

  this->histogram->SetNormalizationBuffer(&norms_buffer);
  this->histogram->SetParameterBuffer(&param_buffer);
  dynamic_cast<pdfz::EvalHist*>(this->histogram)->CreateHistogram();
  if (this->nevents > 0) {
    float rescale = (1.0 * norms_buffer.readOnlyHostPtr()[0] / this->nevents);
    this->efficiency *= rescale;
    this->nexpected *= rescale;
    this->sigma *= rescale;

    std::cout << "Signal::Signal: " << this->nevents - norms_buffer.readOnlyHostPtr()[0] <<
      " of " << this->nevents << " events outside limits (" << rescale << "). Total efficiency: " << this->efficiency << std::endl;
    samples.resize(sample_index * observables.size());
    this->efficiency *= (1.0 * sample_index / this->nevents);
    /*
       if (this->signals[i].name == this->signal_name) {
       this->signal_eff *= rescale;
       std::cout << "FitConfig::FitConfig: Signal efficiency: "
       << this->signal_eff << std::endl;
       }
       */
  }
}
