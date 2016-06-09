#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <TMath.h>
#include <json/value.h>

#include <sxmc/signals.h>
#include <sxmc/pdfz.h>
#include <sxmc/ttree_io.h>
#include <sxmc/pdfz.h>

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


void Signal::set_efficiency(std::vector<Systematic>& systematics) {
  size_t nsys = systematics.size();

  // Determine the total number of systematic parameters
  size_t npars = 0;
  for (size_t i=0; i<nsys; i++) {
    npars += systematics[i].npars;
  }

  // Allocate and fill the parameter buffer
  hemi::Array<double> param_buffer(npars, true);
  param_buffer.writeOnlyHostPtr();

  size_t k = 0;
  for (size_t i=0; i<nsys; i++) {
    for (size_t j=0; j<systematics[i].npars; j++) {
      param_buffer.writeOnlyHostPtr()[k++] = systematics[i].means[j];
    }
  }

  hemi::Array<unsigned> norms_buffer(1, true);
  norms_buffer.writeOnlyHostPtr();
  this->histogram->SetNormalizationBuffer(&norms_buffer);
  this->histogram->SetParameterBuffer(&param_buffer);
  dynamic_cast<pdfz::EvalHist*>(this->histogram)->EvalAsync(false);
  dynamic_cast<pdfz::EvalHist*>(this->histogram)->EvalFinished();

  // Efficiency is the number of events that make it into the histogram
  // over the number of physical events input.
  //
  // Note that this is dependent on the systematics, and for now it is
  // calculated with all systematics at means.
  this->nevents = norms_buffer.readOnlyHostPtr()[0];
  this->efficiency = this->nevents / (double) (this->n_mc);

  // nexpected = physical events expected * efficiency
  // sigma is fractional, does not scale
  this->nexpected *= this->efficiency;

  std::cout << "Signal::set_efficiency: "
            << this->nevents << "/" << this->n_mc << " events remain. "
            << "Total efficiency " << 100.0 * this->efficiency << "%"
            << std::endl;
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

Signal::Signal(std::string _name, std::string _title,
               std::string _filename, unsigned _dataset,
               double _nexpected, double _sigma, bool _fixed,
               std::vector<std::string>& sample_fields,
               std::vector<Observable>& observables,
               std::vector<Observable>& cuts,
               std::vector<Systematic>& systematics)
    : name(_name), title(_title), filename(_filename),
      dataset(_dataset), nexpected(_nexpected), sigma(_sigma), fixed(_fixed) {
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
  set_efficiency(systematics);
}


Observable::Observable(const std::string _name, const Json::Value& config)
      : name(_name) {
  assert(config.isMember("title"));
  this->title = config["title"].asString();

  assert(config.isMember("field"));
  this->field = config["field"].asString();

  assert(config.isMember("bins"));
  this->bins = config["bins"].asInt();

  assert(config.isMember("min"));
  this->lower = config["min"].asFloat();

  assert(config.isMember("max"));
  this->upper = config["max"].asFloat();

  this->units = config.get("units", "").asString();
  this->logscale = config.get("logscale", false).asBool();

  this->yrange.resize(2, -1);

  if (config.isMember("yrange")) {
    this->yrange[0] = config["yrange"][0].asFloat();
    this->yrange[1] = config["yrange"][1].asFloat();
  }
}


Systematic::Systematic(const std::string _name, const Json::Value& config)
      : name(_name) {
  assert(config.isMember("title"));
  this->title = config["title"].asString();

  assert(config.isMember("observable_field"));
  this->observable_field = config["observable_field"].asString();

  assert(config.isMember("type"));
  std::string type_string = config["type"].asString();

  if (type_string == "shift") {
    this->type = pdfz::Systematic::SHIFT;
  }
  else if (type_string == "scale") {
    this->type = pdfz::Systematic::SCALE;
  }
  else if (type_string == "ctscale") {
    this->type = pdfz::Systematic::CTSCALE;
  }
  else if (type_string == "resolution_scale") {
    this->type = pdfz::Systematic::RESOLUTION_SCALE;
    assert(config.isMember("truth_field"));
    this->truth_field = config["truth_field"].asString();
  }
  else {
    std::cerr << "FitConfig::load_pdf_systematics: Unknown systematic type "
              << type_string << std::endl;
    throw(1);
  }

  // Parameter means and standard deviations are an array of coefficients
  // in a power series expansion in the observable.
  assert(config.isMember("mean"));
  short npars = config["mean"].size();

  double* means = new double[npars];
  for (short j=0; j<npars; j++) {
    means[j] = config["mean"][j].asDouble();
  }

  double* sigmas = new double[npars];
  if (config.isMember("sigma")) {
    assert(config["sigma"].size() == (unsigned) npars);
    for (short j=0; j<npars; j++) {
      sigmas[j] = config["sigma"][j].asDouble();
    }
  }
  else {
    for (short j=0; j<npars; j++) {
      sigmas[j] = 0;
    }
  }

  this->npars = npars;
  this->means = means;
  this->sigmas = sigmas;

  // Fixing a systematic is all-or-nothing
  this->fixed = config.get("fixed", false).asBool();
}


//Rate::Rate(const std::string& _name, const Json::Value& params) : name(_name) {
//  // Must specify rate or scale, but not both
//  assert(params.isMember("rate") || params.isMember("scale"));
//  assert(!(params.isMember("rate") && params.isMember("scale")));
//
//  if (params.isMember("rate")) {
//    this->rate = params["rate"].asFloat();
//  }
//  else {
//    // (-) tell signals to scale by the total number of samples
//    this->rate = -1.0 / params["scale"].asFloat();
//  }
//
//  this->sigma = params.get("constraint", 0.0).asFloat();
//  this->fixed = params.get("fixed", false).asBool();
//}

