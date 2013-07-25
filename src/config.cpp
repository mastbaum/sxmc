#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <json/value.h>
#include <json/reader.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include "signals.h"
#include "config.h"
#include "utils.h"
#include "pdfz.h"
#include "hdf5_io.h"

template <class T>
size_t get_index_with_append(std::vector<T>& v, T o) {
  size_t index = std::find(v.begin(), v.end(), o) - v.begin();
  if (index == v.size()) {
    v.push_back(o);
    return v.size() - 1;
  }
  return index;
}


FitConfig::FitConfig(std::string filename) {
  Json::Reader reader;
  Json::Value root;

  std::ifstream t(filename.c_str());
  std::string data((std::istreambuf_iterator<char>(t)),
                   std::istreambuf_iterator<char>());

  bool parse_ok = reader.parse(data, root);
  if (!parse_ok) {
    std::cout  << "FitConfig::FitConfig: JSON parse error:" << std::endl
               << reader.getFormattedErrorMessages();
    throw(1);
  }

  // experiment parameters
  const Json::Value experiment_params = root["experiment"];
  assert(experiment_params.isMember("live_time"));
  this->live_time = experiment_params["live_time"].asFloat();
  this->confidence = experiment_params["confidence"].asFloat();
  this->efficiency = experiment_params.get("efficiency", 1.0).asFloat();

  // pdf parameters
  const Json::Value pdf_params = root["pdfs"];
  std::vector<std::string> hdf5_fields;
  for (Json::Value::const_iterator it=pdf_params["hdf5_fields"].begin();
       it!=pdf_params["hdf5_fields"].end(); ++it) {
    hdf5_fields.push_back((*it).asString());
  }

  std::map<std::string, Observable> all_observables;
  for (Json::Value::const_iterator it=pdf_params["observables"].begin();
       it!=pdf_params["observables"].end(); ++it) {
    Json::Value o_json = pdf_params["observables"][it.key().asString()];
    Observable o;
    o.name = it.key().asString();
    o.title = o_json["title"].asString();
    o.field = o_json["field"].asString();
    o.bins = o_json["bins"].asInt();
    o.lower = o_json["min"].asFloat();
    o.upper = o_json["max"].asFloat();
    o.units = o_json["units"].asString();
    all_observables[it.key().asString()] = o;
  }

  std::map<std::string, Systematic> all_systematics;
  for (Json::Value::const_iterator it=pdf_params["systematics"].begin();
       it!=pdf_params["systematics"].end(); ++it) {
    Json::Value s_json = pdf_params["systematics"][it.key().asString()];
    Systematic s;
    s.name = it.key().asString();
    s.title = s_json["title"].asString();
    s.observable_field = s_json["observable_field"].asString();

    std::string type_string = s_json["type"].asString();
    if (type_string == "scale") {
        s.type = pdfz::Systematic::SCALE;
    }
    else if (type_string == "shift") {
        s.type = pdfz::Systematic::SHIFT;
    }
    else if (type_string == "resolution_scale") {
        s.type = pdfz::Systematic::RESOLUTION_SCALE;
        s.truth_field = s_json["truth_field"].asString();
    }
    else {
      std::cerr << "FitConfig::FitConfig: Unknown systematic type "
                << type_string << std::endl;
      throw(1);
    }

    s.mean = s_json["mean"].asFloat();
    s.sigma = s_json.get("sigma", 0.0).asFloat();
    s.fixed = s_json.get("fixed", false).asBool();

    all_systematics[it.key().asString()] = s;
  }

  // fit parameters
  const Json::Value fit_params = root["fit"];
  this->experiments = fit_params["experiments"].asInt();
  this->steps = fit_params["steps"].asInt();
  this->burnin_fraction = fit_params.get("burnin_fraction", 0.1).asFloat();
  this->signal_name = fit_params["signal_name"].asString();
  this->output_file = fit_params.get("output_file", "fit_spectrum").asString();
  this->debug_mode = fit_params.get("debug_mode", false).asBool();

  for (Json::Value::const_iterator it=fit_params["observables"].begin();
       it!=fit_params["observables"].end(); ++it) {
    this->observables.push_back(all_observables[(*it).asString()]);
  }

  for (Json::Value::const_iterator it=fit_params["cuts"].begin();
       it!=fit_params["cuts"].end(); ++it) {
    this->cuts.push_back(all_observables[(*it).asString()]);
  }

  for (Json::Value::const_iterator it=fit_params["systematics"].begin();
       it!=fit_params["systematics"].end(); ++it) {
    this->systematics.push_back(all_systematics[(*it).asString()]);
  }

  // signal parameters
  std::vector<std::string> signal_names;
  for (Json::Value::iterator it=fit_params["signals"].begin();
       it!=fit_params["signals"].end(); ++it) {
    signal_names.push_back((*it).asString());
  }

  const Json::Value all_signals = root["signals"];
  for (Json::Value::const_iterator it=all_signals.begin();
       it!=all_signals.end(); ++it) {
    if (std::find(signal_names.begin(), signal_names.end(),
                  it.key().asString()) == signal_names.end()) {
      continue;
    }

    const Json::Value signal_params = all_signals[it.key().asString()];

    Signal s;
    s.name = it.key().asString();
    s.title = signal_params.get("title", s.name).asString();
    s.sigma = \
      signal_params.get("sigma", 0.0).asFloat() *
      this->live_time * this->efficiency;
    s.nexpected = \
      signal_params["rate"].asFloat() * this->live_time * this->efficiency;

    // load data
    std::cout << "FitConfig::FitConfig: Loading data for "
              << s.name << std::endl;
    std::vector<std::string> filenames;
    for (Json::Value::const_iterator it=signal_params["files"].begin();
         it!=signal_params["files"].end(); ++it) {
      filenames.push_back((*it).asString());
    }

    // FIXME add offset and rank check to to read_float_vector_hdf5, to enable
    // chaining of multiple files
    std::vector<float> dataset;
    std::vector<unsigned> rank;
    for (size_t i=0; i<filenames.size(); i++) {
      int code = read_float_vector_hdf5(filenames[i], s.name,
                                        dataset, rank);
      assert(code >= 0);
    }

    // build sample dataset with the correct fields and ordering.
    // 1. build a unique list of fields
    // observables
    std::vector<size_t> sample_fields;
    for (size_t i=0; i<this->observables.size(); i++) {
      std::string field_name = this->observables[i].field;
      size_t field = (std::find(hdf5_fields.begin(), hdf5_fields.end(),
                                field_name) -
                      hdf5_fields.begin());
      assert(field < hdf5_fields.size());

      size_t index = get_index_with_append<size_t>(sample_fields, field);
      this->observables[i].field_index = index;
    }

    // systematics
    for (size_t i=0; i<this->systematics.size(); i++) {
      std::string field_name = this->systematics[i].observable_field;
      size_t field = (std::find(hdf5_fields.begin(), hdf5_fields.end(),
                                field_name) -
                      hdf5_fields.begin());
      assert(field < hdf5_fields.size());

      // systematic observable must be an observable
      size_t index = (std::find(sample_fields.begin(), sample_fields.end(),
                                field) -
                      sample_fields.begin());
      assert(index < sample_fields.size());

      this->systematics[i].observable_field_index = index;

      // add non-observable parameters
      if (this->systematics[i].type != pdfz::Systematic::RESOLUTION_SCALE) {
        continue;
      }

      field_name = this->systematics[i].truth_field;
      field = (std::find(hdf5_fields.begin(), hdf5_fields.end(), field_name) -
               hdf5_fields.begin());
      assert(field < hdf5_fields.size());

      index = get_index_with_append<size_t>(sample_fields, field);
      this->systematics[i].truth_field_index = index;
    }

    // 2. copy over the relevant data into sample array
    s.nevents = rank[0];
    std::vector<float> samples(s.nevents * sample_fields.size());

    for (size_t i=0; i<this->cuts.size(); i++) {
      std::string field_name = this->cuts[i].field;
      size_t index = (std::find(hdf5_fields.begin(), hdf5_fields.end(),
                                field_name) -
                      hdf5_fields.begin());

      this->cuts[i].field_index = index;  // index in hdf5 file
    }

    std::vector<bool> field_has_cut(hdf5_fields.size(), false);
    std::vector<double> field_cut_lower(hdf5_fields.size());
    std::vector<double> field_cut_upper(hdf5_fields.size());
    for (size_t i=0; i<hdf5_fields.size(); i++) {
      for (size_t j=0; j<this->cuts.size(); j++) {
        if (this->cuts[j].field_index == i) {
          field_has_cut[i] = true;
          field_cut_lower[i] = this->cuts[j].lower;
          field_cut_upper[i] = this->cuts[j].upper;
        }
      }
    }

    size_t sample_index = 0;
    for (size_t i=0; i<s.nevents; i++) {
      // apply cuts
      if (!this->cuts.empty()) {
        bool event_valid = true;
        for (size_t j=0; j<hdf5_fields.size(); j++) {
          float v = dataset[i * rank[1] + j];
          if (field_has_cut[j] &&
              (v < field_cut_lower[j] || v > field_cut_upper[j])) {
            event_valid = false;
            break;
          }
        }
        if (!event_valid) {
          continue;
        }
      }

      for (size_t j=0; j<sample_fields.size(); j++) {
        samples[sample_index * sample_fields.size() + j] = \
          dataset[i * rank[1] + sample_fields[j]];
      }
      sample_index++;
    }

    float years = 1.0 * s.nevents / (s.nexpected / this->live_time /
                  this->efficiency);

    std::cout << "FitConfig::FitConfig: Initializing PDF for " << s.name
              << " using " << s.nevents << " events (" << years << " y)"
              << std::endl;

    // perhaps we skipped some events due to cuts
    if (sample_index != s.nevents) {
      std::cout << "FitConfig::FitConfig: " << s.nevents - sample_index
                << " events cut" << std::endl;
    }
    samples.resize(sample_index * sample_fields.size());
    s.nexpected *= (1.0 * sample_index / s.nevents);
    s.nevents = sample_index;

    // build bin and limit arrays
    std::vector<double> lower(this->observables.size());
    std::vector<double> upper(this->observables.size());
    std::vector<int> nbins(this->observables.size());
    for (size_t i=0; i<sample_fields.size(); i++) {
      for (size_t j=0; j<this->observables.size(); j++) {
        if (this->observables[j].field_index == i) {
          lower[i] = this->observables[j].lower;
          upper[i] = this->observables[j].upper;
          nbins[i] = this->observables[j].bins;
          break;
        }
      }
    }

    // build the histogram evaluator
    s.histogram = new pdfz::EvalHist(samples, sample_fields.size(),
                                     this->observables.size(),
                                     lower, upper, nbins);

    for (size_t i=0; i<this->systematics.size(); i++) {
      Systematic* syst = &this->systematics[i];

      size_t o_field = syst->observable_field_index;
      size_t t_field = syst->truth_field_index;

      if (syst->type == pdfz::Systematic::SHIFT) {
        s.histogram->AddSystematic(pdfz::ShiftSystematic(o_field, i));
      }
      else if (syst->type == pdfz::Systematic::SCALE) {
        s.histogram->AddSystematic(pdfz::ScaleSystematic(o_field, i));
      }
      else if (syst->type == pdfz::Systematic::RESOLUTION_SCALE) {
        s.histogram->AddSystematic(pdfz::ResolutionScaleSystematic(o_field,
                                                                   t_field,
                                                                   i));
      }
      else {
        std::cerr << "FitConfig::FitConfig: Unknown systematic ID "
                  << (int)syst->type << std::endl;
        assert(false);
      }
    }

    this->signals.push_back(s);
  }

  // rescale expectations using fractions falling inside the pdfs
  hemi::Array<double> param_buffer(this->systematics.size(), true);
  param_buffer.writeOnlyHostPtr();
  for (size_t i=0; i<this->systematics.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = this->systematics[i].mean;                             
  }

  hemi::Array<unsigned> norms_buffer(signals.size(), true);
  norms_buffer.writeOnlyHostPtr();

  for (size_t i=0; i<this->signals.size(); i++) {
    pdfz::Eval* p = this->signals[i].histogram;
    p->SetNormalizationBuffer(&norms_buffer, i);
    p->SetParameterBuffer(&param_buffer);
    dynamic_cast<pdfz::EvalHist*>(p)->CreateHistogram();
    if (this->signals[i].nevents > 0) {
      this->signals[i].nexpected *= \
        (1.0 * norms_buffer.readOnlyHostPtr()[i] / this->signals[i].nevents);
    }
  }
}


void FitConfig::print() const {
  std::cout << "Fit:" << std::endl
            << "  Fake experiments: " << this->experiments << std::endl
            << "  MCMC steps: " << this->steps << std::endl
            << "  Burn-in fraction: " << this->burnin_fraction << std::endl
            << "  Signal name: " << this->signal_name << std::endl
            << "  Output plot: " << this->output_file << std::endl;

  std::cout << "Experiment:" << std::endl
            << "  Live time: " << this->live_time << " y" << std::endl
            << "  Confidence level: " << this->confidence << std::endl;

  std::cout << "Cuts:" << std::endl;
  for (std::vector<Observable>::const_iterator it=this->cuts.begin();
       it!=this->cuts.end(); ++it) {
    std::cout << "  " << it->name << std::endl
              << "    Title: \"" << it->title << "\"" << std::endl
              << "    Lower bound: "<< it->lower << std::endl
              << "    Upper bound: " << it->upper << std::endl;
  }

  std::cout << "Observables:" << std::endl;
  for (std::vector<Observable>::const_iterator it=this->observables.begin();
       it!=this->observables.end(); ++it) {
    std::cout << "  " << it->name << std::endl
              << "    Title: \"" << it->title << "\"" << std::endl
              << "    Lower bound: "<< it->lower << std::endl
              << "    Upper bound: " << it->upper << std::endl
              << "    Bins: " << it->bins << std::endl;
  }

  std::cout << "Signals:" << std::endl;
  for (std::vector<Signal>::const_iterator it=this->signals.begin();
       it!=this->signals.end(); ++it) {
    std::cout << "  " << it->name << std::endl
              << "    Title: \"" << it->title << "\"" << std::endl
              << "    Expectation: "<< it->nexpected << std::endl
              << "    Constraint: ";
    if (it->sigma != 0) {
      std::cout << it->sigma << std::endl;
    }
    else {
      std::cout << "none" << std::endl;
    }
  }

  if (this->systematics.size() > 0) {
    std::cout << "Systematics:" << std::endl;
    for (std::vector<Systematic>::const_iterator it=this->systematics.begin();
         it!=this->systematics.end(); ++it) {
      std::cout << "  " << it->name << std::endl
                << "    Title: \"" << it->title << "\"" << std::endl
                << "    Type: "<< it->type << std::endl
                << "    Observable: "<< it->observable_field << std::endl;
      if (it->type == pdfz::Systematic::RESOLUTION_SCALE) {
        std::cout << "    Truth: " << it->truth_field << std::endl;
      }
      std::cout << "    Mean: "<< it->mean << std::endl
                << "    Constraint: ";
      if (it->sigma != 0) {
        std::cout << it->sigma << std::endl;
      }
      else {
        std::cout << "none" << std::endl;
      }
      std::cout << "    Fixed: ";
      if (it->fixed) {
        std::cout << "yes" << std::endl;
      }
      else {
        std::cout << "no" << std::endl;
      }
    }
  }
}

