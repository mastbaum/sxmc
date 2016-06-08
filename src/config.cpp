#include <iostream>
#include <utility>
#include <fstream>
#include <streambuf>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <json/value.h>
#include <json/reader.h>

#include <sxmc/config.h>
#include <sxmc/ttree_io.h>
#include <sxmc/signals.h>
#include <sxmc/utils.h>

SignalParams::SignalParams(const Json::Value& params, float scale) {
  // Must specify rate or scale, but not both
  assert(params.isMember("rate") || params.isMember("scale"));
  assert(!(params.isMember("rate") && params.isMember("scale")));

  if (params.isMember("rate")) {
    nexpected = params["rate"].asFloat() * scale;
  }
  else {
    // (-) tell signals to scale by the total number of samples
    nexpected = -scale / params["scale"].asFloat();
  }

  sigma = params.get("constraint", 0.0).asFloat() * scale;
  fixed = params.get("fixed", false).asBool();
  title = params.get("title", "Other").asString();
  category = params.get("category", "").asString();

  assert(params.isMember("files"));
  for (Json::Value::const_iterator it=params["files"].begin();
       it!=params["files"].end(); ++it) {
    files.push_back((*it).asString());
  }
}


std::map<std::string, Observable>
FitConfig::get_pdf_observables(const Json::Value& pdf_params) {
  std::map<std::string, Observable> all_observables;
  for (Json::Value::const_iterator it=pdf_params["observables"].begin();
       it!=pdf_params["observables"].end(); ++it) {
    Json::Value o_json = pdf_params["observables"][it.key().asString()];
    Observable o;
    o.name = it.key().asString();
    assert(o_json.isMember("title"));
    o.title = o_json["title"].asString();
    assert(o_json.isMember("field"));
    o.field = o_json["field"].asString();
    assert(o_json.isMember("bins"));
    o.bins = o_json["bins"].asInt();
    assert(o_json.isMember("min"));
    o.lower = o_json["min"].asFloat();
    assert(o_json.isMember("max"));
    o.upper = o_json["max"].asFloat();
    assert(o_json.isMember("units"));
    o.units = o_json["units"].asString();
    o.logscale = o_json.get("logscale", false).asBool();
    o.yrange.resize(2, -1);

    if (o_json.isMember("yrange")) {
      o.yrange[0] = o_json["yrange"][0].asFloat();
      o.yrange[1] = o_json["yrange"][1].asFloat();
    }

    if (o_json.isMember("exclude")) {
      o.exclude = true;
      o.exclude_min = o_json["exclude"][0].asFloat();
      o.exclude_max = o_json["exclude"][1].asFloat();
    }
    else {
      o.exclude = false;
    }

    all_observables[it.key().asString()] = o;
  }

  return all_observables;
}


std::map<std::string, Systematic>
FitConfig::get_pdf_systematics(const Json::Value& pdf_params) {
  std::map<std::string, Systematic> all_systematics;
  for (Json::Value::const_iterator it=pdf_params["systematics"].begin();
       it!=pdf_params["systematics"].end(); ++it) {
    Json::Value s_json = pdf_params["systematics"][it.key().asString()];
    Systematic s;
    s.name = it.key().asString();
    assert(s_json.isMember("title"));
    s.title = s_json["title"].asString();
    assert(s_json.isMember("observable_field"));
    s.observable_field = s_json["observable_field"].asString();

    assert(s_json.isMember("type"));
    std::string type_string = s_json["type"].asString();

    if (type_string == "shift") {
        s.type = pdfz::Systematic::SHIFT;
    }
    else if (type_string == "scale") {
        s.type = pdfz::Systematic::SCALE;
    }
    else if (type_string == "ctscale") {
        s.type = pdfz::Systematic::CTSCALE;
    }
    else if (type_string == "resolution_scale") {
        s.type = pdfz::Systematic::RESOLUTION_SCALE;
        assert(s_json.isMember("truth_field"));
        s.truth_field = s_json["truth_field"].asString();
    }
    else {
      std::cerr << "FitConfig::FitConfig: Unknown systematic type "
                << type_string << std::endl;
      throw(1);
    }

    // Parameter means and standard deviations can be a number or an array,
    // where the latter becomes coefficients in a power series expansion in
    // the observable.
    assert(s_json.isMember("mean"));
    short npars = s_json["mean"].isArray() ? s_json["mean"].size() : 1;
    double* means = NULL;
    double* sigmas = NULL;

    if (s_json["mean"].isArray()) {
      short npars = s_json["mean"].size();
      means = new double[npars];
      for (short i=0; i<npars; i++) {
        means[i] = s_json["mean"][i].asDouble();
      }
    }
    else {
      means = new double[1];
      means[0] = s_json["mean"].asDouble();
    }

    if (s_json.isMember("sigma")) {
      // Constraints for all parameters or no parameters
      assert(s_json["sigma"].isArray() == s_json["mean"].isArray());

      if (s_json["sigma"].isArray()) {
        // Constraints on each of multiple parameters
        assert(s_json["sigma"].size() == (unsigned) npars);
        sigmas = new double[npars];
        for (short i=0; i<npars; i++) {
          sigmas[i] = s_json["sigma"][i].asDouble();
        }
      }
      else {
        // A constraint on a scalar parameter
        sigmas = new double[1];
        sigmas[0] = s_json["sigma"].asDouble();
      }
    }
    else {
      // No constraints specified, set them all to zero
      for (short i=0; i<npars; i++) {
        sigmas[i] = 0;
      }
    }

    s.npars = npars;
    s.means = means;
    s.sigmas = sigmas;

    // Fixing a systematic is all-or-nothing
    s.fixed = s_json.get("fixed", false).asBool();

    all_systematics[it.key().asString()] = s;
  }

  return all_systematics;
}


FitConfig::FitConfig(std::string filename) {
  // Read JSON configuration file
  Json::Reader reader;
  Json::Value root;

  std::ifstream t(filename.c_str());
  std::string data((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

  bool parse_ok = reader.parse(data, root);
  if (!parse_ok) {
    std::cerr  << "FitConfig::FitConfig: JSON parse error:" << std::endl
               << reader.getFormattedErrorMessages();
    throw(1);
  }

  //// Experiment parameters
  const Json::Value experiment_params = root["experiment"];

  assert(experiment_params.isMember("live_time"));
  this->live_time = experiment_params["live_time"].asFloat();

  assert(experiment_params.isMember("confidence"));
  this->confidence = experiment_params["confidence"].asFloat();

  this->efficiency_correction = \
    experiment_params.get("efficiency_correction", 1.0).asFloat();

  this->data = NULL;  // Deferred to signal loading phase

  //// PDF parameters
  const Json::Value pdf_params = root["pdfs"];

  std::map<std::string, Observable> all_observables = \
    get_pdf_observables(pdf_params);

  std::map<std::string, Systematic> all_systematics = \
    get_pdf_systematics(pdf_params);

  //// Fit parameters
  const Json::Value fit_params = root["fit"];

  assert(fit_params.isMember("experiments"));
  this->experiments = fit_params["experiments"].asInt();

  assert(fit_params.isMember("steps"));
  this->steps = fit_params["steps"].asInt();

  this->burnin_fraction = fit_params.get("burnin_fraction", 0.1).asFloat();
  this->debug_mode = fit_params.get("debug_mode", false).asBool();
  this->signal_name = fit_params.get("signal_name", "").asString();
  this->prefix = fit_params.get("prefix", "lspace").asString();
  this->plots = fit_params.get("plots", true).asBool();
  this->seed = fit_params.get("seed", 0).asInt64();

  // Find observables we want to fit for
  for (Json::Value::const_iterator it=fit_params["observables"].begin();
       it!=fit_params["observables"].end(); ++it) {
    this->observables.push_back(all_observables[(*it).asString()]);
  }

  // Find observables we want to cut on
  for (Json::Value::const_iterator it=fit_params["cuts"].begin();
       it!=fit_params["cuts"].end(); ++it) {
    this->cuts.push_back(all_observables[(*it).asString()]);
  }

  // Find systematics we want to use
  for (Json::Value::const_iterator it=fit_params["systematics"].begin();
       it!=fit_params["systematics"].end(); ++it) {
    this->systematics.push_back(all_systematics[(*it).asString()]);
  }

  // Determine the order of data fields in the sampled data.
  // sample_fields is a list of column titles we want in the sampled array,
  // and field_index is an index for sample_fields
  //
  // This is built in the following order:
  //
  //   1. Fields for observables
  //   2. Extra unobserved fields for systematics (e.g. true energy)
  //
  for (size_t i=0; i<this->observables.size(); i++) {
    std::string field_name = this->observables[i].field;
    this->observables[i].field_index = \
      get_index_with_append<std::string>(this->sample_fields, field_name);
  }

  // Add any additional observables used for systematics, and set the indices
  // of observables used in systematics
  for (size_t i=0; i<this->systematics.size(); i++) {
    std::string field_name = this->systematics[i].observable_field;

    // Systematic observable must be an observable or a cut
    size_t index = (std::find(this->sample_fields.begin(),
                              this->sample_fields.end(),
                              field_name) -
                    this->sample_fields.begin());
    assert(index < this->sample_fields.size());

    this->systematics[i].observable_field_index = index;

    // Add non-observable parameters
    if (this->systematics[i].type == pdfz::Systematic::RESOLUTION_SCALE) {
      // Store the truth value in our sample so we can use it to manipulate the
      // observed value correctly when we adjust this systematic
      field_name = this->systematics[i].truth_field;
      this->systematics[i].truth_field_index = \
        get_index_with_append<std::string>(this->sample_fields, field_name);
    }
  }

  //// Signal parameters
  std::vector<std::string> signal_names;
  for (Json::Value::iterator it=fit_params["signals"].begin();
       it!=fit_params["signals"].end(); ++it) {
    signal_names.push_back((*it).asString());
  }

  std::vector<std::string> data_signals;
  if (experiment_params.isMember("data")) {
    for (unsigned k=0; k<experiment_params["data"].size(); k++) {
      std::string dsname = experiment_params["data"][k].asString();
      assert(std::find(signal_names.begin(), signal_names.end(), dsname) ==
             signal_names.end());
      data_signals.push_back(dsname);
    }
    this->data = new std::vector<Signal>;
  }

  // Loop over all possible signals
  const Json::Value all_signals = root["signals"];
  for (Json::Value::const_iterator it=all_signals.begin();
       it!=all_signals.end(); ++it) {
    std::string name = it.key().asString();

    // Check if this signal is being used as fit data (and not as a PDF!)
    if (std::find(data_signals.begin(), data_signals.end(), name) !=
                  data_signals.end()) {
      const Json::Value signal_params = all_signals[name];
      SignalParams sp(signal_params, 1.0);
      std::cout << "FitConfig: Loading data: " << name << std::endl;

      // All active observables are treated as cuts, in order to clip the data
      // to the PDF boundaries.
      std::vector<Observable> cc = this->observables;
      for (size_t ii=0; ii<this->cuts.size(); ii++) {
        cc.push_back(this->cuts[ii]);
      }

      this->data->push_back(
        Signal(name, sp.title, sp.nexpected, sp.sigma, sp.category,
               this->sample_fields, this->observables, cc,
               this->systematics, sp.files));

      continue;
    }

    // Check if we want this signal
    if (std::find(signal_names.begin(), signal_names.end(), name) ==
        signal_names.end()) {
      continue;
    }

    std::cout << "FitConfig: Loading signal: " << name << std::endl;

    const Json::Value signal_params = all_signals[name];

    const float scale = this->live_time * this->efficiency_correction;

    SignalParams sp(signal_params, scale);

    this->signals.push_back(
      Signal(name, sp.title, sp.nexpected, sp.sigma, sp.category,
             this->sample_fields, this->observables, this->cuts,
             this->systematics, sp.files, sp.fixed));
  }
}


void FitConfig::print() const {
  std::cout << "Fit:" << std::endl
    << "  Fake experiments: " << this->experiments << std::endl
    << "  MCMC steps: " << this->steps << std::endl
    << "  Burn-in fraction: " << this->burnin_fraction << std::endl
    << "  Random seed (0=random): " << this->seed << std::endl;

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
      std::cout << "    Means: ";
      for (size_t i=0; i<it->npars; i++) {
        std::cout << "(" << i << ") " << it->means[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "    Constraints: ";
      for (size_t i=0; i<it->npars; i++) {
        std::cout << "(" << i << ") " << it->sigmas[i] << " ";
      }
      std::cout << std::endl;
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

