#include <cassert>
#include <iostream>
#include <utility>
#include <fstream>
#include <streambuf>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <json/value.h>
#include <json/reader.h>

#include <sxmc/config.h>
#include <sxmc/ttree_io.h>
#include <sxmc/signals.h>
#include <sxmc/utils.h>

FitConfig::FitConfig(std::string filename) : data(NULL) {
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

  const Json::Value fit_params = root["fit"];
  const Json::Value obs_params = root["pdfs"]["observables"];
  const Json::Value sys_params = root["pdfs"]["systematics"];

  // Load general fit parameters
  assert(fit_params.isMember("nexperiments"));
  this->nexperiments = fit_params["nexperiments"].asInt();
  assert(this->nexperiments > 0);

  assert(fit_params.isMember("nsteps"));
  this->nsteps = fit_params["nsteps"].asInt();
  assert(this->nsteps > 0);

  this->burnin_fraction = fit_params.get("burnin_fraction", 0.1).asFloat();
  this->debug_mode = fit_params.get("debug_mode", false).asBool();
  this->output_prefix = fit_params.get("output_prefix", "lspace").asString();
  this->plots = fit_params.get("plots", true).asBool();
  this->seed = fit_params.get("seed", 0).asInt64();
  this->confidence = fit_params.get("confidence", 0.683).asFloat();

  this->signal_name = fit_params.get("signal_name", "").asString();

  // Observables
  for (Json::Value::const_iterator it=fit_params["observables"].begin();
       it!=fit_params["observables"].end(); ++it) {
    std::string name = (*it).asString();
    assert(obs_params.isMember(name));
    this->observables.push_back(Observable(name, obs_params[name]));
  }

  // Cuts
  for (Json::Value::const_iterator it=fit_params["cuts"].begin();
       it!=fit_params["cuts"].end(); ++it) {
    std::string name = (*it).asString();
    assert(obs_params.isMember(name));
    for (size_t j=0; j<this->observables.size(); j++) {
      assert(this->observables[j].name != name);
    }
    this->cuts.push_back(Observable(name, obs_params[name]));
  }

  // Systematics
  for (Json::Value::const_iterator it=fit_params["systematics"].begin();
       it!=fit_params["systematics"].end(); ++it) {
    std::string name = (*it).asString();
    assert(sys_params.isMember(name));
    this->systematics.push_back(Systematic(name, sys_params[name]));
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
      // Store the truth value in our sample so we can use it to manipulate
      // the observed value correctly when we adjust this systematic
      field_name = this->systematics[i].truth_field;
      this->systematics[i].truth_field_index = \
        get_index_with_append<std::string>(this->sample_fields, field_name);
    }
  }

  // Add dataset tag to sample fields
  this->sample_fields.push_back("DATASET");

  //// Load signals
  const Json::Value rate_params = root["rates"];
  const Json::Value sig_params = root["signals"];

  for (Json::Value::const_iterator it=fit_params["signals"].begin();
       it!=fit_params["signals"].end(); ++it) {
    std::string name = (*it).asString();

    assert(sig_params.isMember(name));
    Json::Value config = sig_params[name];

    std::cout << "FitConfig: Loading signal: " << name << std::endl;

    assert(config.isMember("dataset"));
    unsigned dataset = config["dataset"].asUInt();

    this->datasets.insert(dataset);

    assert((!config.isMember("rate") &&  config.isMember("scale")) ||
           ( config.isMember("rate") && !config.isMember("scale")));

    float sigma = config.get("constraint", 0.0).asFloat();
    float fixed = config.get("fixed", false).asBool();
    float nexpected = 0;
    if (config.isMember("rate")) {
      nexpected = config["rate"].asFloat();
    }
    else {
      // (-) tells Signal to scale by the total number of samples
      nexpected = -1.0 / config["scale"].asFloat();
    }

    assert(config.isMember("title"));
    std::string title = config["title"].asString();

    assert(config.isMember("filename"));
    std::string filename = config["filename"].asString();

    this->signals.push_back(
      Signal(name, title, filename, dataset,
             nexpected, sigma, fixed,
             this->sample_fields, this->observables,
             this->cuts, this->systematics));
  }

// DATA LOADING FIXME
//
//  this->data = NULL;  // Deferred to signal loading phase
//  
//  std::vector<std::string> data_signals;
//  if (experiment_params.isMember("data")) {
//    for (unsigned k=0; k<experiment_params["data"].size(); k++) {
//      std::string dsname = experiment_params["data"][k].asString();
//      assert(std::find(signal_names.begin(), signal_names.end(), dsname) ==
//             signal_names.end());
//      data_signals.push_back(dsname);
//    }
//    this->data = new std::vector<Signal>;
//  }
//
//  // Loop over all possible signals
//
//    // Check if this signal is being used as fit data (and not as a PDF!)
//    if (std::find(data_signals.begin(), data_signals.end(), name) !=
//                  data_signals.end()) {
//      const Json::Value signal_params = all_signals[name];
//      SignalParams sp(signal_params, 1.0);
//      std::cout << "FitConfig: Loading data: " << name << std::endl;
//
//      // All active observables are treated as cuts, in order to clip the data
//      // to the PDF boundaries.
//      std::vector<Observable> cc = this->observables;
//      for (size_t ii=0; ii<this->cuts.size(); ii++) {
//        cc.push_back(this->cuts[ii]);
//      }
//
//      std::vector<Systematic> syst;  // No systematics applied to data!
//
//      this->data->push_back(
//        Signal(name, sp.title, sp.nexpected, sp.sigma, sp.category,
//               this->sample_fields, this->observables, cc,
//               syst, sp.files));
//
//      continue;
//    }

}


void FitConfig::print() const {
  std::cout << "Fit:" << std::endl
    << "  Number of experiments: " << this->nexperiments << std::endl
    << "  MCMC steps: " << this->nsteps << std::endl
    << "  Burn-in fraction: " << this->burnin_fraction << std::endl
    << "  Random seed (0=random): " << this->seed << std::endl
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

