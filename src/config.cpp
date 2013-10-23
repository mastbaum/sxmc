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
#include <TMath.h>

#include <sxmc/signals.h>
#include <sxmc/config.h>
#include <sxmc/utils.h>

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
  for (Json::Value::const_iterator it=pdf_params["hdf5_fields"].begin();
       it!=pdf_params["hdf5_fields"].end(); ++it) {
    this->hdf5_fields.push_back((*it).asString());
  }

  // create list of possible observables
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

  // create list of possible systematics
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

  // find observables we want to fit for
  for (Json::Value::const_iterator it=fit_params["observables"].begin();
       it!=fit_params["observables"].end(); ++it) {
    this->observables.push_back(all_observables[(*it).asString()]);
  }

  // find observables we want to cut on
  for (Json::Value::const_iterator it=fit_params["cuts"].begin();
       it!=fit_params["cuts"].end(); ++it) {
    this->cuts.push_back(all_observables[(*it).asString()]);
  }

  // find systematics we want to use
  for (Json::Value::const_iterator it=fit_params["systematics"].begin();
       it!=fit_params["systematics"].end(); ++it) {
    this->systematics.push_back(all_systematics[(*it).asString()]);
  }

  // we are now going to determine the order of data fields in our
  // sampled data. sample_fields will map it to the hdf5 fields,
  // and field_index will map each observable/syst to the data field
  for (size_t i=0; i<this->observables.size(); i++) {
    std::string field_name = this->observables[i].field;
    size_t field = (std::find(this->hdf5_fields.begin(), this->hdf5_fields.end(),
          field_name) -
        this->hdf5_fields.begin());
    assert(field < this->hdf5_fields.size());

    size_t index = get_index_with_append<size_t>(this->sample_fields, field);
    this->observables[i].field_index = index;
  }

  // Add any additional observables used for systematics
  for (size_t i=0; i<this->systematics.size(); i++) {
    std::string field_name = this->systematics[i].observable_field;
    size_t field = (std::find(this->hdf5_fields.begin(), this->hdf5_fields.end(),
          field_name) -
        this->hdf5_fields.begin());
    assert(field < this->hdf5_fields.size());

    // systematic observable must be an observable
    size_t index = (std::find(this->sample_fields.begin(), this->sample_fields.end(),
          field) -
        this->sample_fields.begin());
    assert(index < this->sample_fields.size());

    this->systematics[i].observable_field_index = index;

    // add non-observable parameters
    if (this->systematics[i].type != pdfz::Systematic::RESOLUTION_SCALE) {
      continue;
    }

    field_name = this->systematics[i].truth_field;
    field = (std::find(this->hdf5_fields.begin(), this->hdf5_fields.end(), field_name) -
        this->hdf5_fields.begin());
    assert(field < this->hdf5_fields.size());

    index = get_index_with_append<size_t>(this->sample_fields, field);
    this->systematics[i].truth_field_index = index;
  }

  // signal parameters
  std::vector<std::string> signal_names;
  for (Json::Value::iterator it=fit_params["signals"].begin();
       it!=fit_params["signals"].end(); ++it) {
    signal_names.push_back((*it).asString());
  }

  // loop over all possible signals
  const Json::Value all_signals = root["signals"];
  for (Json::Value::const_iterator it=all_signals.begin();
       it!=all_signals.end(); ++it) {
    if (std::find(signal_names.begin(), signal_names.end(),
                  it.key().asString()) == signal_names.end()) {
      continue;
    }

    // if it is one we want to use in our fit, add it
    const Json::Value signal_params = all_signals[it.key().asString()];
    std::string name = it.key().asString(); 
    std::string title = signal_params.get("title", name).asString();
    float sigma = signal_params.get("constraint", 0.0).asFloat()
      * this->live_time * this->efficiency;
    float nexpected = signal_params["rate"].asFloat()
      * this->live_time * this->efficiency;

    // load data
    std::cout << "FitConfig::FitConfig: Loading data for "
              << name << std::endl;
    std::vector<std::string> filenames;
    for (Json::Value::const_iterator jt=signal_params["files"].begin();
         jt!=signal_params["files"].end(); ++jt) {
      filenames.push_back((*jt).asString());
    }

    Signal s(name,title,nexpected,sigma,this->hdf5_fields,
        this->sample_fields,this->observables,this->cuts,this->systematics,
        filenames);

    float years = 1.0 * s.nevents / (s.nexpected / this->live_time /
                  this->efficiency);

    std::cout << "FitConfig::FitConfig: Initialized PDF for " << s.name
              << " using " << s.nevents << " events (" << years << " y)"
              << std::endl;


    this->signals.push_back(s);
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

