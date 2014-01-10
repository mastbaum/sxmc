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
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TMath.h>

#include <sxmc/signals.h>
#include <sxmc/config.h>
#include <sxmc/utils.h>
#include <sxmc/generator.h>

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
  this->efficiency_corr = experiment_params.get("efficiency_corr", 1.0).asFloat();

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

    std::cout << std::endl << std::endl << "NEW SIGNAL: " << it.key().asString() << std::endl;

    // check if we want this signal
    if (std::find(signal_names.begin(), signal_names.end(),
          it.key().asString()) == signal_names.end()) {
      continue;
    }

    const Json::Value signal_params = all_signals[it.key().asString()];

    // check if we are building our pdf out of multiple contributions
    if (signal_params.isMember("pdfs")){
      // loop over contributions
      if (!signal_params.get("chain",true)){
        // its not chained, treat each part as a separate signal
        for (Json::Value::const_iterator jt=signal_params["pdfs"].begin();
            jt!=signal_params["pdfs"].end();++jt){

          const Json::Value contrib_params = signal_params["pdfs"][jt.key().asString()];
          Signal s = CreateSinglePDFSignal(contrib_params,jt.key().asString());
          this->signals.push_back(s);
        }
      }else{
        Signal s = CreateMultiPDFSignal(signal_params,it.key().asString());
        this->signals.push_back(s);
      }
    }else{
      Signal s = CreateSinglePDFSignal(signal_params, it.key().asString());
      this->signals.push_back(s);
    }
  }
}

Signal FitConfig::CreateSinglePDFSignal(const Json::Value signal_params, std::string name)
{
  std::string title = signal_params.get("title", name).asString();
  std::string category = signal_params.get("category","").asString();
  float sigma = signal_params.get("constraint", 0.0).asFloat() * this->live_time * this->efficiency_corr;
  float nexpected = signal_params["rate"].asFloat() * this->live_time * this->efficiency_corr;
  
  std::cout << "FitConfig::CreateSinglePDFSignal: Loading data for "
    << name << std::endl;

  std::vector<std::string> filenames;
  for (Json::Value::const_iterator it=signal_params["files"].begin();
      it!=signal_params["files"].end();++it){
      filenames.push_back((*it).asString());
  }
  Signal s(name,title,nexpected,sigma,category,
      this->hdf5_fields,this->sample_fields,
      this->observables,this->cuts,this->systematics,
      filenames);

  float years = 1.0 * s.nevents / (nexpected / this->live_time /
      this->efficiency_corr);
  std::cout << "FitConfig::CreateSinglePDFSignal: using "
    << s.nevents << " events (" << years << " y) (efficiency: "
    << s.efficiency << ")" << std::endl;

  return s;
}

Signal FitConfig::CreateMultiPDFSignal(const Json::Value signal_params, std::string name)
{
  // it is chained, get each contribution as a pdfz and add them
  std::vector<Signal> pdfs;
  std::vector<double> params;
  int nexpected_total = 0;
  for (Json::Value::const_iterator it=signal_params["pdfs"].begin();
      it!=signal_params["pdfs"].end();++it){
    const Json::Value contrib_params = signal_params["pdfs"][it.key().asString()];
    Signal s = CreateSinglePDFSignal(contrib_params,it.key().asString());
    pdfs.push_back(s);
    params.push_back(s.nexpected*10.0);
    nexpected_total += s.nexpected;
  }
  for (size_t i=0;i<this->systematics.size();i++)
    params.push_back(0);


  // Now sample all these pdfs to get samples with the correct relative weighting
  std::cout << "FitConfig::CreateMultiPDFSignal: Generating combined pdf for " << name << std::endl;
  std::pair<std::vector<float>, std::vector<int> > samples = make_fake_dataset(pdfs,this->systematics,this->observables,params,true);

  // Now create a new pdf from these samples
  std::string title = signal_params.get("title", name).asString();
  std::string category = signal_params.get("category","").asString();
  float sigma = signal_params.get("constraint", 0.0).asFloat() * this->live_time * this->efficiency_corr;
  float nexpected = signal_params.get("rate",nexpected_total).asFloat() * this->live_time * this->efficiency_corr;
  Signal s(name,title,nexpected,sigma,category,this->observables,this->systematics,samples.first,samples.second);

  float years = 1.0 * s.nevents / (s.nexpected / this->live_time /
      this->efficiency_corr);

  std::cout << "FitConfig::CreateMultiPDFSignal: Initialized PDF for " << s.name
    << " using " << s.nevents << " events (" << years << " y)"
    << std::endl;

  return s;
}


void FitConfig::print() const {
  std::cout << "Fit:" << std::endl
            << "  Fake experiments: " << this->experiments << std::endl
            << "  MCMC steps: " << this->steps << std::endl
            << "  Burn-in fraction: " << this->burnin_fraction << std::endl
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

