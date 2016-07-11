#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>
#include <json/value.h>
#include <json/reader.h>

#include <sxmc/config.h>
#include <sxmc/ttree_io.h>
#include <sxmc/signal.h>
#include <sxmc/utils.h>

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

  const Json::Value fit_params = root["fit"];
  const Json::Value obs_params = root["pdfs"]["observables"];
  const Json::Value sys_params = root["pdfs"]["systematics"];
  const Json::Value rate_params = root["rates"];
  const Json::Value sig_params = root["signals"];
  const Json::Value src_params = root["sources"];

  // Load general fit parameters
  assert(fit_params.isMember("nexperiments"));
  this->nexperiments = fit_params["nexperiments"].asInt();
  assert(this->nexperiments > 0);

  assert(fit_params.isMember("nsteps"));
  this->nsteps = fit_params["nsteps"].asInt();
  assert(this->nsteps > 0);

  this->samples = fit_params.get("samples", "").asString();
  std::string err_type = fit_params.get("error_type", "contour").asString();
  if (err_type == "projection") {
    this->error_type = ERROR_PROJECTION;
  }
  else if (err_type == "contour") {
    this->error_type = ERROR_CONTOUR;
  }
  else {
    std::cerr << "FitConfig: Unknown error type \"" << err_type << "\""
              << std::endl;
    assert(false);
  }

  this->burnin_fraction = fit_params.get("burnin_fraction", 0.1).asFloat();
  this->debug_mode = fit_params.get("debug_mode", false).asBool();
  this->output_prefix = fit_params.get("output_prefix", "lspace").asString();
  this->plots = fit_params.get("plots", true).asBool();
  this->sensitivity = fit_params.get("sensitivity", false).asBool();
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

  // Systematics and sources: union of all used in signals
  short sidx = 0;
  short pidx = 0;
  for (Json::Value::const_iterator it=sig_params.begin();
       it!=sig_params.end(); ++it) {
    std::string signal_name = it.key().asString();

    // Systematics
    if (sig_params[signal_name].isMember("systematics")) {
      const Json::Value slist = sig_params[signal_name]["systematics"];
      for (Json::Value::const_iterator it=slist.begin();
           it!=slist.end(); ++it) {
        std::string sys_name = (*it).asString();
        assert(sys_params.isMember(sys_name));
        bool exists = false;
        for (size_t k=0; k<this->systematics.size(); k++) {
          if (this->systematics[k].name == sys_name) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          this->systematics.push_back(
            Systematic(sys_name, sys_params[sys_name]));
          Systematic* s = &this->systematics.back();
          for (size_t k=0; k<s->npars; k++) {
            s->pidx.push_back(pidx++);
          }
        }
      }
    }

    // Sources
    if (sig_params[signal_name].isMember("source")) {
      std::string src_name = sig_params[signal_name]["source"].asString();
      assert(src_params.isMember(src_name));
      bool exists = false;
      for (size_t k=0; k<this->sources.size(); k++) {
        if (this->sources[k].name == src_name) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        Source s(src_name, src_params[src_name]);
        s.index = sidx++;
        this->sources.push_back(s);
      }
    }
    else {
      // The signal is a source for itself
      this->sources.push_back(Source(
        signal_name, sidx++,
        sig_params[signal_name].get("mean", 1.0).asFloat(),
        sig_params[signal_name].get("sigma", 0.0).asFloat(),
        sig_params[signal_name].get("fixed", false).asBool()
      ));
    }
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
  for (Json::Value::const_iterator it=fit_params["signals"].begin();
       it!=fit_params["signals"].end(); ++it) {
    std::string name = (*it).asString();

    assert(sig_params.isMember(name));
    const Json::Value config = sig_params[name];

    std::cout << "FitConfig: Loading signal: " << name << std::endl;

    assert(config.isMember("dataset"));
    unsigned dataset = config["dataset"].asUInt();

    this->datasets.insert(dataset);

    assert((!config.isMember("rate") &&  config.isMember("scale")) ||
           ( config.isMember("rate") && !config.isMember("scale")));

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

    // Systematics
    std::vector<Systematic> systs;
    if (config.isMember("systematics")) {
      const Json::Value slist = config["systematics"];
      for (Json::Value::const_iterator it=slist.begin();
           it!=slist.end(); ++it) {
        std::string sys_name = (*it).asString();
        for (size_t k=0; k<this->systematics.size(); k++) {
          if (this->systematics[k].name == sys_name) {
            systs.push_back(this->systematics[k]);
          }
        }
      }
    }

    // Source
    Source source;
    std::string source_name = config.get("source", name).asString();
    for (size_t k=0; k<this->sources.size(); k++) {
      if (this->sources[k].name == source_name) {
        source = this->sources[k];
        break;
      }
    }

    this->signals.push_back(
      Signal(name, title, filename, dataset, source, nexpected,
             this->sample_fields, this->observables,
             this->cuts, systs));
  }

  //// Load data
  if (root.isMember("data")) {
    const Json::Value data_params = root["data"];

    for (Json::Value::const_iterator it=data_params.begin();
         it!=data_params.end(); ++it) {
      std::string ds_name = it.key().asString();
      unsigned dataset;
      std::istringstream is(ds_name);
      is >> dataset;

      const Json::Value dconfig = data_params[ds_name];

      for (Json::Value::const_iterator jt=dconfig.begin();
           jt!=dconfig.end(); ++jt) {

        const Json::Value drow = *jt;
        std::string title = drow["title"].asString();
        std::string filename = drow["filename"].asString();

        // All active observables are treated as cuts, in order to clip the data
        // to the PDF boundaries.
        std::vector<Observable> cc = this->observables;
        for (size_t ii=0; ii<this->cuts.size(); ii++) {
          cc.push_back(this->cuts[ii]);
        }

        std::vector<Systematic> syst;  // No systematics applied to data!

        std::cout << "FitConfig: Loading data: " << filename << std::endl;
        this->data[dataset].push_back(
          Signal(title, title, filename, dataset, Source(),
                 -1.0, this->sample_fields,
                 this->observables, cc, syst));
      }
    }
  }
}


void FitConfig::print() const {
  std::cout << "Fit:" << std::endl
    << "  Number of experiments: " << this->nexperiments << std::endl
    << "  MCMC steps: " << this->nsteps << std::endl
    << "  Burn-in fraction: " << this->burnin_fraction << std::endl
    << "  Random seed (0=random): " << this->seed << std::endl
    << "  Confidence level: " << this->confidence << std::endl;

  if (this->samples != "") {
    std::cout << "  Samples file: " << this->samples << std::endl;
  }

  std::cout << "Signals:" << std::endl;
  for (std::vector<Signal>::const_iterator it=this->signals.begin();
      it!=this->signals.end(); ++it) {
    it->print();
  }

  std::cout << "Sources:" << std::endl;
  for (std::vector<Source>::const_iterator it=this->sources.begin();
      it!=this->sources.end(); ++it) {
    it->print();
  }

  std::cout << "Cuts:" << std::endl;
  for (std::vector<Observable>::const_iterator it=this->cuts.begin();
      it!=this->cuts.end(); ++it) {
    it->print();
  }

  std::cout << "Observables:" << std::endl;
  for (std::vector<Observable>::const_iterator it=this->observables.begin();
      it!=this->observables.end(); ++it) {
    it->print();
  }

  if (this->systematics.size() > 0) {
    std::cout << "Systematics:" << std::endl;
    for (std::vector<Systematic>::const_iterator it=this->systematics.begin();
         it!=this->systematics.end(); ++it) {
      it->print();
    }
  }
}

