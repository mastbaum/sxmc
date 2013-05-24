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

FitConfig::FitConfig(std::string filename) {
  Json::Value root;
  Json::Reader reader;

  std::ifstream t(filename.c_str());
  std::string data((std::istreambuf_iterator<char>(t)),
                   std::istreambuf_iterator<char>());

  bool parse_ok = reader.parse(data, root);
  if (!parse_ok) {
    std::cout  << "FitConfig: JSON parse error:" << std::endl
               << reader.getFormattedErrorMessages();
    assert(false);
  }

  // experiment parameters
  const Json::Value experiment_params = root["experiment"];
  assert(experiment_params.isMember("live_time"));
  this->live_time = experiment_params["live_time"].asFloat();
  this->confidence = experiment_params["confidence"].asFloat();
  this->efficiency = experiment_params.get("efficiency", 1.0).asFloat();

  // fit parameters
  const Json::Value fit_params = root["fit"];
  this->mode = fit_params["mode"].asInt();
  this->signal_name = fit_params["signal_name"].asString();
  this->e_range.min = fit_params["energy_range"][0].asFloat();
  this->e_range.max = fit_params["energy_range"][1].asFloat();
  this->r_range.min = fit_params["radius_range"][0].asFloat();
  this->r_range.max = fit_params["radius_range"][1].asFloat();
  this->experiments = fit_params["experiments"].asInt();
  this->steps = fit_params["steps"].asInt();
  this->burnin_fraction = fit_params["burnin_fraction"].asFloat();
  this->output_file = fit_params.get("output_file", "fit_spectrum").asString();
  this->rebin_e = fit_params.get("rebin_e", 1).asInt();
  this->rebin_r = fit_params.get("rebin_r", 1).asInt();

  std::vector<std::string> fit_signal_names;
  for (Json::Value::iterator it=fit_params["signals"].begin(); it!=fit_params["signals"].end(); ++it) {
    fit_signal_names.push_back((*it).asString());
  }

  // signal parameters
  const Json::Value signal_names = root["signals"];
  for (Json::Value::const_iterator it=signal_names.begin(); it!=signal_names.end(); ++it) {
    if (std::find(fit_signal_names.begin(), fit_signal_names.end(), it.key().asString()) == fit_signal_names.end()) {
      continue;
    }

    const Json::Value signal_params = root["signals"][it.key().asString()];

    Signal s;
    s.name = it.key().asString();
    s.title = signal_params.get("title", s.name).asString();
    s.constraint = signal_params.get("constraint", 0.0).asFloat();
    s.nexpected = signal_params["rate"].asFloat() * this->live_time * this->efficiency;
    std::string filename = signal_params["filename"].asString();

    TH2F* h2d = load_histogram(filename, "pdf");
    h2d->Sumw2();
    h2d->Scale(1.0/h2d->Integral());

    if (this->mode == 0) {
      s.histogram = dynamic_cast<TH1*>(project1d(h2d, &this->r_range));
      int x1 = s.histogram->GetXaxis()->FindBin(e_range.min);
      int x2 = s.histogram->GetXaxis()->FindBin(e_range.max);
      s.nexpected *= (s.histogram->Integral(x1, x2) / h2d->Integral());
      s.histogram->Rebin(rebin_e);
    }
    else if (this->mode == 1) {
      s.histogram = dynamic_cast<TH1*>(h2d);
      int x1 = h2d->GetXaxis()->FindBin(r_range.min);
      int x2 = h2d->GetXaxis()->FindBin(r_range.max);
      int y1 = h2d->GetYaxis()->FindBin(e_range.min);
      int y2 = h2d->GetYaxis()->FindBin(e_range.max);
      double integral = h2d->Integral(x1, x2, y1, y2);
      s.nexpected *= integral / h2d->Integral();
      dynamic_cast<TH2F*>(s.histogram)->RebinY(this->rebin_e);
      dynamic_cast<TH2F*>(s.histogram)->RebinX(this->rebin_r);
    }
    else {
      std::cerr << "Unknown fit mode " << static_cast<int>(this->mode) << std::endl;
      assert(false);
    }

    s.histogram->Scale(1.0 / s.histogram->Integral());
    this->signals.push_back(s);
  }
}


TH2F* FitConfig::load_histogram(std::string const filename, std::string const objname) {
  TFile f(filename.c_str());
  assert(!f.IsZombie());
  TH2F* h = dynamic_cast<TH2F*>(f.Get(objname.c_str()));
  h->SetDirectory(0);
  return h;
}


TH1D* FitConfig::project1d(TH2F* const hist2d, Range<float>* const r_range) {
  int first_bin = 0;
  int last_bin = -1;
  if (r_range != nullptr) {
    first_bin = hist2d->GetXaxis()->FindBin(r_range->min);
    last_bin = hist2d->GetXaxis()->FindBin(r_range->max);
  }
  TH1D* h = hist2d->ProjectionY("pdf_proj", first_bin, last_bin);
  h->SetDirectory(0);
  return h;
}


void FitConfig::print() const {
  std::cout << "Fit:" << std::endl
            << "  Mode: " << static_cast<int>(this->mode) << std::endl
            << "  Fake experiments: " << this->experiments << std::endl
            << "  MCMC steps: " << this->steps << std::endl
            << "  Burn-in fraction: " << this->burnin_fraction << std::endl
            << "  Signal name: " << this->signal_name << std::endl
            << "  Energy: (" << this->e_range.min << ", " << this->e_range.max << ") MeV" << std::endl
            << "  Radius: (" << this->r_range.min << ", " << this->r_range.max << ") mm" << std::endl
            << "  Rebinning: E/" << this->rebin_e << ", R/" << this->rebin_r << std::endl
            << "  Output plot:" << this->output_file << std::endl
            << "Experiment:" << std::endl
            << "  Live time: " << this->live_time << " y" << std::endl
            << "  Confidence level: " << this->confidence << std::endl
            << "Signals:" << std::endl;

  for (std::vector<Signal>::const_iterator it=this->signals.begin(); it!=this->signals.end(); ++it) {
    std::cout << "  " << it->name << std::endl;
    std::cout << "    Title: \"" << it->title << "\"" << std::endl;
    std::cout << "    Expectation:  "<< it->nexpected << std::endl;
    std::cout << "    Constraint: ";
    if (it->constraint != 0) {
      std::cout << it->constraint << std::endl;
    }
    else {
      std::cout << "none" << std::endl;
    }
  }
}

