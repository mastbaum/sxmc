#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <json/value.h>

#include <sxmc/systematic.h>
#include <sxmc/pdfz.h>

Systematic::Systematic(const std::string _name, const Json::Value& config)
      : name(_name), means(NULL), sigmas(NULL) {
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
  else if (type_string == "oscillation") {
    this->type = pdfz::Systematic::OSCILLATION;
    assert(config.isMember("osc_lut"));
    this->osc_lut = config["osc_lut"].asString();
    assert(config.isMember("truth_field"));
    this->truth_field = config["truth_field"].asString();
    assert(config.isMember("pid_field"));
    this->pid_field = config["pid_field"].asString();
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


void Systematic::print() const {
  std::cout << "  " << this->name << std::endl
            << "    Title: \"" << this->title << "\"" << std::endl
            << "    Type: " << this->type << std::endl
            << "    Observable: " << this->observable_field
            << " (index " << this->observable_field_index << ")" << std::endl;
  if (this->type == pdfz::Systematic::RESOLUTION_SCALE ||
      this->type == pdfz::Systematic::OSCILLATION) {
    std::cout << "    Truth: " << this->truth_field
            << " (index " << this->truth_field_index << ")" << std::endl;
  }
  if (this->type == pdfz::Systematic::OSCILLATION) {
    std::cout << "    PID: " << this->pid_field
            << " (index " << this->pid_field_index << ")" << std::endl;
  }
  std::cout << "    Means: ";
  for (size_t i=0; i<this->npars; i++) {
    std::cout << "(" << i << ") " << this->means[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "    Constraints: ";
  for (size_t i=0; i<this->npars; i++) {
    std::cout << "(" << i << ") " << this->sigmas[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "    Fixed: " << (this->fixed ? "yes" : "no") << std::endl;
}

