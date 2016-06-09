#include <iostream>
#include <string>
#include <json/value.h>

#include <sxmc/source.h>

Source::Source(const std::string& _name, const Json::Value& params)
    : name(_name) {
  this->mean = params.get("mean", 1.0).asFloat();
  this->sigma = params.get("sigma", 0.0).asFloat();
  this->fixed = params.get("fixed", false).asBool();
}


void Source::print() const {
  std::cout << "  " << this->name << std::endl
    << "    Mean: " << this->mean << std::endl
    << "    Constraint: ";
  if (this->sigma != 0) {
    std::cout << this->sigma << std::endl;
  }
  else {
    std::cout << "none" << std::endl;
  }
  std::cout << "    Fixed: " << (this->fixed ? "yes" : "no") << std::endl;
}

