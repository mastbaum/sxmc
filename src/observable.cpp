#include <cassert>
#include <iostream>
#include <string>
#include <json/value.h>

#include <sxmc/observable.h>

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


void Observable::print() const {
  std::cout << "  " << this->name << std::endl
    << "    Title: \"" << this->title << "\"" << std::endl
    << "    Field: \"" << this->field << "\""
                       << " (index " << this->field_index << ")" << std::endl
    << "    Units: \"" << this->units << "\"" << std::endl
    << "    Lower bound: " << this->lower << std::endl
    << "    Upper bound: " << this->upper << std::endl
    << "    Bins: " << this->bins << std::endl;
}

