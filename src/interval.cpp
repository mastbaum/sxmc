#include <string>
#include <sstream>

#include <sxmc/interval.h>

std::string Interval::str() {
  float lower_error = this->point_estimate - this->lower;
  float upper_error = this->upper - this->point_estimate;

  std::stringstream ss;
  ss << this->point_estimate;
  if (this->one_sided) {
    ss << " <" << this->upper << " (" << 100 * this->cl << "\% CL)";
  }
  else {
    ss << " -" << lower_error << " +" << upper_error;
  }

  return ss.str();
}

