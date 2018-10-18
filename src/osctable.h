#ifndef __OSCTABLE_H__
#define __OSCTABLE_H__

/**
 * \file osctable.h
 *
 * Handlers for oscillation tables
 */


class OscTable {
 public:
  /** Constructor **/
  OscTable();
  
  /** Load from root file **/
  OscTable(std::string _filename, std::string _nutype,
	   std::vector<double>& lut_pee,
	   std::vector<double>& lut_pars);
  
  /** Destructor. **/
  virtual ~OscTable() {}

  std::string filename;
  std::string nutype;

  /**
   * Calculates the minimum, maximum of an input vector and the number of steps and step size.
   *
   * \param inputVector
   * \return min max nsteps and step passed by reference
   */
  void GetParLimits(std::vector<double> inputVector, double &min, double &max, double &nsteps, double &step);
  
};

#endif  // __OSCTABLE_H__
