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

  void GetParLimits(std::vector<double> inputVector, double &min, double &max, double &nsteps, double &step);

  // These functions will not be used... just here for reference?

  /**
   * Calculates the Index that matches the pulled tan2 theta12.
   *
   * \param inTheta - the input pull of tan2 theta 12
   * \return index of bin for inTheta in tan2 theta 12 array.
   */
  int FindThetaIndex(double inTheta);
  
  /**
   * Calculates the Index that matches the pulled Delta m2 21
   *
   * \param inMass - the input pull of Delta m2 21
   * \return index of bin for inMass in Delta m2 21 Array.
   */
  int FindMassIndex(double inMass);
  
  /**
   * Calculates the Index that matches the Energy
   *
   * \param inE = the energy to lookup
   * \return index of bin for inE in Energy array
   */
  int FindEnergyIndex(double inE);
  
  /**
   * Calculates the global array Index that matches the 
   * tan2 theta 12, Delta m2 21 and energy indices
   *
   * \param iThetaIndex - the index of tan2 theta 12 in the tan2 theta 12 array
   * \param iMassIndex - the index of Delta m2 21 in the Delta m2 21 array.
   * \param iEIndex - the index of the energy in the energy array
   * \return global index of full array
   */
  int CalculateTotalIndex(int iThetaIndex, int iMassIndex, int iEIndex);

  double Reweight(double dTheta, double dMass, double dEnergy);
  
};

#endif  // __OSCTABLE_H__
