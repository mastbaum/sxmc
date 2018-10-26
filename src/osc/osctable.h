#ifndef __OSCTABLE_H__
#define __OSCTABLE_H__

/**
 * \file osctable.h
 *
 * Handlers for oscillation tables
 */

namespace sxmc {
  namespace osc {
    /**
     * Opens a ROOT file and reads the given dataset of oscillation values to a float
     * vector.
     *
     * \param filename - Name of file to read from
     * \param nutype - If a file is for hep or 8B
     * \param lut_pee - a vector of probability variables (returned by reference)
     * \param lut_pars - a vector of array parameters (returned by reference)
     * \return nothing
     */
    void load_oscillation_table(std::string filename, std::string nutype,
				std::vector<double>& lut_pee,
				std::vector<double>& lut_pars);
  
    /**
     * Calculates the minimum, maximum of an input vector and the number of steps and step size.
     *
     * \param inputVector - Vector to extract values from
     * \param min - Minimum entry of vector (returned by reference)
     * \param max - Maximum entry of vector (returned by reference)
     * \param nsteps - The number of steps between minimum and maximum of vector (returned by reference)
     * \param step - The size of the step between one entry and the next (returned by reference)
     * \return nothing
     */
    void get_parameter_limits(std::vector<double> inputVector, double &min, double &max, double &nsteps, double &step);
    
  } // namespace osc
} // namespace sxmc

#endif  // __OSCTABLE_H__
