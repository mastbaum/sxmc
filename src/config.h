#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>
#include <vector>
#include "utils.h"
#include "signals.h"
#include "pdfz.h"

class TH1D;
class TH2F;

/**
 * Get the index of an object in a vector.
 *
 * If the object isn't found, add it to the end and then return the index.
 * Useful for creating unique ordered lists.
 *
 * \param v The vector
 * \param o The object to locate
 * \return The index of the object
 */
template <class T>
size_t get_index_with_append(std::vector<T>& v, T o);


/**
 * \class FitConfig
 * \brief Manages the configuration of a fit
 *
 * A utility class for reading in configuration parameters from a JSON file
 * and storing the parameters of the fit.
 */
class FitConfig {
  public:
    /**
     * Constructor
     *
     * \param filename Name of the JSON file to load the configuration from
     */
    FitConfig(std::string filename);

    virtual ~FitConfig() {}

    /** Pretty-print the fit parameters */
    void print() const;

    unsigned experiments;  //!< number of ensemble experiments
    unsigned steps;  //!< number of mcmc steps
    float confidence;  //!< confidence level for results (e.g. 0.9)
    float live_time;  //!< experiment live time in years
    float efficiency;  //!< overall efficiency correction
    float burnin_fraction;  //!< fraction of steps to use for burn-in period
    float signal_eff;  //!< signal efficiency, i.e. fraction of signal in fit
    bool debug_mode;  //!< enable/disable debugging mode (accept/save all)
    std::string output_file;  //!< base filename for output
    std::string signal_name;  //!< name of the signal that is the signal
    std::vector<Signal> signals;  //!< signal histograms and metadata
    std::vector<Systematic> systematics;  //!< Systematics used in PDFs
    std::vector<Observable> observables;  //!< Observables used in PDFs
    std::vector<Observable> cuts;  //!< Cuts applied before fit
};

#endif  // __CONFIG_H__

