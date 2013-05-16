/**
 * \class FitConfig
 * \brief Manages the configuration of a fit
 *
 * A utility class for reading in configuration parameters from a JSON file.
 */

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>
#include <vector>
#include "utils.h"
#include "signals.h"

class TH1D;
class TH2F;

class FitConfig {
  public:
    /**
     * Constructor
     *
     * \param filename Name of the JSON file to load
     */
    FitConfig(std::string filename);

    virtual ~FitConfig() {}

    /** Pretty-print the fit parameters */
    void print() const;

    /**
     * Project a TH2F down to an energy-only TH1F, optionally cutting on
     * radius.
     *
     * \param[in] hist2d the TH2F object to project
     * \param[in] r_range a Range defining the radial range to project
     * \return a pointer to the resulting TH1D
     */
    static TH1D* project1d(TH2F* const hist2d,
                           Range<float>* const r_range=NULL);

    unsigned mode;  //!< fit type (energy, energy+radius, ...)
    unsigned experiments;  //!< number of ensemble experiments
    unsigned steps;  //!< number of mcmc steps
    unsigned rebin_e;  //!< rebinning factor for energy
    unsigned rebin_r;  //!< rebinning factor for radius
    float confidence;  //!< confidence level for results (e.g. 0.9)
    float live_time;  //!< experiment live time in years
    float efficiency;  //!< overall efficiency correction
    float burnin_fraction;  //!< fraction of steps to use for burn-in period
    std::vector<Signal> signals;  //!< signal histograms and metadata
    std::string output_file;  //!< base filename for output
    std::string signal_name;  //!< name of the signal that is the signal
    Range<float> e_range;  //!< range of energies to include in fit
    Range<float> r_range;  //!< range of radii to include in fit

  protected:
    /**
     * Load an energy/radius ROOT TH2F from a file
     *
     * \param filename name of ROOT file containing the PDF
     * \param objname name of the PDF ROOT object ("pdf")
     * \return the PDF (a TH2F pointer)
     */
    static TH2F* load_histogram(std::string const filename,
                                std::string const objname);
};

#endif  // __CONFIG_H__

