#ifndef __LIKELIHOOD_H__
#define __LIKELIHOOD_H__

/**
 * \file likelihood.h
 *
 * Utilities for working with likelihood functions.
*/

#include <map>
#include <string>

#include <sxmc/error_estimator.h>
#include <sxmc/interval.h>

class TNtuple;
class TH1F;

/**
 * Samples from a likelihood function.
 *
 * Wraps a TNtuple containing samples from the likelihood function, providing
 * statistics functions.
*/
class LikelihoodSpace {
  public:
    /**
     * Constructor.
     *
     * Note: The instance takes over ownership of the samples TNtuple!
     *
     * \param samples - A set of samples of the likelihood space
     * \param _cl - Confidence level for error estimation
     * \param _error_type - Type of calculation for intervals
    */
    LikelihoodSpace(TNtuple* samples, float _cl=0.683,
                    ErrorType _error_type=ERROR_CONTOUR,
                    bool _own_samples=true);

    /** Destructor. */
    virtual ~LikelihoodSpace();

    /** Get the parameters and uncertainties for the maximum L point. */
    std::map<std::string, Interval> get_best_fit() const { return ml_params; }

    /**
     * Get the parameters and uncertainties for the maximum L point,
     * along with the NLL.
    */
    std::map<std::string, Interval> get_best_fit(float& _nll) const {
      _nll = this->nll;
      return ml_params;
    }

    /** Print the parameters for the maximum-likelihood point. */
    void print_best_fit() const;

    /** Print the correlation matrix for all of the parameters. */
    void print_correlations();

    /**
     * Get a projection.
     *
     * \param name Name of dimension to project out
     * \returns ROOT TH1F with projected histogram
    */
    TH1F* get_projection(std::string name);

    /**
     * Get points within a given distance of the maximum.
     *
     * \param delta Number of likelihood units from max to include
     * \returns TNtuple with requested samples
    */
    TNtuple* get_contour(float delta);

    /**
     * Extract the best-fit parameters and uncertainties.
     *
     * \param[out] ml The likelihood at the maximum (negative for NLL)
     * \parm cl - Confidence level for the error estimate
     * \param error_type Type of uncertainty calculation to use
     * \returns A map from parameter names to Intervals
    */
    std::map<std::string, Interval>
    extract_best_fit(float& ml, float cl, ErrorType error_type=ERROR_CONTOUR);

    /**
     * Get a pointer to the TNtuple of samples.
     *
     * \returns A pointer to the samples
    */
    const TNtuple* get_samples() const { return samples; }

  private:
    std::map<std::string, Interval> ml_params;  //!< Likelihood-maximizing pars
    ErrorType error_type;  //!< Error calculation for intervals
    float nll;  //!< The maximum likelihood (negative for NLL)
    float cl;  //!< Confidence level for intervals
    TNtuple* samples;  //!< Samples of the likelihood function
    bool own_samples;  //!< We own the samples object
};

#endif  // __LIKELIHOOD_H__

