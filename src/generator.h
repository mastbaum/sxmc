#ifndef __GENERATOR_H__
#define __GENERATOR_H__

/**
 * \file generator.h
 *
 * Fake data generation.
 */

#include <utility>
#include <vector>
#include <TH1.h>

#include <sxmc/signal.h>
#include <sxmc/observable.h>
#include <sxmc/systematic.h>

/**
 * Make a fake data set.
 *
 * Create a fake data set by extracting ROOT histograms from each signal's PDF
 * and sampling. Only works for signals with 1, 2, or 3 observables.
 * Systematics are applied to PDFs before sampling.
 *
 * The output array of sampled observations is in a row-major format like:
 *
 *         obs.0 obs.1 obs.2
 *     ev0   0     1     2
 *     ev1   3     4     5
 *
 * \param signals List of Signals defining the PDFs to sample
 * \param systematics List of Systematics to be applied to PDFs
 * \param observables List of Observables common to all PDFs
 * \param poisson If true, Poisson-distribute the signal rates
 * \return Array with samples
 */
std::vector<float>
make_fake_dataset(std::vector<Signal>& signals,
                  std::vector<Systematic>& systematics,
                  std::vector<Observable>& observables,
                  bool poisson=true, std::string signal="");

#endif  // __GENERATOR_H__

