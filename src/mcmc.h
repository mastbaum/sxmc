#ifndef __MCMC_H__
#define __MCMC_H__

#include <vector>
#include <TNtuple.h>
#include "signals.h"

/** Sample the space, returning an ntuple of normalizations and likelihoods. */
TNtuple* mcmc(std::vector<Signal> signals, TNtuple* data, unsigned nsteps,
              float burnin_fraction);

#endif  // __MCMC_H__

