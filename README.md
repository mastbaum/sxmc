Signal Extraction with an MCMC {#mainpage}
==============================
`sxmc` is a GPU-accelerated unbinned maximum likelihood fitter based on a
Markov Chain Monte Carlo, intended for calculating confidence intervals or
limits.

Documentation
=============
`sxmc` includes both a User's Guide and thorough documentation of the code.
To build (requires [Doxygen](http://doxygen.org) 1.8.3+):

    $ make doc

The output is placed into the `doc/html` directory.

Building
========
`sxmc` requires the following libraries:

* [ROOT](http://root.cern.ch)
* [CUDA Runtime](https://developer.nvidia.com/cuda-downloads)

Set the environment variable `CUDART_ROOT` to the path to the CUDA runtime,
where CUDA headers can be found.

It also uses [hemi](https://github.com/harrism/hemi), which is included as a
git submodule. After cloning `sxmc`, run:

    $ git submodule init
    $ git submodule update

`sxmc` can run without any special hardware, but runs much faster with the
help of a CUDA-enabled GPU. To build with GPU support, set `CUDA_ROOT` to
point to your installation of the CUDA tools, for example:

    $ CUDA_ROOT=/usr/local/cuda make

If no GPU is available, `sxmc` will simply loop instead of running things in
parallel. You still need to have the CUDA headers installed, but no libraries
or hardware are required. To build run:

    $ make

By default, `sxmc` is built in debug mode.  For better performance (especially
when using GPU support), pass the `OPTIMIZE=1` flag to make:

    $ make OPTIMIZE=1

Basic Usage
===========
1. Create ROOT data files: The data used to build the PDFs is stored in
   TNtuples. The branch names match those used in the configuration file. Use
   one TNtuple per ROOT file; `sxmc` will use the first one it finds.

2. Configure fit: Set up the fit parameters and signal PDFs using a JSON-format
   configuration file. An example is provided in `config/`, and documentation
   is given below.

3. To run fits:
   `$ ./bin/sxmc config/your_file.json output_dir`

Configuration
=============
The fit is configured entirely through a JSON-format file. The JSON parser
in `sxmc` supports C-style comments.

An example is given in the `config` directory, and the following sections
describe the parameters for each section.

Fit
---
The fit section describes the parameters of the MCMC fit itself.

* `nexperiments` - The number of fake experiments to run, for ensemble
  testing.
* `nsteps` - Number of steps in the MCMC
* `burnin_fraction` - In order to reduce bias due to the choice of starting
  parameters, the first `2 * burnin_fracion * steps` steps are thrown out,
  (burn-in phase). After each set of `burnin_fraction * steps` steps, the
  jump distributions are recalculated based on the spread in each the
  parameter (default: 0.1).
* `output_prefix` - Prefix for output file names
* `debug_mode` - Accept every MCMC step
* `signal_name` - The name of the signal of interest, if any
* `signals` - A list of signal names to include in the fit. These should
  match up with the keys from `signals` section of the configuration.
* `observables` - Observable dimensions to use in the fit. These should
  correspond to keys in the `observables` subsection of the `pdfs` section.
* `cuts` - Cuts placed on the data. These should correspond to keys in the
  `observables` subsection of the `pdfs` section.

The general philosophy is that signals and observables (and cuts) are
defined in other sections of the configuration file, and those which will
actually be used in the fit are invoked in the `fit` section. This minimizes
rewriting of configuration files when trying out different combinations of
parameters.

PDFs
----
The PDFs sections defines the axes (observables) and shape parameters
(systematics) for building probability distributions from the data.

PDFs will be built from Monte Carlo data, loaded from ROOT files. The
input ROOT files should contain a TNtuple (or otherwise simple TTree), and
sxmc will load the first object it finds in the file. The names of branches
are used to identify fields for use in the fit.

* `observables` - A list of quantities observed in data. The observable
  has a string key used for identification, and the value is an object with
  the following parameters:
  * `title` - A title used for plotting (ROOT LaTeX)
  * `units` - Units, also used for plotting
  * `field` - The name of a branch containing the observable
  * `bins` - Number of bins for this PDF dimension
  * `min` - Minimum value for this PDF dimension
  * `max` - Maximum value for this PDF dimension
  * `logscale` - Show plots in log scale
  * `yrange` - Manually specify a y-axis range for plots
* `systematics` - A list of systematic parameters. Similar to observables,
  these have a string key and an object value:
  * `title` - A title used for plotting
  * `type` - The type of systematic. Currently supported are `shift`, `scale`,
and `resolution_scale`.
  * `observable_field` - The observable field (i.e. the branch name) affected
by the systematic transformation
  * `true_field` - The resolution shift moves the field toward or away from a
true value, with a branch name defined here (for resolution scalings only)
  * `mean` - An array of expected values for the parameter, which are
coefficients in a series expansion in the observable.
  * `sigma` - The (Gaussian) uncertainties for the parameters in the expansion
  * `fixed` - Fixed or floating (boolean true or false). Fixing a systematic
means all terms in the series or none.

Note that the elements in ``observables`` and ``systematics`` are not
necessarily used in the fit; they must be explicitly invoked in the
corresponding fields in the ``fit`` section.

Example:

    "pdfs": {
      "observables": {
        "energy": {
          "title": "Energy (MeV)",
          "units": "MeV",
          "field": "energy",  // This is the Ntuple field
          "bins": 10,
          "min": 5.0,
          "max": 15.0,
          "logscale": true,
          "yrange": [0, 500]
        }
      }
      "systematics": {
        "energy_scale": {
          "title": "Energy scale",
          "type": "scale",
          "observable_field": "energy",  // This is also the Ntuple field!
          "mean":  [0.0, 0.0],
          "sigma": [1e-3, 1e-5]
        }
      }
    }

This defines one observable, `energy`, with 10 bins from 5 to 15 MeV. A
systematic parameter `energy_scale` will apply a nonlinear scaling to the
energy; with two means defined, the parameterization is
`(a0 + a1 * energy) * energy`.

Signals
-------
The signals section defines the properties of signals that could be included
in the fit. To include a signal, put its name in the `signals` field in the
`fit` section.

Each signal has a string name and is defined by an object with the following
fields:

* `title` - A title used for plotting
* `filename` - Path to a ROOT file containing the data set
* `rate` - The expected rate of events per unit live time
* `scale` - A scale factor for the MC, e.g. if the MC represents 100 times
the expected rate per unit live time
* `constraint` - A Gaussian constraint on the rate (optional)
* `source` - The name of the associated source (see next section) (optional)
* `dataset` - The index of the dataset, corresponding to data defined in the
`data` section (see notes below). If unsure, set to `0`.
* `systematics` - A list of systematics to apply to this signal. The names in
this list are the keys in the `systematics` section.

The number of events per unit live time may be specified either as a rate or
using a scale factor. The rate applies to the entire data set used to build
the PDF, before cuts due to PDF extents (that efficiency is calculated by
sxmc, with any systematic parameters fixed to their mean values). The scale
option is used for the case where the MC events represent some multiple of
the expected rate, for example solar neutrino events generated with 500 times
the model flux. sxmc will abort if both ``rate`` and ``scale`` are set.

For example, here we define a signal `signal1` with some systematic
`observable1_scale`, normalized by scaling the Monte Carlo by a factor of
1/500, and with a 25% Gaussian constraint:

    "signals": {
      "signal1": {
        "title": "Signal 1",
        "filename": "/data/signal1.root",
        "systematics": ["observable1_scale"],
        "dataset": 0,
        "scale": 500.0,
        "constraint": 0.25
      }
    }

Sources
-------
The rates (normalizations) of signals can be linked together by giving them a
common "source." In this case, the source rate (scaling relative to 1) is what
is floated in the fit, and the scale factors for each signal define the signal
rates. The `sources` section contains a dictionary mapping string
source names to source definitions with the following fields:

* `title` - A title used for output and plotting
* `mean` - The mean value to use in the fit (default: 1.0)
* `sigma` - A gaussian constraint (default: 0.0, meaning no constraint)
* `fixed` - Fix the parameter in the fit.

For example, here we assert that two signals are related by a common flux:

    "sources": {
      "flux": {
        "title": "Some Flux",
        "sigma" 0.1
      }
    },
    "signals": {
      "cc": {
        "title": "Some Signal, CC interactions",
        "filename": "/data/cc.root",
        "dataset": 0,
        "scale": 500.0,
        "source": "flux"
      },
      "nc": {
        "title": "Some Signal, NC interactions",
        "filename": "data/nc.root",
        "dataset": 0,
        "scale": 1000.0,
        "source": "flux"
      }
    }

Data
----
If no data is specified, fake data will be sampled from the PDFs (supported
only for three or fewer dimensions). You can specify a specific dataset to fit
using the `data` section. The data should be ROOT files with the same Ntuple
format as the PDF data.

You can fit multiple data sets (for example, from different configurations of
the same experiment) in the same fit, by specifying dataset-specific PDFs and
associating them with the right data via the `dataset` ID tag.

There is yet one more dimension to the data: when running an ensemble of fits
(i.e. `nexperiments` in the `fit` section is greater than 1) with
explicitly-defined datasets, you must specify data for each of the experiments.

The `data` section consists of a set of key-value pairs, where the keys are
a sequential integer dataset IDs expressed as strings ("0", "1", etc.) and
the values are lists of dataset definitions. The list index corresponds to
ensemble experiment 0, 1, ..., and the dataset definition has the following
fields:

* `title`: A string title
* `filename`: Path to a ROOT file containing the data

For example, this defines a two-dataset fit with two fake experiments:

    "data": {
      "0": [
        {
          "title": "Run 1, Fake dataset 1",
          "filename": "/data/run1_fake1.root"
        },
        {
          "title": "Run 1, Fake dataset 2",
          "filename": "/data/run1_fake2.root"
        }
      ],
      "1": [
        {
          "title": "Run 2, Fake dataset 1",
          "filename": "/data/run2_fake1.root"
        },
        {
          "title": "Run 2, Fake dataset 2",
          "filename": "/data/run2_fake2.root"
        }
      ]
    }

Tests and Benchmarks
====================
``sxmc`` includes a suite of tests, focused on the GPU-based PDF evaluation
code at the core of the likelihood calculation. To run the tests:

    $ make test
    $ ./bin/test_sxmc

This will also build a benchmark utility, which determines how many events
per second the PDF code can put in a histogram. This is useful for estimating
the performance of ``sxmc`` on different GPU hardware. To run the benchmarking:

    $ make test
    $ ./bin/bench_sxmc pdfz

Examples of output on various processors:

Device                   | samples/s   | Notes 
------------------------ | :---------: | :----------------:
Intel Core i7 920        | 1.84914e+07 | CPU mode, 2.67 GHz
Nvidia GeForce GT 650M   | 5.71602e+08 | MacBookPro10,1
Nvidia GeForce GTX 580   | 1.60766e+09 | --
Nvidia Tesla K40         | 2.99546e+09 | --

Authors
=======
`sxmc` was created by Andy Mastbaum. On-GPU histogramming was originally
developed by S. Seibert.

See `LICENSE.txt` for license information.

