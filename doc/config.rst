Fit Configuration
=================
The fit is configured entirely through a JSON-format file. The JSON parser
in sxmc supports C-style comments.

The format is the following, with explanations provided below::

    {
      "fit": {
        "experiments": 1,
        "steps": 300000,
        "burnin_fraction": 0.2,
        "output_file": "fit_example",
        "debug_mode": false,
        "signals": [
          "hepcc", "b8cc", "atm_nu_no_osc", "dsnbcc"
        ],
        "observables": ["energy"],
        "cuts": ["radius"],
        "systematics": []
      },
    
      "experiment": {
        "live_time": 1.0,
        "confidence": 0.9
      },
    
      "pdfs": {
        "observables": {
          "energy": {
            "title": "Teff (MeV)",
            "units": "MeV",
            "field": "Enekp",
            "bins": 48,
            "min": 6.0,
            "max": 30.0
          },
          "radius": {
            "title": "R (cm)",
            "units": "cm",
            "field": "Rfp",
            "bins": 10,
            "min": 0.0,
            "max": 550.0
          }
        }
      },
      "signals": {
        "hepcc": {
          "title": "hep CC",
          "files": ["pdfs/hepcc.root"],
          "rate": 34.2
        },
        "b8cc": {
          "title": "B8 CC",
          "files": ["pdfs/b8cc.root"],
          "rate": 5743.7
        },
      }
    }

Fit
---
The fit section describes the parameters of the MCMC fit itself.

* ``experiments`` - The number of fake experiments to run, for ensemble
  testing.
* ``steps`` - Number of steps in the MCMC
* ``burnin_fraction`` - In order to reduce bias due to the choice of starting
  parameters, the first `2 * burnin_fracion * steps` steps are thrown out,
  in what is known as the burn-in phase. After each set of
  `burnin_fraction * steps` steps, the jump distributions are recalculated
  based on the spread in each the parameter.
* ``output_file`` - Prefix for output file names
* ``debug_mode`` - Accept every MCMC step
* ``signals`` - A list of signal names to include in the fit. These should
  match up with the keys from ``signals`` section of the configuration.
* ``observables`` - Observable dimensions to use in the fit. These should
  correspond to keys in the ``observables`` subsection of the ``pdfs`` section.
* ``cuts`` - Cuts placed on the data. These should correspond to keys in the
  ``observables`` subsection of the ``pdfs`` section.
* ``systematics`` - Systematics to float in the fit. These should correspond
  to keys in the ``systematics`` subsection of the ``pdfs`` section. Note that
  floating systematics will make the fit much slower.

The general philosophy is that signals, observables, and systematics are
defined in other sections of the configuration file, and those which will
actually be used in the fit are invoked in the ``fit`` section. This
minimizes rewriting of configuration files when trying out different
combinations of parameters.

Experiment
----------
The experiment section describes the parameters of a single measurement,
real or fake.

* ``live_time`` - A scale factor applied to the expected rate
* ``efficiency_correction`` - A scale factor applied to the expected rate
* ``confidence`` - A confidence level for a limit. Currently unused.
* ``data`` - The name of a signal containing data to fit, or empty to randomly
  sample from the PDFs. Sampling only works in >= 3 observable dimensions.

PDFs
----
The PDFs sections defines the axes (observables) and shape parameters
(systematics) for building probability distributions from the data.

PDFs will be built from Monte Carlo data, loaded from ROOT files. The
input ROOT files should contain a TNtuple (or otherwise simple TTree), and
sxmc will load the first object it finds in the file. The names of branches
are used to identify fields for use in the fit.

* ``observables`` - A list of quantities observed in data. The observable
  has a string key used for identification, and the value is an object with
  the following parameters:

  * ``title`` - A title used for plotting (ROOT LaTeX)
  * ``units`` - Units, also used for plotting
  * ``field`` - The name of a branch containing the observable
  * ``bins`` - Number of bins for this PDF dimension
  * ``min`` - Minimum value for this PDF dimension
  * ``max`` - Maximum value for this PDF dimension
  * ``logscale`` - Show plots in log scale
  * ``yrange`` - Manually specify a y-axis range for plots
* ``systematics`` - A list of systematic parameters. Similar to observables,
  these have a string key and an object value:

  * ``title`` - A title used for plotting
  * ``type`` - The type of systematic. Currently supported are "shift", "scale", and "resolution_scale"
  * ``observable_field`` - The observable field (i.e. the branch name) affected by the systematic transformation
  * ``true_field`` - The resolution shift moves the field toward or away from a true value, with a branch name defined here
  * ``mean`` - The expected value for the parameter
  * ``sigma`` - The (Gaussian) uncertainty in the parameter
  * ``fixed`` - Fixed or floating (boolean true or false)

Note that the elements in ``observables`` and ``systematics`` are not
necessarily used in the fit; they must be explicitly invoked in the
corresponding fields in the ``fit`` section.

Signals
-------
The signals section defines the properties of signals that could be included in the fit.

Each signal has a string name and is defined by an object with the following
fields:

* ``title`` - A title used for plotting
* ``files`` - A list of ROOT files containing the data set
* ``rate`` - The expected rate of events per unit live time
* ``scale`` - A scale factor for the MC, e.g. if the MC represents 100 times
  the expected rate per unit live time
* ``constraint`` `(optional)` - A Gaussian constraint on the rate

The number of events per unit live time may be specified either as a rate or
using a scale factor. The rate applies to the entire data set used to build
the PDF, before cuts due to PDF extents (that efficiency is calculated by
sxmc, with any systematic parameters fixed to their mean values). The scale
option is used for the case where the MC events represent some multiple of
the expected rate, for example solar neutrino events generated with 500 times
the SSM flux. sxmc will abort if both ``rate`` and ``scale`` are set.

