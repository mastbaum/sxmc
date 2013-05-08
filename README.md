Signal Extraction with MCMC
===========================
A GPU-accelerated fit program which calculates fully frequentist confidence
intervals or limits, computing the likelihood space with a Markov Chain
Monte Carlo.

Building
--------
`sxmc` requires the following libraries:

* CUDA (and the nvcc compiler)
* JsonCpp
* ROOT
* RAT (to make PDFs from MC files)

It also requires a CUDA-enabled Nvidia GPU.

To build, run `make` and specify a `CUDA_ROOT` environment variable. E.g.,

    $ make CUDA_ROOT=/opt/cuda-5.0

Usage
-----
1. Create PDFs: PDFs are ROOT TH2Fs with event energy and radius dimensions.
   Each PDF is stored in its own ROOT file and named "pdf".

2. Configure fit: Set up the fit parameters and signal PDFs using a JSON-format
   configuration file. An example is provided in `config/`.

3. To calculate signal sensitivity, run: `$ ./bin/sensitivity config/your_file.json`

Configuration Files
-------------------
The fit is controlled via a JSON-format configuration file. See
`config/example.json` for an annotated example.

