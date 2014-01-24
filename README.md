Signal Extraction with MCMC
===========================
A GPU-accelerated unbinned maximum likelihood fit using a Markov Chain Monte
Carlo, intended for calculating confidence intervals or limits.

Building
--------
`sxmc` requires the following libraries:

* [ROOT](http://root.cern.ch)
* [Doxygen](http://doxygen.org) (if building documentation)

It also uses [hemi](https://github.com/harrism/hemi), which is included as a
git submodule. After cloning `sxmc`, run:

    $ git submodule init
    $ git submodule update

`sxmc` runs much faster with the help of a CUDA-enabled GPU. To build with GPU
support, set `CUDA_ROOT` to point to your installation of the CUDA tools:

    $ CUDA_ROOT=/usr/local/cuda make

If no GPU is available, `sxmc` will simply loop instead of running things in
parallel. To build without GPU support:

    $ make

You still need to have the CUDA headers installed, but no libraries or hardware
are required.

By default, `sxmc` is built in debug mode.  For much higher performance when
using GPU support, pass the `OPTIMIZE=1` flag to make:

    $ make OPTIMIZE=1

Documentation
-------------
The code is fully documented for Doxygen. To view HTML documentation online,
see [here](http://mastbaum.github.io/sxmc/doc/html). To build HTML and LaTeX
documentation yourself, run

    $ make doc

The output is placed into the `doc` directory.

Usage
-----
1. Create ROOT data files: The data used to build the PDFs is stored in
   TNtuples. The branch names match those used in the configuration file.

2. Configure fit: Set up the fit parameters and signal PDFs using a JSON-format
   configuration file. An example is provided in `config/`.

3. To run fits:
   `$ ./bin/sxmc config/your_file.json output_dir`

