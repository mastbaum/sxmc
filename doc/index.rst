.. sxmc documentation master file, created by
   sphinx-quickstart on Sat Mar  5 15:08:27 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to sxmc's documentation!
================================
sxmc is a GPU-accelerated unbinned maximum likelihood fitter using a Markov
Chain Monte Carlo, intended for calculating confidence intervals or limits.

Contents:

.. toctree::
   :maxdepth: 2

   index
   config

Documentation
-------------
In addition to this User's Guide, the code is thoroughly documented with
Doxygen. To build HTML and LaTeX documentation, run::

    $ make doc

The output is placed into the `doc` directory, with HTML output
`here </Users/mastbaum/projects/sxmc/doc/html/index.html>`_.

Basic Usage
-----------
1. Create ROOT data files: The data used to build the PDFs is stored in
   TNtuples. The branch names match those used in the configuration file.

2. Configure fit: Set up the fit parameters and signal PDFs using a JSON-format
   configuration file. An example is provided in `config`.

3. To run fits::

   $ make
   $ ./bin/sxmc config/your_file.json output_dir

Building
--------
``sxmc`` requires the following libraries:

* `ROOT <http://root.cern.ch>`_
* `Doxygen <http://doxygen.org>`_ (if building documentation)
* `CUDA Runtime <https://developer.nvidia.com/cuda-downloads>`_

It also uses `hemi <https://github.com/harrism/hemi>`_, which is included as a
git submodule. After cloning ``sxmc``, run::

    $ git submodule init
    $ git submodule update

``sxmc`` can be compiled to run either entirely on a CPU, or with the PDF
building and likelihood calculation accelerated with a CUDA-enabled GPU.
The latter is much faster, particularly when PDF shape parameters may vary
in the fit. CUDA headers are always required, but no libraries or GPU hardware
are necessary for building in CPU mode.

Set the path to the CUDA runtime with the environment variable ``CUDART_ROOT``.
To build with GPU support, set ``CUDA_ROOT`` to point to your installation of
the CUDA tools. For example::

    $ CUDA_ROOT=/usr/local/cuda make

If no GPU is available, build without GPU support::

    $ make

By default, ``sxmc`` is built in debug mode.  For much higher performance when
using GPU support, pass the ``OPTIMIZE=1`` flag to make::

    $ make OPTIMIZE=1

Unit Tests and Benchmarks
-------------------------
``sxmc`` includes a suite of tests, focused on the GPU-based PDF evaluation
code at the core of the likelihood calculation. To run the tests::

    $ make test
    $ ./bin/test_sxmc

This will also build a benchmark utility, which determines how many events
per second the PDF code can put in a histogram. This is useful for estimating
the performance of ``sxmc`` on different GPU hardware. To run the benchmark::

    $ make test
    $ ./bin/bench_sxmc pdfz

Examples of output on various processors:

================= =========== ==================
Device            samples/s   Notes 
================= =========== ==================
Intel Core i7 920 1.84914e+07 CPU mode, 2.67 GHz
GeForce GT 650M   5.71602e+08 MacBookPro10,1
GeForce GTX 580   1.60766e+09
Tesla K40         2.99546e+09
================= =========== ==================



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

