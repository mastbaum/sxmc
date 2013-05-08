#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <cuda.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH1.h>
#include "signals.h"
#include "nll.h"

// macro to print error and abort on cuda errors, c/o stan
#define CUDA_CHECK_ERROR(call) do { \
  cudaError err = call; \
  if (cudaSuccess != err) { \
    fprintf(stderr, "Cuda error in file '%s' in line %i: %s.\n", \
            __FILE__, __LINE__, cudaGetErrorString(err)); \
    exit(EXIT_FAILURE); \
  } \
} while (0)


// compute the log likelihood given signal rates and normalizations
__global__ void ll(const float* lut, const double* pars, const size_t ne,
                   const size_t ns, double* sums) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  for (int i=idx; i<(int)ne; i+=gridDim.x*blockDim.x) {
    double s = 0;
    for (size_t j=0; j<ns; j++) {
      s += pars[j] * lut[i*ns+j];
    }
    sums[idx] += log(s);
  }
}


NLL::NLL(const std::vector<Signal>& signals, TNtuple* data) {
  this->nsignals = signals.size();
  this->nevents = data->GetEntries();

  this->expectations = new double[nsignals];
  this->constraints = new double[nsignals];

  for (size_t i=0; i<nsignals; i++) {
    this->expectations[i] = signals[i].nexpected;
    this->constraints[i] = signals[i].constraint;
  }

  this->lut = NLL::build_lut(signals, data);

  CUDA_CHECK_ERROR(cudaMalloc(&this->lut_device,
                              this->nevents * this->nsignals * sizeof(float)));

  CUDA_CHECK_ERROR(cudaMemcpy(this->lut_device, lut,
                              this->nevents * this->nsignals * sizeof(float),
                              cudaMemcpyHostToDevice));
}


NLL::~NLL() {
  delete this->lut;
  delete[] this->expectations;
  delete[] this->constraints;

  CUDA_CHECK_ERROR(cudaFree(this->lut_device));
}


float* NLL::build_lut(const std::vector<Signal>& signals, TNtuple* data) {
  int nevents = data->GetEntries();
  float* lut = new float[signals.size() * nevents];

  std::vector<float> minima;
  for (size_t i=0; i<signals.size(); i++) {
    minima.push_back(signals[i].histogram->GetMinimum(0) * 0.0001);
  }

  float e;
  float r;
  data->SetBranchAddress("e", &e);
  data->SetBranchAddress("r", &r);

  for (int i=0; i<nevents; i++) {
    data->GetEntry(i);
    for (size_t j=0; j<signals.size(); j++) {
      double v = 0;
      if (signals[j].histogram->IsA() == TH2F::Class()) {
        v = dynamic_cast<TH2F*>(signals[j].histogram)->Interpolate(r, e);
      }
      else if (signals[j].histogram->IsA() == TH1D::Class()) {
        v = dynamic_cast<TH1D*>(signals[j].histogram)->Interpolate(e);
      }
      else {
        std::cerr << "build_lut: Unknown histogram class "
                  << signals[j].histogram->ClassName() << std::endl;
        assert(false);
      }

      if (v <= 0) {
        v = minima[j];
      }
      lut[i * signals.size() + j] = v;
    }
  }

  return lut;
}


double NLL::operator()(double* norms) {
  double result = 0;

  // N + fractional constraints
  for (size_t i=0; i<this->nsignals; i++) {
    if (norms[i] < 0) {
      return 1e10;
    }
    result += norms[i];
    if (this->constraints[i] > 0) {
      result += 0.5 * pow((norms[i]/this->expectations[i] - 1) /
                          this->constraints[i], 2);
    }
  }

  // loop over events, on the gpu
  size_t blocksize = 256;
  size_t nblocks = 16;
  size_t nthreads = nblocks * blocksize;

  double* norms_device;
  CUDA_CHECK_ERROR(cudaMalloc(&norms_device, this->nsignals * sizeof(double)));

  CUDA_CHECK_ERROR(cudaMemcpy(norms_device, norms,
                              this->nsignals * sizeof(double),
                              cudaMemcpyHostToDevice));

  double* sums_device;
  CUDA_CHECK_ERROR(cudaMalloc(&sums_device, nthreads * sizeof(double)));
  CUDA_CHECK_ERROR(cudaMemset(sums_device, 0, nthreads * sizeof(double)));

  ll<<<nblocks, blocksize>>>(this->lut_device, norms_device, this->nevents,
                             this->nsignals, sums_device);

  CUDA_CHECK_ERROR(cudaThreadSynchronize());

  double* sums = new double[blocksize * nblocks];

  CUDA_CHECK_ERROR(cudaMemcpy(sums, sums_device, nthreads * sizeof(double),
                              cudaMemcpyDeviceToHost));

  CUDA_CHECK_ERROR(cudaFree(norms_device));
  CUDA_CHECK_ERROR(cudaFree(sums_device));

  double sum = 0;
  for (size_t i=0; i<nthreads; i++) {
    sum += sums[i];
  }

  delete[] sums;

  result -= sum;

  return result;
}

