#ifndef __TEST_PDFZ_FIXTURES__
#define __TEST_PDFZ_FIXTURES__

#include <gtest/gtest.h>
#include <sxmc/pdfz.h>

class EvalHistConstructor : public ::testing::Test {
protected:
  virtual void SetUp() {
    nobservables = 1;
    nfields = 1;
    samples.resize(7);
    samples[0] = 0.1;
    samples[1] = 0.2;
    samples[2] = 0.3;
    samples[3] = 0.4;
    samples[4] = 0.5;
    samples[5] = 1.1;
    samples[6] = -0.1;

    weights.resize(7, 1);

    lower.resize(1);
    lower[0] = 0.0;
    upper.resize(1);
    upper[0] = 1.0;

    nbins.resize(1);
    nbins[0] = 2;
  }

  int nobservables;
  int nfields;
  std::vector<float> samples;
  std::vector<int> weights;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<int> nbins;
};

class EvalHistMethods : public EvalHistConstructor {
protected:
  virtual void SetUp() {
    EvalHistConstructor::SetUp();
    evaluator = \
      new pdfz::EvalHist(samples, weights, nfields, nobservables,
                         lower, upper, nbins);
    eval_points.resize(6);
    eval_points[0] = -0.1;
    eval_points[1] = 0.0;
    eval_points[2] = 0.25;
    eval_points[3] = 0.5;
    eval_points[4] = 0.75;
    eval_points[5] = 1.0;

    pdf_values = new hemi::Array<float>(20, true);
    norm = new hemi::Array<unsigned int>(3, true);
    params = new hemi::Array<double>(5, true);
    params->writeOnlyHostPtr(); // Force memory allocation
  }

  virtual void TearDown() {
    delete evaluator;
    delete pdf_values;
    delete norm;
    delete params;
    EvalHistConstructor::TearDown();
  }

  pdfz::EvalHist* evaluator;
  std::vector<float> eval_points;
  hemi::Array<float>* pdf_values;
  hemi::Array<unsigned int>* norm;
  hemi::Array<double>* params;
};

#endif // __TEST_PDFZ_FIXTURES__

