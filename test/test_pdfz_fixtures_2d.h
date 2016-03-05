#ifndef __TEST_PDFZ_FIXTURES_2D__
#define __TEST_PDFZ_FIXTURES_2D__

#include <gtest/gtest.h>
#include <sxmc/pdfz.h>

class EvalHist2DConstructor : public ::testing::Test {
protected:
  virtual void SetUp() {
    nobservables = 2;
    nfields = 2;
    samples.resize(14);
    samples[0] = 0.4; samples[1] = 10.5;
    samples[2] = 0.5; samples[3] = 11.0;
    samples[4] = 0.75; samples[5] = 11.0;
    samples[6] = 0.6; samples[7] = 11.5;
    samples[8] = 0.6; samples[9] = 11.8;
    samples[10] = 0.9; samples[11] = 11.5;
    samples[12] = 0.4; samples[13] = 12.0;

    lower.resize(2);
    lower[0] = 0.0; lower[1] = 10.0;

    upper.resize(2);
    upper[0] = 1.0; upper[1] = 12.0;

    nbins.resize(2);
    nbins[0] = 2; nbins[1] = 3;
  }

  // virtual void TearDown() {}
  int nobservables;
  int nfields;
  std::vector<float> samples;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<int> nbins;
};

class EvalHist2DMethods : public EvalHist2DConstructor {
protected:
    virtual void SetUp() {
        EvalHist2DConstructor::SetUp();
        evaluator = new pdfz::EvalHist(samples, nfields, nobservables, lower, upper, nbins);
        eval_points.resize(16);
        eval_points[0] = 0.2; eval_points[1] = 10.2;
        eval_points[2] = 0.7; eval_points[3] = 10.4;
        eval_points[4] = 0.5; eval_points[5] = 11.0;
        eval_points[6] = 0.25; eval_points[7] = 11.8;
        eval_points[8] = 0.9; eval_points[9] = 11.9;
        eval_points[10] = 0.3; eval_points[11] = 12.0;
        eval_points[12] = 0.3; eval_points[13] = 13.0;
        eval_points[14] = 0.3; eval_points[15] = 5.0;

        pdf_values = new hemi::Array<float>(40, true);
        norm = new hemi::Array<unsigned int>(3, true);
        params = new hemi::Array<double>(5, true);
        params->writeOnlyHostPtr(); // Force memory allocation
    }

    virtual void TearDown() {
        delete evaluator;
        delete pdf_values;
        delete norm;
        delete params;
        EvalHist2DConstructor::TearDown();
    }


    pdfz::EvalHist *evaluator;
    std::vector<float> eval_points;
    hemi::Array<float> *pdf_values;
    hemi::Array<unsigned int> *norm;
    hemi::Array<double> *params;
};



#endif // __TEST_PDFZ_FIXTURES_2D__
