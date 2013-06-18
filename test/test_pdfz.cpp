#include <gtest/gtest.h>
#include "pdfz.h"
#include <cmath>

#ifndef __CUDACC__
#define isnan std::isnan
#endif

TEST(PdfzError, Constructor)
{
    pdfz::Error err("test");
    EXPECT_EQ(err.msg, "test");
}

TEST(SystematicObjects, ShiftSystematicConstructor)
{
    pdfz::ShiftSystematic syst(1, 0);
    EXPECT_EQ(pdfz::Systematic::SHIFT, syst.type);
    EXPECT_EQ(1, syst.obs);
    EXPECT_EQ(0, syst.par);
}

TEST(SystematicObjects, ScaleSystematicConstructor)
{
    pdfz::ScaleSystematic syst(0, 3);
    EXPECT_EQ(pdfz::Systematic::SCALE, syst.type);
    EXPECT_EQ(0, syst.obs);
    EXPECT_EQ(3, syst.par);
}

TEST(SystematicObjects, ResolutionSystematicConstructor)
{
    pdfz::ResolutionSystematic syst(0, 2, 1);
    EXPECT_EQ(pdfz::Systematic::RESOLUTION, syst.type);
    EXPECT_EQ(0, syst.obs);
    EXPECT_EQ(2, syst.true_obs);
    EXPECT_EQ(1, syst.par);
}

////////////////////////////////

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

    lower.resize(1);
    lower[0] = 0.0;
    upper.resize(1);
    upper[0] = 1.0;

    nbins.resize(1);
    nbins[0] = 2;
  }

  // virtual void TearDown() {}
  int nobservables;
  int nfields;
  std::vector<float> samples;
  std::vector<float> lower;
  std::vector<float> upper;
  std::vector<int> nbins;
};

TEST_F(EvalHistConstructor, WrongSampleSize)
{
    ASSERT_THROW(pdfz::EvalHist(samples, 2 /* nfields */, nobservables, lower, upper, nbins), pdfz::Error);
}

TEST_F(EvalHistConstructor, NobsLargerThanNfields)
{
    ASSERT_THROW(pdfz::EvalHist(samples, nfields , 7 /* nobservables */, lower, upper, nbins), pdfz::Error);
}

TEST_F(EvalHistConstructor, WrongLowerSize)
{
    lower.resize(2);
    ASSERT_THROW(pdfz::EvalHist(samples, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}

TEST_F(EvalHistConstructor, WrongUpperSize)
{
    upper.resize(2);
    ASSERT_THROW(pdfz::EvalHist(samples, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}

TEST_F(EvalHistConstructor, WrongNbinsSize)
{
    nbins.resize(2);
    ASSERT_THROW(pdfz::EvalHist(samples, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}

TEST_F(EvalHistConstructor, ZeroBins)
{
    nbins[0] = 0;
    ASSERT_THROW(pdfz::EvalHist(samples, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}

///////////////

class EvalHistMethods : public EvalHistConstructor {
  protected:
    virtual void SetUp() {
        EvalHistConstructor::SetUp();
        evaluator = new pdfz::EvalHist(samples, nfields, nobservables, lower, upper, nbins);
        eval_points.resize(6);
        eval_points[0] = -0.1;
        eval_points[1] = 0.0;
        eval_points[2] = 0.25;
        eval_points[3] = 0.5;
        eval_points[4] = 0.75;
        eval_points[5] = 1.0;

        pdf_values = new hemi::Array<float>(20, true);
        norm = new hemi::Array<unsigned int>(3, true);
        params = new hemi::Array<float>(5, true);
    }

    virtual void TearDown() {
        delete evaluator;
        delete pdf_values;
        delete norm;
        delete params;
        EvalHistConstructor::TearDown();
    }


    pdfz::EvalHist *evaluator;
    std::vector<float> eval_points;
    hemi::Array<float> *pdf_values;
    hemi::Array<unsigned int> *norm;
    hemi::Array<float> *params;
};


TEST_F(EvalHistMethods, Evaluation)
{
    evaluator->SetEvalPoints(eval_points);
    evaluator->SetPDFValueBuffer(pdf_values);
    evaluator->SetNormalizationBuffer(norm);
    evaluator->SetParameterBuffer(params);
    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 5, *norm->readOnlyHostPtr());

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[0]));
    ASSERT_FLOAT_EQ(1.6, results[1]);
    ASSERT_FLOAT_EQ(1.6, results[2]);
    ASSERT_FLOAT_EQ(0.4, results[3]);
    ASSERT_FLOAT_EQ(0.4, results[4]);
    ASSERT_TRUE(isnan(results[5]));
}

TEST_F(EvalHistMethods, EvaluationOffsetStride)
{
    evaluator->SetEvalPoints(eval_points);
    evaluator->SetPDFValueBuffer(pdf_values, 3, 2);
    evaluator->SetNormalizationBuffer(norm, 1);
    evaluator->SetParameterBuffer(params);

    // Detect incorrect writes
    norm->writeOnlyHostPtr()[0] = 77;
    norm->writeOnlyHostPtr()[1] = 88;
    norm->writeOnlyHostPtr()[2] = 99;
    // Force flush to device
    norm->readOnlyDevicePtr();

    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 77, norm->readOnlyHostPtr()[0]);
    EXPECT_EQ((unsigned int) 5, norm->readOnlyHostPtr()[1]);
    EXPECT_EQ((unsigned int) 99, norm->readOnlyHostPtr()[2]);

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[3]));
    ASSERT_FLOAT_EQ(1.6, results[5]);
    ASSERT_FLOAT_EQ(1.6, results[7]);
    ASSERT_FLOAT_EQ(0.4, results[9]);
    ASSERT_FLOAT_EQ(0.4, results[11]);
    ASSERT_TRUE(isnan(results[13]));
}
