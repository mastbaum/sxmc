#include "test_pdfz_fixtures.h"

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

TEST(SystematicObjects, ResolutionScaleSystematicConstructor)
{
    pdfz::ResolutionScaleSystematic syst(0, 2, 1);
    EXPECT_EQ(pdfz::Systematic::RESOLUTION_SCALE, syst.type);
    EXPECT_EQ(0, syst.obs);
    EXPECT_EQ(2, syst.true_obs);
    EXPECT_EQ(1, syst.par);
}

////////////////////////////////

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

TEST_F(EvalHistMethods, CreateHistogram1D)
{
    evaluator->SetNormalizationBuffer(norm);
    evaluator->SetParameterBuffer(params);

    TH1 *hist = evaluator->CreateHistogram();

    EXPECT_EQ(2, hist->GetNbinsX());
    ASSERT_FLOAT_EQ(1.0, hist->Integral("width"));
    ASSERT_FLOAT_EQ(1.6, hist->GetBinContent(hist->FindBin(0.25)));
    ASSERT_FLOAT_EQ(0.4, hist->GetBinContent(hist->FindBin(0.75)));

    delete hist;
}
