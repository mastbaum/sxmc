#include <cmath>
#include "test_pdfz_fixtures_2d.h"

#ifndef __CUDACC__
#define isnan std::isnan
#endif

TEST_F(EvalHist2DConstructor, WrongSampleSize) {
  ASSERT_THROW(pdfz::EvalHist(samples, weights, 3 /* nfields */, nobservables, lower, upper, nbins), pdfz::Error);
}


TEST_F(EvalHist2DConstructor, NobsLargerThanNfields) {
  ASSERT_THROW(pdfz::EvalHist(samples, weights, nfields , 7 /* nobservables */, lower, upper, nbins), pdfz::Error);
}


TEST_F(EvalHist2DConstructor, WrongLowerSize) {
  lower.resize(1);
  ASSERT_THROW(pdfz::EvalHist(samples, weights, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}


TEST_F(EvalHist2DConstructor, WrongUpperSize) {
  upper.resize(1);
  ASSERT_THROW(pdfz::EvalHist(samples, weights, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}


TEST_F(EvalHist2DConstructor, WrongNbinsSize) {
  nbins.resize(1);
  ASSERT_THROW(pdfz::EvalHist(samples, weights, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}


TEST_F(EvalHist2DConstructor, ZeroBins) {
  nbins[1] = 0;
  ASSERT_THROW(pdfz::EvalHist(samples, weights, nfields, nobservables, lower, upper, nbins), pdfz::Error);
}


///////////////


TEST_F(EvalHist2DMethods, Evaluation) {
  evaluator->SetEvalPoints(eval_points);
  evaluator->SetPDFValueBuffer(pdf_values);
  evaluator->SetNormalizationBuffer(norm);
  evaluator->SetParameterBuffer(params);

  evaluator->EvalAsync();
  evaluator->EvalFinished();

  EXPECT_EQ((unsigned int) 6, *norm->readOnlyHostPtr());

  // Number of samples in boundary * bin area
  const double norm = 6 * (0.5 * (2.0/3.0));

  float* results = pdf_values->hostPtr();
  ASSERT_FLOAT_EQ(1/norm, results[0]);
  ASSERT_FLOAT_EQ(0/norm, results[1]);
  ASSERT_FLOAT_EQ(2/norm, results[2]);
  ASSERT_FLOAT_EQ(0/norm, results[3]);
  ASSERT_FLOAT_EQ(3/norm, results[4]);
  ASSERT_TRUE(isnan(results[6]));
  ASSERT_TRUE(isnan(results[7]));
  ASSERT_TRUE(isnan(results[8]));
}


TEST_F(EvalHist2DMethods, EvaluationOffsetStride) {
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
  EXPECT_EQ((unsigned int) 6, norm->readOnlyHostPtr()[1]);
  EXPECT_EQ((unsigned int) 99, norm->readOnlyHostPtr()[2]);

  // Number of samples in boundary * bin area
  const double norm = 6 * (0.5 * (2.0/3.0));

  float* results = pdf_values->hostPtr();
  ASSERT_FLOAT_EQ(1/norm, results[3]);
  ASSERT_FLOAT_EQ(0/norm, results[5]);
  ASSERT_FLOAT_EQ(2/norm, results[7]);
  ASSERT_FLOAT_EQ(0/norm, results[9]);
  ASSERT_FLOAT_EQ(3/norm, results[11]);
  ASSERT_TRUE(isnan(results[13]));
  ASSERT_TRUE(isnan(results[15]));
  ASSERT_TRUE(isnan(results[17]));
}


TEST_F(EvalHist2DMethods, CreateHistogram2D) {
  evaluator->SetNormalizationBuffer(norm);
  evaluator->SetParameterBuffer(params);

  TH1* hist = evaluator->CreateHistogram();

  EXPECT_EQ(2, hist->GetNbinsX());
  EXPECT_EQ(3, hist->GetNbinsY());
  ASSERT_FLOAT_EQ(1.0, hist->Integral("width"));

  // Number of samples in boundary * bin area
  const double norm = 6 * (0.5 * (2.0/3.0));

  ASSERT_FLOAT_EQ(1 / norm, hist->GetBinContent(hist->FindBin(0.25, 10.5)));
  ASSERT_FLOAT_EQ(0 / norm, hist->GetBinContent(hist->FindBin(0.75, 10.5)));
  ASSERT_FLOAT_EQ(0 / norm, hist->GetBinContent(hist->FindBin(0.25, 11.0)));
  ASSERT_FLOAT_EQ(2 / norm, hist->GetBinContent(hist->FindBin(0.75, 11.0)));
  ASSERT_FLOAT_EQ(0 / norm, hist->GetBinContent(hist->FindBin(0.25, 11.5)));
  ASSERT_FLOAT_EQ(3 / norm, hist->GetBinContent(hist->FindBin(0.75, 11.5)));

  delete hist;
}

