#include "test_pdfz_fixtures.h"
#include <cmath>

#ifndef __CUDACC__
#define isnan std::isnan
#endif

class EvalHistSystematics : public EvalHistMethods {
protected:
    virtual void SetUp() {
        EvalHistMethods::SetUp();
        evaluator->SetEvalPoints(eval_points);
        evaluator->SetPDFValueBuffer(pdf_values);
        evaluator->SetNormalizationBuffer(norm);
        evaluator->SetParameterBuffer(params);
    }
};

////////////// Shift Systematics

class EvalShiftSystematics : public EvalHistSystematics {
protected:
    virtual void SetUp() {
        EvalHistSystematics::SetUp();
        shift = new pdfz::ShiftSystematic(0, 0);
        evaluator->AddSystematic(*shift);
    }

    virtual void TearDown() {
        delete shift;
        EvalHistSystematics::TearDown();
    }

    pdfz::ShiftSystematic *shift;
};

TEST_F(EvalShiftSystematics, ZeroShift)
{
    params->writeOnlyHostPtr()[0] = 0.0;
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

TEST_F(EvalShiftSystematics, NegShift)
{
    params->writeOnlyHostPtr()[0] = -0.25;
    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 4, *norm->readOnlyHostPtr());

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[0]));
    ASSERT_FLOAT_EQ(1.5, results[1]);
    ASSERT_FLOAT_EQ(1.5, results[2]);
    ASSERT_FLOAT_EQ(0.5, results[3]);
    ASSERT_FLOAT_EQ(0.5, results[4]);
    ASSERT_TRUE(isnan(results[5]));
}

TEST_F(EvalShiftSystematics, PosShift)
{
    params->writeOnlyHostPtr()[0] = 0.25;
    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 6, *norm->readOnlyHostPtr());

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[0]));
    ASSERT_FLOAT_EQ(1.0, results[1]);
    ASSERT_FLOAT_EQ(1.0, results[2]);
    ASSERT_FLOAT_EQ(1.0, results[3]);
    ASSERT_FLOAT_EQ(1.0, results[4]);
    ASSERT_TRUE(isnan(results[5]));
}

////////////// Shift Systematics

class EvalScaleSystematics : public EvalHistSystematics {
protected:
    virtual void SetUp() {
        EvalHistSystematics::SetUp();
        scale = new pdfz::ScaleSystematic(0, 0);
        evaluator->AddSystematic(*scale);
    }

    virtual void TearDown() {
        delete scale;
        EvalHistSystematics::TearDown();
    }

    pdfz::ScaleSystematic *scale;
};

TEST_F(EvalScaleSystematics, ZeroScale)
{
    params->writeOnlyHostPtr()[0] = 0.0;
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

TEST_F(EvalScaleSystematics, NegScale)
{
    params->writeOnlyHostPtr()[0] = -0.1;
    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 6, *norm->readOnlyHostPtr());

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[0]));
    ASSERT_FLOAT_EQ(5.0/3, results[1]);
    ASSERT_FLOAT_EQ(5.0/3, results[2]);
    ASSERT_FLOAT_EQ(1.0/3, results[3]);
    ASSERT_FLOAT_EQ(1.0/3, results[4]);
    ASSERT_TRUE(isnan(results[5]));
}

TEST_F(EvalScaleSystematics, PosScale)
{
    params->writeOnlyHostPtr()[0] = 1.0; // Remember scale is 1 + 1.0!
    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 4, *norm->readOnlyHostPtr());

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[0]));
    ASSERT_FLOAT_EQ(1.0, results[1]);
    ASSERT_FLOAT_EQ(1.0, results[2]);
    ASSERT_FLOAT_EQ(1.0, results[3]);
    ASSERT_FLOAT_EQ(1.0, results[4]);
    ASSERT_TRUE(isnan(results[5]));
}
