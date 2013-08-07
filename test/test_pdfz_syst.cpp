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

////////////// Scale Systematics

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

////////////// Resolution Scale Systematics

class EvalResolutionScaleSystematics : public EvalHistConstructor {
protected:
    virtual void SetUp() {
        EvalHistConstructor::SetUp();

        // Change sampes to have a true energy field
        nfields = 2;
        samples.resize(7*2);
        samples[0] = 0.1; samples[1] = 0.7;
        samples[2] = 0.2; samples[3] = 0.7;
        samples[4] = 0.3; samples[5] = 0.7;
        samples[6] = 0.4; samples[7] = 0.7;
        samples[8] = 0.5; samples[9] = 0.7;
        samples[10] = 1.1; samples[11] = 0.7;
        samples[12] = -0.1; samples[13] = 0.7;

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
        params = new hemi::Array<double>(5, true);
        params->writeOnlyHostPtr(); // Force memory allocation

        evaluator->SetEvalPoints(eval_points);
        evaluator->SetPDFValueBuffer(pdf_values);
        evaluator->SetNormalizationBuffer(norm);
        evaluator->SetParameterBuffer(params);

        res = new pdfz::ResolutionScaleSystematic(0, 1, 0);
        evaluator->AddSystematic(*res);
    }

    virtual void TearDown() {
        delete res;
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
    hemi::Array<double> *params;
    pdfz::ResolutionScaleSystematic *res;
};

TEST_F(EvalResolutionScaleSystematics, ZeroScale)
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

TEST_F(EvalResolutionScaleSystematics, NegScale)
{
    params->writeOnlyHostPtr()[0] = -0.30;
    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 7, *norm->readOnlyHostPtr());

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[0]));
    ASSERT_FLOAT_EQ(2.0*5/7, results[1]);
    ASSERT_FLOAT_EQ(2.0*5/7, results[2]);
    ASSERT_FLOAT_EQ(2.0*2/7, results[3]);
    ASSERT_FLOAT_EQ(2.0*2/7, results[4]);
    ASSERT_TRUE(isnan(results[5]));
}

TEST_F(EvalResolutionScaleSystematics, PosScale)
{
    params->writeOnlyHostPtr()[0] = 0.30;
    evaluator->EvalAsync();
    evaluator->EvalFinished();

    EXPECT_EQ((unsigned int) 4, *norm->readOnlyHostPtr());

    float *results = pdf_values->hostPtr();
    ASSERT_TRUE(isnan(results[0]));
    ASSERT_FLOAT_EQ(2.0, results[1]);
    ASSERT_FLOAT_EQ(2.0, results[2]);
    ASSERT_FLOAT_EQ(0.0, results[3]);
    ASSERT_FLOAT_EQ(0.0, results[4]);
    ASSERT_TRUE(isnan(results[5]));
}


