#include <iostream>
#include <string>

#include <TRandom.h>
#include <TStopwatch.h>

#include "pdfz.h"

using namespace std;

void fill_gaussian(std::vector<float> &samples)
{
    for (unsigned int i=0; i < samples.size(); i++)
        samples[i] = gRandom->Gaus();
}

void fill_clamped_gaussian(std::vector<float> &samples, float lower, float upper)
{
    for (unsigned int i=0; i < samples.size(); i++) {
        while (true) {
            samples[i] = gRandom->Gaus();
            if (samples[i] >= lower && samples[i] < upper)
                break;
        }
    }
}


void bench_pdfz()
{
    const int nsamples = 10000000;
    const int neval_points = 100000;
    const int nbins = 1000;

    cout << "pdfz benchmark\n"
            "--------------\n"
            "Config: # of samples = " << nsamples << "\n"
         << "        # of evaluation points = " << neval_points << "\n"
         << "        # of bins = " << nbins << "\n";


    std::vector<float> lower(1);
    std::vector<float> upper(1);
    std::vector<int> nbins_vec(1);

    lower[0] = -3.0f;
    upper[0] = 3.0f;
    nbins_vec[0] = nbins;

    std::vector<float> samples(nsamples);
    fill_gaussian(samples);

    pdfz::EvalHist evaluator(samples, 1, 1, lower, upper, nbins_vec);

    // Setup for evaluation
    vector<float> eval_points(neval_points);
    fill_clamped_gaussian(eval_points, lower[0], upper[0]);

    hemi::Array<float> pdf_values(neval_points, true);
    hemi::Array<unsigned int> norm (1, true);
    hemi::Array<float> params(1, true);

    evaluator.SetEvalPoints(eval_points);
    evaluator.SetPDFValueBuffer(&pdf_values);
    evaluator.SetNormalizationBuffer(&norm);
    evaluator.SetParameterBuffer(&params);

    // Warmup
    evaluator.EvalAsync();
    evaluator.EvalFinished();

    const int nreps = 100;
    TStopwatch timer;
    timer.Start();
    for (int i=0; i < nreps; i++) {
        evaluator.EvalAsync();
        evaluator.EvalFinished();
    }
    timer.Stop();

    double samples_per_second = nsamples * nreps / timer.RealTime();
    cout << "# of samples histogrammed per second: " << samples_per_second << "\n";
}


void bench_pdfz_group()
{
    // Benchmark parameters
    const int neval_points = 100000;
    const int nbins = 1000;

    std::vector<float> lower(1);
    std::vector<float> upper(1);
    std::vector<int> nbins_vec(1);

    lower[0] = -3.0f;
    upper[0] = 3.0f;
    nbins_vec[0] = nbins;

    const int nsignals = 29;
    const int nsamples[nsignals] = { 
        1e3 /* 0vbb */,
        2e5 /* 2vbb */,
        1e4 /* b8 */,
        1e3 /* int_bi214 */,
        1e3 /* int_tl208 */,
        3e6 /* av_bi214 */,
        5e5 /* av_tl208 */,
        1e6 /* water_bi214 */,
        8e4 /* water_tl208 */,
        2e4 /* pmt_bg */,
        1e3 /* Sc-44 */,
        1e3 /* TI-44 */,
        1e3 /* Ga-68 */,
        1e3 /* Al-26 */,
        1e3 /* Rb-82 */,
        1e3 /* Sr-82 */,
        1e3 /* Y-88  */,
        1e3 /* K-42  */,
        1e3 /* Ar-42 */,
        1e3 /* Co-56 */,
        1e3 /* Co-60 */,
        1e3 /* Ag-110m */,
        1e3 /* Rh-106 */,
        1e3 /* Sn/Sb-126*/,
        1e3 /* Sb-124 */,
        1e3 /* Na-22 */,
        1e3 /* Rb-84 */,
        1e3 /* Sr-90 */,
        1e3 /* Rh-102 */,
    };

    // Banner
    cout << "pdfz group benchmark\n"
            "--------------------\n"
            "Config: # of samples =";
    int nsamples_total = 0;
    for (int i=0; i < nsignals; i++) {
        cout << " " << nsamples[i];
        nsamples_total += nsamples[i];
    }

    cout << "\n        # of evaluation points = " << neval_points << "\n"
         <<   "        # of bins = " << nbins << "\n";

    // Setup arrays for evaluators
    vector<float> eval_points(neval_points);
    fill_clamped_gaussian(eval_points, lower[0], upper[0]);

    hemi::Array<float> pdf_values(neval_points * nsignals, true);
    hemi::Array<unsigned int> norm (nsignals, true);
    hemi::Array<float> params(1, true);

    // Initialize evaluators
    pdfz::EvalHist *evaluators[nsignals];
    std::vector<float> samples;
    for (int i = 0; i < nsignals; i++) {
        samples.resize(nsamples[i]);
        fill_gaussian(samples);
        pdfz::EvalHist *evaluator = new pdfz::EvalHist(samples, 1, 1, lower, upper, nbins_vec);

        evaluator->SetEvalPoints(eval_points);
        evaluator->SetPDFValueBuffer(&pdf_values, neval_points * i);
        evaluator->SetNormalizationBuffer(&norm, i);
        evaluator->SetParameterBuffer(&params);

        // Warmup
        evaluator->EvalAsync();
        evaluator->EvalFinished();
        evaluators[i] = evaluator;
    }


    // Evaluate all the PDFs
    const int nreps = 100;
    TStopwatch timer;
    timer.Start();
    for (int i=0; i < nreps; i++) {
        for (int isig=0; isig < nsignals; isig++) {
            evaluators[isig]->EvalAsync();
        }
        for (int isig=0; isig < nsignals; isig++) {
            evaluators[isig]->EvalFinished();
        }
    }
    timer.Stop();

    double samples_per_second = nsamples_total * nreps / timer.RealTime();
    cout << "Avg # of samples histogrammed per second: " << samples_per_second << "\n";

}


int main(int argc, char **argv)
{
    if (argc != 2) {
        cerr << "Usage: bench_sxmc [benchmark_name]\n";
        cerr << "  Available benchmarks: pdfz pdfz_group\n";
        return 1;
    }

    if (string("pdfz") == argv[1])
        bench_pdfz();
    else if (string("pdfz_group") == argv[1])
        bench_pdfz_group();
    else {
        cerr << "Unknown benchmark name: " << argv[1] << "\n";
        return 1;
    }

    return 0;
}
