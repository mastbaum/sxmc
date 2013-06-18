#include <iostream>
#include <string>

#include <TRandom.h>
#include <TStopwatch.h>

#include "pdfz.h"

using namespace std;

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
    for (int i=0; i < nsamples; i++)
        samples[i] = gRandom->Gaus();

    pdfz::EvalHist evaluator(samples, 1, 1, lower, upper, nbins_vec);

    // Setup for evaluation
    vector<float> eval_points(neval_points);
    for (int i=0; i < neval_points; i++) {
        while (true) {
            eval_points[i] = gRandom->Gaus();
            if (eval_points[i] >= lower[0] && eval_points[i] < upper[0])
                break;
        }
    }

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

    double evals_per_second = neval_points * nreps / timer.RealTime();
    cout << "PDF evaluations/second: " << evals_per_second << "\n";
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        cerr << "Usage: bench_sxmc [benchmark_name]\n";
        cerr << "  Available benchmarks: pdfz\n";
        return 1;
    }

    if (string("pdfz") == argv[1])
        bench_pdfz();

    return 0;
}