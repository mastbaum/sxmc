#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TNtuple.h>
#include <TLegend.h>

#include <sxmc/plots.h>

const int colors[28] = {kRed,      kGreen,    kBlue,      kMagenta,
                        kCyan,     kYellow,   kOrange,    kViolet+2,
                        kRed+2,    kGreen+2,  kBlue+2,    kMagenta+2,
                        kCyan+2,   kYellow+2, kOrange+2,  kRed-7,
                        kGreen-7,  kBlue-7,   kMagenta-7, kCyan-7,
                        kYellow-7, kOrange-7, kRed-6,     kAzure+1,
                        kTeal+1,   kSpring-9, kAzure-9};



SpectralPlot::SpectralPlot(int _line_width, float _xmin, float _xmax,
                           float _ymin, float _ymax, bool _logy,
                           std::string _title,
                           std::string _xtitle,
                           std::string _ytitle)
    : logy(_logy), line_width(_line_width),
      xmin(_xmin), xmax(_xmax), ymin(_ymin), ymax(_ymax),
      title(_title), xtitle(_xtitle), ytitle(_ytitle) {
  this->c = new TCanvas();

  if (this->logy) {
    this->c->SetLogy();
  }

  this->legend = new TLegend(0.85, 0.15, 0.985, 0.95);
  this->legend->SetFillColor(kWhite);
}


SpectralPlot::SpectralPlot(const SpectralPlot& o) {
  this->logy = o.logy;
  this->line_width = o.line_width;
  this->xmin = o.xmin;
  this->xmax = o.xmax;
  this->ymin = o.ymin;
  this->ymax = o.ymax;
  this->title = o.title;
  this->xtitle = o.xtitle;
  this->ytitle = o.ytitle;
  for (size_t i=0; i<o.histograms.size(); i++) {
    this->histograms.push_back(o.histograms[i]);
  }
  this->c = new TCanvas();

  if (o.logy) {
    this->c->SetLogy();
  }
  this->legend = (TLegend*) o.legend->Clone("");
}


SpectralPlot::~SpectralPlot() {
  for (size_t i=0; i<this->histograms.size(); i++) {
    delete this->histograms[i];
  }
  delete this->legend;
}


void SpectralPlot::add(TH1* _h, std::string title, std::string options) {
  std::string name = "__" + std::string(title);
  TH1* h = dynamic_cast<TH1*>(_h->Clone(name.c_str()));
  h->SetDirectory(NULL);

  h->SetLineWidth(this->line_width);
  h->SetTitle(this->title.c_str());
  h->SetXTitle(xtitle.c_str());
  h->SetYTitle(ytitle.c_str());

  this->legend->AddEntry(h, title.c_str());

  if (!h || h->Integral() == 0) {
    return;
  }

  this->histograms.push_back(h);

  if (this->histograms.size() == 1) {
    h->SetAxisRange(this->ymin, this->ymax, "Y");
    h->SetAxisRange(this->xmin, this->xmax, "X");
    h->GetXaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    if (this->logy) {
      this->c->SetLogy();
    }
    this->c->cd();
    h->DrawClone(options.c_str());
  }
  else {
    this->c->cd();
    h->DrawClone(("same " + options).c_str());
  }
  this->c->Update();
}


void SpectralPlot::save(std::string filename) {
  this->c->cd();
  this->legend->SetTextFont(132);
  this->legend->Draw();
  this->c->Update();
  this->c->SaveAs(filename.c_str(), "q");
}


TH1* SpectralPlot::make_like(TH1* h, std::string name) {
  TH1* hnew = dynamic_cast<TH1*>(h->Clone(name.c_str()));
  hnew->Reset();
  return hnew;
}


void plot_fit(std::map<std::string, Interval> best_fit, float live_time,
              std::vector<Signal> signals,
              std::vector<Systematic> systematics,
              std::vector<Observable> observables,
              std::vector<float> data,
              std::string output_path) {
  std::vector<SpectralPlot> plots_full;
  std::vector<SpectralPlot> plots_external;
  std::vector<SpectralPlot> plots_cosmogenic;

  // Set up plots for each observable
  for (size_t i=0; i<observables.size(); i++) {
    Observable* o = &observables[i];
    std::stringstream ytitle;
    ytitle << "Counts/" << std::setprecision(3)
           << (o->upper - o->lower) / o->bins << " " << o->units
           << "/" << live_time << " y";
    plots_full.push_back(SpectralPlot(2, o->lower, o->upper, 1e-2, 1e6,
                         true, "", o->title, ytitle.str().c_str()));
    plots_external.push_back(SpectralPlot(2, o->lower, o->upper, 1e-2, 1e6,
                             true, "", o->title, ytitle.str().c_str()));
    plots_cosmogenic.push_back(SpectralPlot(2, o->lower, o->upper, 1e-2, 1e6,
                               true, "", o->title, ytitle.str().c_str()));
  }

  std::vector<TH1D*> external_total(observables.size(), NULL);
  std::vector<TH1D*> cosmogenic_total(observables.size(), NULL);
  std::vector<TH1D*> fit_total(observables.size(), NULL);

  // Extract best-fit parameter values
  std::vector<float> params;
  for (size_t i=0; i<signals.size(); i++) {
    params.push_back(best_fit[signals[i].name].point_estimate);
  }
  for (size_t i=0; i<systematics.size(); i++) {
    params.push_back(best_fit[systematics[i].name].point_estimate);
  }

  hemi::Array<unsigned> norms_buffer(signals.size(), true);
  norms_buffer.writeOnlyHostPtr();

  hemi::Array<double> param_buffer(params.size(), true);
  for (size_t i=0; i<params.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = params[i];
  }

  for (size_t i=0; i<signals.size(); i++) {
    pdfz::EvalHist* phist = \
      dynamic_cast<pdfz::EvalHist*>(signals[i].histogram);

    phist->SetParameterBuffer(&param_buffer, signals.size());
    phist->SetNormalizationBuffer(&norms_buffer, i);

    TH1* hpdf_nd = phist->CreateHistogram();
    hpdf_nd->Scale(params[i] / hpdf_nd->Integral());

    std::vector<TH1D*> hpdf(observables.size(), NULL);
    if (hpdf_nd->IsA() == TH1D::Class()) {
      hpdf[0] = dynamic_cast<TH1D*>(hpdf_nd);
    }
    else if (hpdf_nd->IsA() == TH2D::Class()) {
      hpdf[0] = dynamic_cast<TH2D*>(hpdf_nd)->ProjectionX("hpdf_x");
      hpdf[1] = dynamic_cast<TH2D*>(hpdf_nd)->ProjectionY("hpdf_y");
    }
    else if (hpdf_nd->IsA() == TH3D::Class()) {
      hpdf[0] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionX("hpdf_x");
      hpdf[1] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionY("hpdf_y");
      hpdf[2] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionZ("hpdf_y");
    }

    std::string n = signals[i].name;

    for (size_t j=0; j<observables.size(); j++) {
      hpdf[j]->SetLineColor(colors[i % 28]);
      if (fit_total[j] == NULL) {
        std::string hfname = "fit_total_" + signals[i].name;
        fit_total[j] = (TH1D*) hpdf[j]->Clone(hfname.c_str());
      }
      else {
        if (hpdf[j] && hpdf[j]->Integral() > 0) {
          fit_total[j]->Add(hpdf[j]);
        }
      }

      if (n == "av_tl208" || n == "av_bi214" ||
          n == "water_tl208" || n == "water_bi214" ||
          n == "int_ropes_tl208" || n == "int_ropes_bi214" ||
          n == "hd_ropes_tl208" || n == "hd_ropes_bi214" ||
          n == "pmt_bg") {
        plots_external[j].add(hpdf[j], signals[i].title, "hist");
        std::string hname = "et" + signals[i].name + observables[j].name;
        if (external_total[j] == NULL) {
          external_total[j] = (TH1D*) hpdf[j]->Clone(hname.c_str());
        }
        else {
          if (hpdf[j] && hpdf[j]->Integral() > 0) {
            external_total[j]->Add((TH1D*) hpdf[j]->Clone(hname.c_str()));
          }
        }
      }
      else if (n != "zeronu" && n != "twonu" && n != "int_tl208" &&
               n != "int_bi214" && n != "b8") {
        plots_cosmogenic[j].add(hpdf[j], signals[i].title, "hist");
        std::string hname = "ct" + signals[i].name + observables[j].name;
        if (cosmogenic_total[j] == NULL) {
          cosmogenic_total[j] = (TH1D*) hpdf[j]->Clone(hname.c_str());
        }
        else {
          cosmogenic_total[j]->Add((TH1D*) hpdf[j]->Clone(hname.c_str()));
        }
      }
      else {
        if (hpdf[j] && hpdf[j]->Integral() > 0) {
          plots_full[j].add(hpdf[j], signals[i].title, "hist");
        }
      }
    }
  }

  for (size_t i=0; i<observables.size(); i++) {
    TH1D* hdata = \
      (TH1D*) SpectralPlot::make_like(plots_full[i].histograms[0], "hdata");

    hdata->SetMarkerStyle(20);
    hdata->SetLineColor(kBlack);

    for (size_t idata=0; idata<data.size() / observables.size(); idata++) {
      hdata->Fill(data[idata * observables.size() + i]);
    }

    if (external_total[i] != NULL) {
      external_total[i]->SetLineColor(kOrange + 1);
      TH1D* et = (TH1D*) external_total[i]->Clone("et");
      et->SetLineStyle(2);
      plots_external[i].add(et, "Total", "hist");
      plots_full[i].add(external_total[i], "External", "hist");
    }
    if (cosmogenic_total[i] != NULL) {
      cosmogenic_total[i]->SetLineColor(kAzure + 1);
      TH1D* ct = (TH1D*) cosmogenic_total[i]->Clone("ct");
      ct->SetLineStyle(2);
      plots_cosmogenic[i].add(ct, "Total", "hist");
      plots_full[i].add(cosmogenic_total[i], "Cosmogenic", "hist");
    }
    if (fit_total[i] != NULL) {
      fit_total[i]->SetLineColor(kRed);
      plots_full[i].add(fit_total[i], "Fit", "hist");
    }

    plots_full[i].add(hdata, "Fake Data");

    plots_full[i].save(output_path + observables[i].name + "_spectrum_full.pdf");
    plots_external[i].save(output_path + observables[i].name + "_spectrum_external.pdf");
    plots_cosmogenic[i].save(output_path + observables[i].name +
                             "_spectrum_cosmogenic.pdf");
  }
}

