#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TDirectory.h>
#include "utils.h"

SpectralPlot::SpectralPlot(int _line_width, float _xmin, float _xmax,
                           float _ymin, float _ymax, bool _logy,
                           std::string _title,
                           std::string _xtitle,
                           std::string _ytitle)
  : first(true), logy(_logy), line_width(_line_width), xmin(_xmin),
    xmax(_xmax), ymin(_ymin), ymax(_ymax), title(_title), xtitle(_xtitle),
    ytitle(_ytitle) {
  if (this->logy) {
    this->c.SetLogy();
  }

  this->legend = new TLegend(0.825, 0.15, 0.985, 0.95);
  this->legend->SetFillColor(kWhite);
}


SpectralPlot::~SpectralPlot() {
  for (size_t i=0; i<this->histograms.size(); i++) {
    delete this->histograms[i];
  }
  delete this->legend;
}


void SpectralPlot::add(TH1* _h, std::string title, std::string options) {
  TH1* h = (TH1*) _h->Clone(("__" + std::string(_h->GetName())).c_str());

  h->SetDirectory(NULL);
  h->SetLineWidth(this->line_width);
  h->SetTitle(this->title.c_str());
  h->SetXTitle(xtitle.c_str());
  h->SetYTitle(ytitle.c_str());

  this->legend->AddEntry(h, title.c_str());

  this->histograms.push_back(h);
  this->c.cd();

  if (this->first) {
    h->Draw(options.c_str());
    this->first = false;
  }
  else {
    h->Draw(("same " + options).c_str());
  }
}


void SpectralPlot::save(std::string filename) {
  this->c.cd();
  this->legend->Draw();
  this->c.SaveAs(filename.c_str());
}


TH1* SpectralPlot::make_like(TH1* h, std::string name) {
  TH1* hnew = (TH1*) h->Clone(name.c_str());
  hnew->Reset();
  return hnew;
}

