#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>

#include <sxmc/plots.h>

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
  this->c->SaveAs(filename.c_str());
}


TH1* SpectralPlot::make_like(TH1* h, std::string name) {
  TH1* hnew = dynamic_cast<TH1*>(h->Clone(name.c_str()));
  hnew->Reset();
  return hnew;
}

