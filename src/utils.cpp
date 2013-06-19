#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TDirectory.h>
#include <TNtuple.h>
#include <TMath.h>
#include "utils.h"


float get_ntuple_entry(TNtuple* nt, int i, std::string field) {
  float v;
  nt->SetBranchAddress(field.c_str(), &v);
  nt->GetEvent(i);
  nt->ResetBranchAddresses();
  return v;
}


std::vector<float> get_correlation_matrix(TNtuple* nt) {
  int nentries = nt->GetEntries();

  // get list of branch names
  std::vector<std::string> names;
  for (size_t i=0; i<nt->GetListOfBranches()->GetEntries(); i++) {
    names.push_back(nt->GetListOfBranches()->At(i)->GetName());
  }

  std::vector<float> matrix(names.size() * names.size());

  // convert the ntuple to a vector, calculating means as we go
  std::vector<float> table(names.size() * nentries);
  std::vector<float> means(names.size(), 0);
  for (size_t i=0; i<nentries; i++) {
    for (size_t j=0; j<names.size(); j++) {
      float v = get_ntuple_entry(nt, i, names[j]);
      table[j + i * names.size()] = v;
      means[j] += v;
    }
  }

  // sums to means
  for (size_t i=0; i<names.size(); i++) {
    means[i] /= nentries;
  }

  // compute correlations
  for (size_t i=0; i<names.size(); i++) {
    for (size_t j=i; j<names.size(); j++) {
      float t = 0;
      float dx2 = 0;
      float dy2 = 0;
      for (int k=0; k<nentries; k++) {
        float x1 = table[k * names.size() + i] - means[i];
        float x2 = table[k * names.size() + j] - means[j];
        t += x1 * x2;
        dx2 += x1 * x1;
        dy2 += x2 * x2;
      }
      matrix[i * names.size() + j] = t / TMath::Sqrt(dx2 * dy2);
    }
  }

  return matrix;
}


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

