#include <math.h>
#include <string>
#include <TH1.h>
#include <TFile.h>

#include <sxmc/utils.h>
#include <sxmc/partialfill.h>

PartialFill::PartialFill(double _start_fill, double _end_fill,
    double _live_time, double _detector_radius, int _time_bins, double _scint_leaching_rate, double _water_leaching_rate) : 
  start_fill(_start_fill), end_fill(_end_fill), live_time(_live_time), detector_radius(_detector_radius), time_bins(_time_bins), scint_leaching_rate(_scint_leaching_rate), water_leaching_rate(_water_leaching_rate)
{
  this->start_volume = GetVolumeAtFill(this->start_fill);
  this->end_volume = GetVolumeAtFill(this->end_fill);
  this->fill_rate = (this->end_volume-this->start_volume)/this->live_time;

  double delta_fill = 0.001;

  std::vector<double> x, y;
  for (double fill=this->start_fill;fill < this->end_fill + delta_fill; fill+= delta_fill){
    y.push_back(fill);
    x.push_back(GetTimeAtFill(fill));
  }

  this->GetFillAtTime.Setup(x,y);

  TH1D* tscint = new TH1D("th1_scint","th1_scint",this->time_bins,0,this->live_time);
  double full_volume = GetVolumeAtFill(1.0);
  for (int i=1;i <= tscint->GetNbinsX();i++){
    double fill = GetFillAtTime(tscint->GetBinCenter(i));
    double volume = GetVolumeAtFill(fill);
    tscint->SetBinContent(i,volume/full_volume*tscint->GetBinWidth(i));
  }
  this->time_profiles["scintillator"] = tscint;

  TH1D* tinwater = new TH1D("th1_inwater","th1_inwater",this->time_bins,0,this->live_time);
  for (int i=1;i <= tinwater->GetNbinsX();i++)
    tinwater->SetBinContent(i,tscint->GetBinWidth(i)-tscint->GetBinContent(i));
  this->time_profiles["inwater"] = tinwater;

  std::vector<double> times;
  std::vector<double> waterrate;
  std::vector<double> scintrate;
  std::vector<double> avrate;
  
  avrate.push_back(1.0);
  waterrate.push_back(0.0);
  scintrate.push_back(0.0);

  double deltat = this->live_time/this->time_bins/10.0;
  double time = deltat;
  while (time <= this->live_time+deltat){
    times.push_back(time);
    double scint_frac = GetSurfaceAreaAtFill(GetFillAtTime(time))/GetSurfaceAreaAtFill(1.0);
    double water_frac = 1.0-scint_frac;
    double scint_leached = scint_leaching_rate * scint_frac * deltat * avrate.back();
    double water_leached = water_leaching_rate * water_frac * deltat * avrate.back();
    avrate.push_back(avrate.back()-scint_leached-water_leached);
    scintrate.push_back(scintrate.back()+scint_leached);
    waterrate.push_back(waterrate.back()+water_leached);
    time += deltat;
  }
  LinearInterpolator iav(times,avrate);
  LinearInterpolator iwater(times,waterrate);
  LinearInterpolator iscint(times,scintrate);
  TH1D* tavleaching = new TH1D("th1_avleaching","th1_avleaching",this->time_bins,0,this->live_time);
  for (int i=1;i <= tavleaching->GetNbinsX();i++)
    tavleaching->SetBinContent(i,iav(tavleaching->GetBinCenter(i))*tavleaching->GetBinWidth(i));
  this->time_profiles["avleaching"] = tavleaching;
  TH1D* tscintleaching = new TH1D("th1_scintleaching","th1_scintleaching",this->time_bins,0,this->live_time);
  for (int i=1;i <= tscintleaching->GetNbinsX();i++)
    tscintleaching->SetBinContent(i,iscint(tscintleaching->GetBinCenter(i))*tscintleaching->GetBinWidth(i));
  this->time_profiles["scintleaching"] = tscintleaching;
  TH1D* twaterleaching = new TH1D("th1_waterleaching","th1_waterleaching",this->time_bins,0,this->live_time);
  for (int i=1;i <= twaterleaching->GetNbinsX();i++)
    twaterleaching->SetBinContent(i,iwater(twaterleaching->GetBinCenter(i))*twaterleaching->GetBinWidth(i));
  this->time_profiles["waterleaching"] = twaterleaching;
  
  TFile f("time_profiles.root","RECREATE");
  tscint->Write();
  tinwater->Write();
  tavleaching->Write();
  twaterleaching->Write();
  tscintleaching->Write();
  f.Close();

}

double PartialFill::GetVolumeAtFill(double fill){
  return 4/3.0*3.14159265*pow(this->detector_radius,3)*pow(fill,2)*(3-2*fill);
}

double PartialFill::GetTimeAtFill(double fill){
  double volume = GetVolumeAtFill(fill);
  return (volume-this->start_volume)/this->fill_rate;
}

double PartialFill::GetSurfaceAreaAtFill(double fill){
  return 4.0*3.14159265*pow(this->detector_radius,2)*fill;
}
