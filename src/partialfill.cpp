#include <math.h>
#include <string>
#include <TH1.h>

#include <sxmc/utils.h>
#include <sxmc/partialfill.h>

PartialFill::PartialFill(double _start_fill, double _end_fill,
    double _live_time, double _detector_radius, int _time_bins) : 
  start_fill(_start_fill), end_fill(_end_fill), live_time(_live_time), detector_radius(_detector_radius), time_bins(_time_bins)
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

  TH1D* timeh1 = new TH1D("th1_scint","th1_scint",this->time_bins,0,this->live_time);
  double full_volume = GetVolumeAtFill(1.0);
  for (int i=1;i <= timeh1->GetNbinsX();i++){
    double fill = GetFillAtTime(timeh1->GetBinCenter(i));
    double volume = GetVolumeAtFill(fill);
    timeh1->SetBinContent(i,volume/full_volume*timeh1->GetBinWidth(i));
  }
  this->time_profiles["scintillator"] = timeh1;
}

double PartialFill::GetVolumeAtFill(double fill){
  return 4/3.0*3.14159265*pow(this->detector_radius,3)*pow(fill,2)*(3-2*fill);
}

double PartialFill::GetTimeAtFill(double fill){
  double volume = GetVolumeAtFill(fill);
  return (volume-this->start_volume)/this->fill_rate;
}

