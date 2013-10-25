#ifndef __PARTIALFILL_H__
#define __PARTIALFILL_H__

#include <string>
#include <map>
#include <TH1.h>

#include <sxmc/utils.h>

/**
 * \class PartialFill
 * \brief Manages the configuration of the partial fill related parameters
 *
 * Could be part of the main FitConfig but I figured I'd isolate the
 * parts only suited for the partial fill analysis
 */

class PartialFill {
  public:
    PartialFill(double _start_fill, double _end_fill,
        double _live_time, double _detector_radius, int _time_bins);

    double GetVolumeAtFill(double fill);
    double GetTimeAtFill(double fill);

    LinearInterpolator GetFillAtTime;

    TH1D* GetTimeProfile(std::string location){return time_profiles[location];};

    int GetTimeBins(){return time_bins;};

  protected:
    double start_fill, end_fill, live_time, start_volume, end_volume;
    double fill_rate;
    double detector_radius;
    int time_bins;
    std::map<std::string, TH1D*> time_profiles;
};

#endif
