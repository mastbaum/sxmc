#include <stdio.h>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <sxmc/osctable.h>
#include <sxmc/ttree_io.h>

OscTable::OscTable(std::string _filename, std::string _nutype, 
		   std::vector<double>& lut_pee,
		   std::vector<double>& lut_pars
		   )
  : filename(_filename), nutype(_nutype) {

  // Load data from look up table file.
  std::vector<float> data;
  std::vector<unsigned int> rank;
  std::vector<std::string> ttree_fields;

  int rc = \
    sxmc::io::read_float_vector_ttree(_filename, data, rank, ttree_fields);
  assert(rc >= 0);

  // Find the index relating to the survival probability
  // This is ugly, but works!
  int lut_index = -1;
  int theta_index = -1;
  int mass_index = -1;
  int ene_index = -1;

  for (unsigned int i = 0; i < ttree_fields.size(); i++) {
    std::size_t isHep = ttree_fields[i].find("dHepPee");
    std::size_t isB8 = ttree_fields[i].find("dB8Pee");
    std::size_t isTheta = ttree_fields[i].find("dTanSqT12");
    std::size_t isMass = ttree_fields[i].find("dDeltaMSq21");
    std::size_t isEnergy = ttree_fields[i].find("dEnu");
    if (isHep != std::string::npos && nutype == "hep")
      lut_index = i;
    if (isB8 != std::string::npos && nutype == "b8")
      lut_index = i;

    if (isTheta != std::string::npos) theta_index = i;
    if (isMass != std::string::npos) mass_index = i;
    if (isEnergy != std::string::npos) ene_index = i;
  }

  if (lut_index < 0)
    std::cerr << "OscTable: Expected survival probability is missing from tree." << std::endl;

  std::vector<double> lut_theta;
  std::vector<double> lut_mass;
  std::vector<double> lut_energy;

  // Save lut information of interest into an array.
  for (unsigned int i = 0; i < data.size(); i += ttree_fields.size() ) {
    lut_pee.push_back(data[i+lut_index]);

    lut_theta.push_back(data[i+theta_index]);
    lut_mass.push_back(data[i+mass_index]);
    lut_energy.push_back(data[i+ene_index]);
  }

  // Note that order matters!
  std::vector< std::vector<double> > lut_pars_vectors;
  lut_pars_vectors.push_back(lut_theta);
  lut_pars_vectors.push_back(lut_mass);
  lut_pars_vectors.push_back(lut_energy);

  double min = 0.;
  double max = 0.;
  double nsteps = 0.;
  double step = 0.;

  for (unsigned int i = 0; i < lut_pars_vectors.size(); i++) {
    GetParLimits(lut_pars_vectors[i], min, max, nsteps, step);
    lut_pars.push_back(min);
    lut_pars.push_back(max);
    lut_pars.push_back(nsteps);
    lut_pars.push_back(step);
  }

  std::cout << "OscTable: Loading " << nutype << " P_ee values: " << lut_pee.size() << std::endl;
}

void OscTable::GetParLimits(std::vector<double> inputVector, double &min, double &max, double &nsteps, double &step) {
  std::sort(inputVector.begin(), inputVector.end());
  inputVector.erase(std::unique(inputVector.begin(), inputVector.end()), inputVector.end());

  min = inputVector.front();
  max = inputVector.back();
  nsteps = inputVector.size();
  step = (max - min)/ (nsteps - 1.);
  //std::cout << min << " " << max << " " << nsteps << " " << step << std::endl;
}
