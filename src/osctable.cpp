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

int OscTable::FindThetaIndex(double inTheta) {
  // Fixed values
  // Based on input file information
  double tansqt12 = 0.446;
  double tansqt12_up = 0.030;
  double tansqt12_down = 0.029;
  
  // Number of steps in lookup table
  int iTanSteps = 100;
  
  // Define bin limits
  double lowTheta = tansqt12 - 5* tansqt12_down;
  double stepSize = (5 * (tansqt12_up + tansqt12_down) / iTanSteps);
  
  // Calculate the index
  double stepNum = (inTheta - lowTheta) / stepSize;
  
  int index = (int)stepNum;
  if (stepNum - index > 0.5) index++;
  
  // Special cases for indices beyond the lookup table
  if (index < 0) index = 0;
  if (index > iTanSteps) index = iTanSteps;
  
  return index;
}

int OscTable::FindMassIndex(double inMass) {
  // Fixed values
  // Based on input file information
  double dm21sq = 7.41e-05;
  double dm21sq_up = 0.21e-05;
  double dm21sq_down = 0.18e-05;
  
  // Number of steps in lookup table
  int iMassSteps = 100;
  
  // Define bin limits
  double lowMass = dm21sq - 5* dm21sq_down;
  double stepSize = (5 * (dm21sq_up + dm21sq_down) / iMassSteps);
  
  // Calculate the index
  double stepNum = (inMass - lowMass) / stepSize;
  
  int index = (int)stepNum;
  if (stepNum - index > 0.5) index++;
  
  // Special cases for indices beyond the lookup table
  if (index < 0) index = 0;
  if (index > iMassSteps) index = iMassSteps;
  
  return index;
}

int OscTable::FindEnergyIndex(double inE) {
  
  // Number of steps in lookup table
  // Based on input file information
  int iESteps = 100;
  
  // Define bin limits
  double lowEnergy = 10.;
  double stepSize = 10. / iESteps;
  
  // Calculate the index
  double stepNum = (inE - lowEnergy) / stepSize;
  
  int index = (int)stepNum;
  if (stepNum - index > 0.5) index++;
  
  // Special cases for indices beyond the lookup table
  if (index < 0) index = 0;
  if (index > iESteps) index = iESteps;
  
  return index;
}

int OscTable::CalculateTotalIndex(int iThetaIndex, int iMassIndex, int iEIndex) {
  // Based on input file information (101 bins for each axis)
  int outIndex = iThetaIndex*101*101 + iMassIndex*101 + iEIndex;
  return outIndex;
}

double OscTable::Reweight(double dTheta, double dMass, double dEnergy) {
  int iTheta  = OscTable::FindThetaIndex(dTheta);
  int iMass  = OscTable::FindMassIndex(dMass);
  int iEnergy  = OscTable::FindEnergyIndex(dEnergy);

  int iIndex = OscTable::CalculateTotalIndex(iTheta, iMass, iEnergy);

  return (double)iIndex;
}
