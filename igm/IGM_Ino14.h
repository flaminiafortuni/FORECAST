#include <cmath>
#include <armadillo>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

class IGM
{

public:
  IGM() = default;
  
  double tau_IGM(long double lobs, float z);

  double tau_LAF_LS(long double lobs, float z);

  double tau_DLA_LS(long double lobs, float z);

  double tau_DLA_LC(long double lobs, float z);

  double tau_LAF_LC(long double lobs, float z);

  //double *igm_absorption(double *lambda, int n_lambda, double z);

  void igm_absorption(float zr, std::vector<long double> &wave, std::vector<long double> &sed);

  // Destructor
  ~IGM() = default;

private:


  std::vector<long double> tsed;


};
