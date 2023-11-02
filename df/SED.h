#include <cmath>
#include <armadillo>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

class SED
{

public:
  SED() = default;

  void extract_spec(std::vector<std::vector<double> > &full_table, std::vector<double> &time_grid, int a_indx, float ages4, std::vector<long double> &spe);
  
  double getY(std:: vector<double> x, std:: vector<double> y,double xi);

  double compute_mab(float zr, std::vector <long double> &waves,std::vector <long double> &sed, std::vector <long double> &fwaves,std::vector <long double> &fresp);

  void z_evol(float zr, std::vector<long double> &wave, std::vector<long double> &sed, std::vector<double> &zl, std::vector<double> &dlum);//double dlum

  void dust_attenuation( std::vector <float> &ex_l,std::vector <long double> &dustAtt, std::vector <long double> &wavesc,std::vector <long double> &spe);

  template <class T> int locate(const std::vector<T> &v, const T x);


  // Destructor
  ~SED() = default;

private:

//interp_2spec variables
  double a_1;
  double a_2;
  double t_1;
  double t_2;
  double m_1;
  double m_2;
  std::vector<double> f_1;
  std::vector<double> f_2;
  const double ergsa = 3.9e+33;
 
//getY function variables
  int nn;
  int i;
  double f;
  double a0;
  double a1; 
  double a2;
  double a3;
  double f2;
  size_t n;
  int jl;
  int ju;
  bool as;
  int jm;

//z_evol variables
  double dlsmpc;
  double dls;
  double mpctocm=3.086e+24;
  const double speedcunitas = 2.9979e+18;
  int in_Lmin=0, in_Lmax=0;
  double h_step=0.;
  std::vector<double> wave_i;
  std::vector<double>::iterator x;
  double val;
  typedef std::vector<long double> stdvec;
  typedef std::vector<double> stdvecd;
  arma::vec nsed;
  arma::vec nwave;
  arma::vec sed_interp;
  arma::vec wave_interp;
  std::vector<long double> twave;
  std::vector<long double> tsed;
  std::vector<long double> tespo;

//compute_mab variables
  arma::vec cfresp;
  arma::vec cfwaves;
  double lminf;
  double lmaxf;
  std::vector<long double> newwaves;
  std::vector<long double> newsed;
  arma::vec csed;
  arma::vec cwaves;
  arma::vec filter_interp;
  arma::vec filterSpec;
  arma::vec i1;
  arma::vec i2; 
  arma::mat I1;  				  
  arma::mat I2;
  double flambda;
  double fnu;
  double mAB;

  std::vector<double> newwaves_i;
  std::vector<double>::iterator x_2;
  double val_2;
  arma::vec nsed_2;
  arma::vec nwave_2;
  arma::vec sed_interp_2;
  arma::vec wave_interp_2;

  //dust_attenuation variables
  arma::vec dustAtt_interp;
  arma::vec ex_l_interp;
  arma::vec a_spe;
  arma::vec a_ex_l;
  arma::vec a_dustAtt;
  arma::vec out_spe;
  typedef std::vector<float> stdvec2;

};
