#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include"IGM_Ino14.h"

//*# Part of Credits to Brandt Robertson

//takes in input waves [AA] already redshifted, and redshift;
//gives in output

double IGM::tau_IGM(long double lobs, float z)
{
  double tau;
  //eqn 15 of inoue et al. 2014
  tau = tau_LAF_LS(lobs, z) + tau_DLA_LS(lobs, z) + tau_LAF_LC(lobs, z) + tau_DLA_LC(lobs, z);
  //printf("%e\t%e\t%e\n",lobs,z,tau);
  return tau;
}

double IGM::tau_LAF_LS(long double lobs, float z) //checked
{
  int j;
  double tau  = 0;
  double tauj = 0;
  int nj = 39;
  //table 2 of inoue et al. 2014
  double lj[39]  = {1215.67, 1025.72, 972.537, 949.743, 937.803, 930.748, 926.226, 923.150, 920.963, 919.352, 918.129, 917.181, 916.429, 915.824, 915.329, 914.919, 914.576, 914.286, 914.039, 913.826, 913.641, 913.480, 913.339, 913.215, 913.104, 913.006, 912.918, 912.839, 912.768, 912.703, 912.645, 912.592, 912.543, 912.499, 912.458, 912.420, 912.385, 912.353, 912.324};
  double Aj1[39] = {1.690e-2,4.692e-3,2.239e-3,1.319e-3,8.707e-4,6.178e-4,4.609e-4,3.569e-4,2.843e-4,2.318e-4,1.923e-4,1.622e-4,1.385e-4,1.196e-4,1.043e-4,9.174e-5,8.128e-5,7.251e-5,6.505e-5,5.868e-5,5.319e-5,4.843e-5,4.427e-5,4.063e-5,3.738e-5,3.454e-5,3.199e-5,2.971e-5,2.766e-5,2.582e-5,2.415e-5,2.263e-5,2.126e-5,2.000e-5,1.885e-5,1.779e-5,1.682e-5,1.593e-5,1.510e-5};
  double Aj2[39] = {2.354e-3,6.536e-4,3.119e-4,1.837e-4,1.213e-4,8.606e-5,6.421e-5,4.971e-5,3.960e-5,3.229e-5,2.679e-5,2.259e-5,1.929e-5,1.666e-5,1.453e-5,1.278e-5,1.132e-5,1.010e-5,9.062e-6,8.174e-6,7.409e-6,6.746e-6,6.167e-6,5.660e-6,5.207e-6,4.811e-6,4.456e-6,4.139e-6,3.853e-6,3.596e-6,3.364e-6,3.153e-6,2.961e-6,2.785e-6,2.625e-6,2.479e-6,2.343e-6,2.219e-6,2.103e-6};
  double Aj3[39] = {1.026e-4,2.849e-5,1.360e-5,8.010e-6,5.287e-6,3.752e-6,2.799e-6,2.167e-6,1.726e-6,1.407e-6,1.168e-6,9.847e-7,8.410e-7,7.263e-7,6.334e-7,5.571e-7,4.936e-7,4.403e-7,3.950e-7,3.563e-7,3.230e-7,2.941e-7,2.689e-7,2.467e-7,2.270e-7,2.097e-7,1.943e-7,1.804e-7,1.680e-7,1.568e-7,1.466e-7,1.375e-7,1.291e-7,1.214e-7,1.145e-7,1.080e-7,1.022e-7,9.673e-8,9.169e-8};

  for(j=0;j<nj;j++)
  {
    tauj = 0;

    //eqn 21 of inoue et al. 2014
    if((lj[j]<lobs)&&(lobs<lj[j]*(1+z)))
    {
      if(lobs<2.2*lj[j]){
        tauj = Aj1[j]*pow(lobs/lj[j],1.2);
      }else if(2.2*lj[j]<=lobs && lobs<5.7*lj[j]){
        tauj = Aj2[j]*pow(lobs/lj[j],3.7);
      }else if(5.7*lj[j]<=lobs){
        tauj = Aj3[j]*pow(lobs/lj[j],5.5);
      }
    }
    tau += tauj;
  }

  //return combined tau
  return tau;
}
double IGM::tau_DLA_LS(long double lobs, float z) //checked
{
  int j;
  double tau  = 0;
  double tauj = 0;
  int nj = 39;
  //table 2 of inoue et al. 2014
  double lj[39]  = {1215.67, 1025.72, 972.537, 949.743, 937.803, 930.748, 926.226, 923.150, 920.963, 919.352, 918.129, 917.181, 916.429, 915.824, 915.329, 914.919, 914.576, 914.286, 914.039, 913.826, 913.641, 913.480, 913.339, 913.215, 913.104, 913.006, 912.918, 912.839, 912.768, 912.703, 912.645, 912.592, 912.543, 912.499, 912.458, 912.420, 912.385, 912.353, 912.324};
  double Aj1[39] = {1.617e-4,1.545e-4,1.498e-4,1.460e-4,1.429e-4,1.402e-4,1.377e-4,1.355e-4,1.335e-4,1.316e-4,1.298e-4,1.281e-4,1.265e-4,1.250e-4,1.236e-4,1.222e-4,1.209e-4,1.197e-4,1.185e-4,1.173e-4,1.162e-4,1.151e-4,1.140e-4,1.130e-4,1.120e-4,1.110e-4,1.101e-4,1.091e-4,1.082e-4,1.073e-4,1.065e-4,1.056e-4,1.048e-4,1.040e-4,1.032e-4,1.024e-4,1.017e-4,1.009e-4,1.002e-4};
  double Aj2[39] = {5.390e-5,5.151e-5,4.992e-5,4.868e-5,4.763e-5,4.672e-5,4.590e-5,4.516e-5,4.448e-5,4.385e-5,4.326e-5,4.271e-5,4.218e-5,4.168e-5,4.120e-5,4.075e-5,4.031e-5,3.989e-5,3.949e-5,3.910e-5,3.872e-5,3.836e-5,3.800e-5,3.766e-5,3.732e-5,3.700e-5,3.668e-5,3.637e-5,3.607e-5,3.578e-5,3.549e-5,3.521e-5,3.493e-5,3.466e-5,3.440e-5,3.414e-5,3.389e-5,3.364e-5,3.339e-5};
  for(j=0;j<nj;j++)
  {
    tauj = 0;

    //eqn 22 of inoue et al. 2014
    if((lj[j]<lobs)&&(lobs<lj[j]*(1+z)))
    {
      if(lobs<3.0*lj[j]){
        tauj = Aj1[j]*pow(lobs/lj[j],2.0);
      }else{
        tauj = Aj2[j]*pow(lobs/lj[j],3.0);
      }
    }
    tau += tauj;
  }

  //return combined tau
  return tau;
}
double IGM::tau_DLA_LC(long double lobs, float z) //checked
{
  int j;
  //eqn 28 and 29 of Inoue et al. 2014
  double tau = 0;
  double lL = 911.8; //Lyman limit
  if(z<2.0)
  {
    if(lobs<lL*(1+z))
    {
      tau = 0.211*pow(1+z,2) - 7.66e-2 * pow(1+z,2.3)*pow(lobs/lL,-0.3) - 0.135*pow(lobs/lL,2.);
    }
    else{
      tau = 0.;
    }
  }else if (z>=2.){ //POTENZIALE
    //z>=2.
    tau = 4.70e-2*pow(1+z,3)-1.78e-2*pow(1+z,3.3)*pow(lobs/lL,-0.3);
    if(lobs<3*lL)
    {
      tau += 0.634 -0.135*pow(lobs/lL,2.0) - 0.291*pow(lobs/lL,-0.3);
    }else if(lobs>=3.*lL && lobs<lL*(1+z)){
      tau += -2.92e-2*pow(lobs/lL,3);
    }else{
      tau = 0;
    }
  }
  return tau;
}
double IGM::tau_LAF_LC(long double lobs, float z) //checked
{
  //eqn 25, 26 and 27 of Inoue et al. 2014
  double tau = 0;
  double lL = 911.8; //Lyman limit
  if(z>=0 && z<1.2)
  {
    if(lobs<lL*(1+z))
    {
      tau = 0.325*(pow(lobs/lL,1.2)-pow(1+z,-0.9)*pow(lobs/lL,2.1));
    }
    else if(lobs>=lL*(1.+z)){
      tau = 0.;
    }
  }else if(z>=1.2 && z<4.7){
    //1.2<=z<4.7
    if(lobs<2.2*lL)
    {
      tau = 2.55e-2*pow(1+z,1.6)*pow(lobs/lL,2.1) + 0.325*pow(lobs/lL,1.2) - 0.250*pow(lobs/lL,2.1);
    }else if(lobs>=2.2*lL && lobs<lL*(1+z)){
      tau = 2.55e-2*(pow(1+z,1.6)*pow(lobs/lL,2.1) - pow(lobs/lL,3.7));
    }
    else if(lobs>=lL*(1.+z)){
      tau = 0.;
    }
  }else{
    //z>4.7
    if(lobs<2.2*lL)
    {
      tau = 5.22e-4*pow(1+z,3.4)*pow(lobs/lL,2.1) + 0.325*pow(lobs/lL,1.2) - 3.14e-2*pow(lobs/lL,2.1);
    }else if(lobs>=2.2*lL && lobs<5.7*lL){
      tau = 5.22e-4*pow(1+z,3.4)*pow(lobs/lL,2.1) + 0.218*pow(lobs/lL,2.1) - 2.55e-2*pow(lobs/lL,3.7);
    }else if(lobs>=5.7*lL && lobs<lL*(1+z)){
      tau = 5.22e-4*(pow(1+z,3.4)*pow(lobs/lL,2.1) - pow(lobs/lL,5.5));
    }
    else if(lobs>=lL*(1.+z)){
      tau = 0.;
    }    
  }
  return tau;
}

void IGM::igm_absorption(float zr, std::vector<long double> &wave, std::vector<long double> &sed)
{
  //lambda is in angstrom
  //sed is in erg/s/cm2/angstrom


  for (int l=0; l<wave.size(); l++){
    tsed.push_back(sed[l]*exp(-1.0*tau_IGM(wave[l],zr)));
    //cout << sed[l] << " " << wave[l] << " " << *exp(-1.0*tau_IGM(wave[l],zr)) << " " << tau_IGM(wave[l],zr) << endl; 
  }
  sed.clear();
  sed.shrink_to_fit();

  for (int l=0; l<wave.size(); l++){
    sed.push_back(tsed[l]);
  }

  tsed.clear();
  tsed.shrink_to_fit();  
}



/*! \fn double *igm_absorption(double *lambda, int n_lambda, double z);
 *  \brief Wavelength-dependent absorption owing to neutral H in the IGM.
 *
 *
 *   Provide an array of wavelengths in Angstrom, with the length of the 
 *   array n_lambda, at a redshift z and return an array of the attenuation 
 *   owing to neutral H in the IGM along the line of sight.
 *
 *   See Inoue et al. 2014
 */
/*
double IGM::*igm_absorption(double *lambda, int n_lambda, float z)
{
  //lambda is in angstroms
  //n_lambda is lambda.size()
  double *attenuation;
  attenuation  = (double *) calloc(n_lambda, sizeof(double));

  for(int i=0;i<n_lambda;i++)
    attenuation[i] = 1.0;

  double l;

  for(int j=0;j<n_lambda;j++)
  {
    l = lambda[j]; //leave in angstroms
    attenuation[j] *= exp(-1.0*tau_IGM(l,z));
  }

  return attenuation;
}
*/

/*
void print_igm_absorption(double *lambda, int n_lambda, double z)
{
  int i;
  FILE *fp;
  char fname[200];
  sprintf(fname,"igm_absorption.%3.1f.txt",z);
  double *absorption = igm_absorption(lambda, n_lambda, z);
  fp = fopen(fname,"w");
  for(i=0;i<n_lambda;i++)
    fprintf(fp,"%e\t%e\n",lambda[i],absorption[i]);
  fclose(fp);
  free(absorption);
}
*/
