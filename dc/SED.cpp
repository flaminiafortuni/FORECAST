#include <cmath>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "SED.h"
//#include "IGM.h"
using namespace arma;
using namespace std;
//IGM igmy;


void SED::extract_spec(std::vector <vector<double> > &full_table, std::vector<double> &time_grid, int a_indx, float ages4, std::vector<long double> &spe){
  
//old_UNUSED                                                                                                                                     
  if(ages4<=0.){
    for (int i=0; i<13391; i++){
      spe.push_back(full_table[0][i]*ergsa);
    }
  } else if(ages4>=20.){
    for (int i=0; i<13391; i++){
      spe.push_back(full_table[220][i]*ergsa);
    }
  }
  else if (ages4>0. & ages4<20.){
    for (int i=0; i<13391; i++){
      spe.push_back(full_table[a_indx][i]*ergsa);
    }
  }
}




double SED::getY(std:: vector<double> x, std:: vector<double> y,double xi){
  nn = x.size();
  if(x[0]<x[nn-1]){         
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }  
  else{
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];  
  }
  i= locate (x,xi);
  i = std::min (std::max (i,0), int (nn)-2);
  f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    f2 = f*f;
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];
    a1 = y[i-1] - y[i] - a0;
    a2 = y[i+1] - y[i-1];
    a3 = y[i];
    return a0*f*f2+a1*f2+a2*f+a3;
  }                                                                      
  else return f*y[i+1]+(1-f)*y[i];           
}

template <class T> int SED::locate(const std::vector<T> &v, const T x){
  n = v.size ();
  jl = -1;
  ju = n;
  as = (v[n-1] >= v[0]);
  while (ju-jl > 1){
    jm = (ju+jl)/2;
    if ((x >= v[jm]) == as)
      jl=jm;
    else
      ju=jm;
  }
  if (x == v[0])
    return 0;
  else if (x == v[n-1])
    return n-2;
  else
    return jl;
}


////F(lambda) original algorithm

void SED::z_evol(float zr, std::vector<long double> &wave, std::vector<long double> &sed, std::vector<double> &zl, std::vector<double> &dlum ){ 

  //z_evol for observed-frame apparent magnitude (i.e. m with d=DL and z>0)
  //it becomes in rest-frame apparent magnitude when z=0
  
  dls=0.;
  dlsmpc=0;
   
  if (zr>0.){
    dlsmpc=SED::getY(zl,dlum,zr); //get Dl knowing SSP's z 
    dls=dlsmpc*mpctocm/0.6774; //Dl in cm 
  }
  else if (zr==0.){
    dls= 3.086e+19; //10 pc to cm
  }
 
 
  //redshifting and computing F(lambda/AA) from the input L(lambda/AA)
  for (int l=0; l<wave.size(); l++){
    twave.push_back(wave[l]*(1+zr));
    tsed.push_back(sed[l]/((1+zr)*pow(dls,2)*4*M_PI));
  }
  wave.clear();
  sed.clear();
  
  //put again in original vectors
  for (int l=0; l<twave.size(); l++){
    wave.push_back(twave[l]);
    sed.push_back(tsed[l]);
  }

  tsed.clear();
  twave.clear();

}



void SED::z_evolABS(float zr, std::vector<long double> &wave, std::vector<long double> &sed ){ 

  //z_evol for observed-frame absolute magnitude (i.e. M with z>0 (obs-frame) but dL=10 pc (abs mag) )
  //it becomes for rest-frame absolute magnitude when z=0
  
  dls= 3.086e+19; //10 pc to cm
  
 
  //redshifting and computing F(lambda/AA) from the input L(lambda/AA)
  for (int l=0; l<wave.size(); l++){
    twave.push_back(wave[l]*(1+zr));
    tsed.push_back(sed[l]/((1+zr)*pow(dls,2)*4*M_PI));
  }
 
  
  //rebin the spectrum in 1e6 points evenly distributed [14,2.000014e6] AA

  wave_i.resize(1000000);
  for (x = wave_i.begin(), val = 14.; x != wave_i.end(); ++x, val += h_step) {
    *x = val;
  }
 
   
  wave_interp  =  arma::conv_to< vec >::from(wave_i);
  nsed  =  arma::conv_to< vec >::from(tsed);
  nwave =  arma::conv_to< vec >::from(twave);

  
  arma::interp1(nwave,nsed,wave_interp,sed_interp,"*linear",0);

  tsed.clear();
  twave.clear();
  
  nsed.clear();
  nwave.clear();
  
  //back to std::vec to modify the vectors
  tsed  =  arma::conv_to< stdvec >::from(sed_interp);
  twave =  arma::conv_to< stdvec >::from(wave_interp);
  
  
  wave_i.clear();
  sed_interp.clear();
  wave_interp.clear();
  

  //for(int l=0;l<twave.size();l++)
  //tespo.push_back(igmy.tau_i14(zr,twave[l])); //Inoue+14
  
  wave.clear();
  sed.clear();

  for (int l=0; l<twave.size(); l++){
    wave.push_back(twave[l]);
    sed.push_back(tsed[l]);//*tespo[l]);
  }

  
  tespo.clear();
  twave.clear();
  tsed.clear();


}



double SED::compute_mab(float zr, std::vector <long double> &waves,std::vector <long double> &sed, std::vector <long double> &fwaves,std::vector <long double> &fresp){

  //GALAXEV 2003 paper; Blanton et al. 2003 -> m_AB integral with F(lambda)

  //find lambda_min and lambda_max of the filter with a treshold (1e-3)
  for (auto i=0;i<fresp.size();i++){
    if (fresp[i]>=1.e-4){
      lminf=fwaves[i];
      in_Lmin=i;
      break;
    }
  }
  
  for (auto i=fresp.size();i--!=0;){
    if (fresp[i]>=1.e-4){
      lmaxf=fwaves[i];
      in_Lmax=i;
      break;
    }
  }

  
  //filter response and wavelength (AA) in arma::vec
  cfresp  =  arma::conv_to< vec >::from(fresp);
  cfwaves =  arma::conv_to< vec >::from(fwaves);


  //select F(lambda/AA) ONLY in the range of wavelength of the filter response
  for (int i=0; i<sed.size(); i++){
    if ((waves[i]>=lminf) &&  (waves[i]<=lmaxf)){
      newwaves.push_back(waves.at(i));
      newsed.push_back(sed.at(i));
    }
  }


  ////NEW REBIN: rebin flux density (newwaves,newsed) between lminf and lmaxf with a step typically of 1 AA
 
  //resize wave_i by size_wave 
  int size_wave=(int)(2*(lmaxf-lminf))+1;
  //find h_step
  h_step=(lmaxf-lminf)/(size_wave);
  
  //rebin the spectrum in size_wave (lambda [AA]) interval
  wave_i.resize(size_wave);
  for (x = wave_i.begin(), val = lminf; x != wave_i.end(); ++x, val += h_step) {
    *x = val;
  }
   
  wave_interp  =  arma::conv_to< vec >::from(wave_i);
  nsed  =  arma::conv_to< vec >::from(newsed);
  nwave =  arma::conv_to< vec >::from(newwaves);

  
  arma::interp1(nwave,nsed,wave_interp,sed_interp,"*linear",0);

  newsed.clear();
  newwaves.clear();
  
  nsed.clear();
  nwave.clear();
  ////END NEW REBIN 

 
  //back to std::vec to modify the vectors
  //F(lambda/AA) and wavelength in erg/s/cm2/AA and in AA. Possibly z-evolved
  csed    =  arma::conv_to< stdvecd >::from(sed_interp);
  cwaves  =  arma::conv_to< stdvecd >::from(wave_interp);
  
  sed_interp.clear();
  wave_interp.clear();

  //Interpolating the filter response in the same wavelengths as F(lambda/AA)
  arma::interp1(cfwaves,cfresp,cwaves,filter_interp,"*linear",cfresp[0]);
  
  
  //integrand kernel of the numerator of the m_AB integral
  filterSpec = filter_interp % csed; //"%" in arma:: is equivalent to "*"
  i1 = filterSpec % cwaves;
  
  //integrand kernel of the denominator of the m_AB integral
  i2 = filter_interp % (1./cwaves);
  
  //numerator and denominator integrals of m_AB
  I1 = arma::trapz(cwaves,i1);
  I2 = arma::trapz(cwaves,i2);
 
  flambda = arma::as_scalar(I1/I2);
  fnu= flambda/speedcunitas;

  mAB = -2.5*log10(fnu)-48.6;
  /*
  cout << "I1 " <<  I1 << endl;
  cout << "I2 " <<  I2 << endl;
  cout << "fl " <<  flambda << endl;
  cout << "fn " <<  fnu << endl;
  cout << "mab" <<  mAB << endl;
  cout << "c" <<  speedcunitas << endl;
  */

  
  flambda=0.;
  fnu=0.;
  cfresp.clear();
  cfwaves.clear();
  csed.clear();
  cwaves.clear();
  filter_interp.clear();
  filterSpec.clear();
  i1.clear();
  i2.clear();
  I1.clear();
  I2.clear();

  return mAB;
}


void SED::dust_attenuation( std::vector <float> &ex_l,std::vector <long double> &dustAtt, std::vector <long double> &wavesc,std::vector <long double> &spe){
 
  //rebin the dust attenuation in 1221 waves points
  ex_l_interp  =  arma::conv_to< vec >::from(wavesc);
  a_spe  =  arma::conv_to< vec >::from(spe);
  
  a_ex_l =  arma::conv_to< vec >::from(ex_l);
  a_dustAtt  =  arma::conv_to< vec >::from(dustAtt);
 
  arma::interp1(a_ex_l,a_dustAtt,ex_l_interp,dustAtt_interp,"*linear",1);

  dustAtt.clear();
  ex_l.clear();
  dustAtt.shrink_to_fit();
  ex_l.shrink_to_fit();

  spe.clear();
  spe.shrink_to_fit();
  
 
  a_dustAtt.clear();
  a_ex_l.clear();

  out_spe= dustAtt_interp % a_spe;

  //back to std::vec to modify the vectors -> 1221 component as spe
  spe  =  arma::conv_to< stdvec >::from(out_spe);

  
  a_spe.clear();
  out_spe.clear();

  dustAtt=  arma::conv_to< stdvec >::from(dustAtt_interp);
  ex_l =  arma::conv_to< stdvec2 >::from(ex_l_interp);

  
  dustAtt_interp.clear();
  ex_l_interp.clear();


}


void SED::modelC( float zr, std::vector <long double> &ex_curve,std::vector <long double> &ex_l,  std::vector<float> &ex_ls, std::vector <long double> &modC, double NHImean, double ZHImean){

//model C Nelson+19

  //std::vector<float> ex_ls(15);
  std::copy ( ex_l.begin(), ex_l.end(), ex_ls.begin() );
  for ( auto &item : ex_ls ) item *= 1.+zr;

  // taua= pow((1.+zr),-0.5)*NHImean; //this was the previous z-scaling, consistnt with Fortuni+23
  taua= pow(1+zr,-1.95)*NHImean;

  
  //albedo: Omega_lambda in Calzetti+94
  for (auto k=0;k<ex_curve.size();k++){
    double y=log10(ex_ls[k]);
    if((ex_ls[k]<=3460.) && (ex_ls[k]>=1000.)){
      omegal.push_back(0.43+0.366*(1.-exp(-(pow((y-3.),2.)/0.2))));
    }
    else if((ex_ls[k]<=7000.) && (ex_ls[k]>3460.)){
      omegal.push_back(-0.48*y+2.41);
    }
    else{
      omegal.push_back(0.);
    }
  }
  
  
  //weight parameter describing anisotropy scattering: h_lambda in Calzetti+94
  for (auto k=0;k<ex_curve.size();k++){
    double y=log10(ex_ls[k]);	
    if((ex_ls[k]<=7000.) && (ex_ls[k]>=1200.)){
      hl.push_back(1.-0.561*exp(-(pow(abs(y-3.3112),2.2)/0.17)));
    }
    else{
      hl.push_back(0.); // combined with omegal for l>7000 gives tau_scatt=1
    }
  }
  
  
  //final tau in Nelson+19
  for (auto k=0;k<ex_curve.size();k++){  
    if(ex_ls[k]<2000.){
      long double tau_abs=taua*pow((ZHImean),1.35);
      long double tau_s= hl[k]*pow((1.-omegal[k]),0.5)+(1-omegal[k])*(1-hl[k]);
      long double tauT = tau_abs*tau_s;
      modC.push_back(1./(tauT)*(1.-exp(-tauT))); //in the Fortuni+23 paper version, ex_curve[k] multiplied tauT	  
    }
    else if(ex_ls[k]>2000.){
      long double tau_abs=taua*pow((ZHImean),1.6);
      long double tau_s= hl[k]*pow((1.-omegal[k]),0.5)+(1-omegal[k])*(1-hl[k]);
      long double tauT = tau_abs*tau_s;
      modC.push_back(1./(tauT)*(1.-exp(-tauT))); //in the Fortuni+23 paper version, ex_curve[k] multiplied tauT	  	
    }
  }

  omegal.clear();
  hl.clear();
  omegal.shrink_to_fit();
  hl.shrink_to_fit();
  
}
