/*
 * @file src/ReadTNGPartcle.cpp
 * @date 25/09/20
 * @author Erik Romelli - INAF-OATs
 */

#include <vector>
#include <iostream>
#include "/usr/include/hdf5/serial/H5Cpp.h"
#include </usr/include/eigen3/Eigen/Dense>
#include "readTNGParticle.h"
using namespace std;



bool readTNGParticle::pathExists(hid_t id, const std::string& path)
{
  
  return H5Lexists( id, path.c_str(), H5P_DEFAULT ) > 0;
  
}


void readTNGParticle::Initialize(std::string workdir, int particleID){

  snapNameTemplate = workdir+"/snapdir_0"+std::to_string(particleID)+"/snap_0"+std::to_string(particleID);
  fofNameTemplate = workdir+"/groups_0"+std::to_string(particleID)+"/fof_subhalo_tab_0"+std::to_string(particleID);
  offsetFileName = workdir+"/offsets_0"+std::to_string(particleID)+".hdf5";

}


template <class T>  T readTNGParticle::getHeaderValue(std::string keyword){

  T tmp;
  H5::Attribute attr = snapHeader.openAttribute(keyword);
  H5::DataType attrType = attr.getDataType();
  attr.read(attrType, &tmp);
  return tmp;

}


template <class T>  std::vector<T> readTNGParticle::getHeaderValueNonScalar(std::string keyword){

  T tmpArray[6];
  H5::Attribute attr = snapHeader.openAttribute(keyword);
  H5::DataType attrType = attr.getDataType();
  attr.read(attrType, &tmpArray);

  // Trick to convert array to std::vector
  int n = sizeof(tmpArray)/sizeof(tmpArray[0]);
	std::vector<T> tmp(tmpArray, tmpArray+n);

  return tmp;

}


std::vector<double> readTNGParticle::getSnapVectorValue(std::string keyword){

  H5::DataSet dataset = snapStars.openDataSet(keyword);
  H5::DataSpace dataspace = dataset.getSpace();

  hsize_t naxes[1];

  dataspace.getSimpleExtentDims(naxes, NULL);

  std::vector<double> tmp(naxes[0], 0.);
  dataset.read(tmp.data(), H5::PredType::NATIVE_DOUBLE);

  return tmp;

}


void readTNGParticle::readCoordinates(){

  H5::DataSet dataset = snapStars.openDataSet("Coordinates");
  H5::DataSpace dataspace = dataset.getSpace();
  hsize_t naxes[2];
  dataspace.getSimpleExtentDims(naxes, NULL);

  coordinates = Eigen::MatrixXd::Zero(naxes[1], naxes[0]);

  dataset.read(coordinates.data(), H5::PredType::NATIVE_DOUBLE);

  x = std::vector<double>(naxes[0], 0.);
  y = std::vector<double>(naxes[0], 0.);
  z = std::vector<double>(naxes[0], 0.);

  for (int i = 0; i < naxes[0]; i++){
    x[i] = coordinates(0,i);
    y[i] = coordinates(1,i);
    z[i] = coordinates(2,i);
  }

}


void readTNGParticle::readHeader(int cutID){

  std::string snapFileName = snapNameTemplate + "." + std::to_string(cutID) + ".hdf5";
  snap = H5::H5File(snapFileName, H5F_ACC_RDONLY );

  snapHeader = snap.openGroup("Header");

  boxSize = readTNGParticle::getHeaderValue<double>("BoxSize");
  omegaZero = readTNGParticle::getHeaderValue<double>("Omega0");
  omegaLambda = readTNGParticle::getHeaderValue<double>("OmegaLambda");
  time = readTNGParticle::getHeaderValue<double>("Time");
  redshift = readTNGParticle::getHeaderValue<double>("Redshift");
  massTable = readTNGParticle::getHeaderValueNonScalar<double>("MassTable");
  numPartTotal = readTNGParticle::getHeaderValueNonScalar<int>("NumPart_Total");
  numPartTotal_HW = readTNGParticle::getHeaderValueNonScalar<int>("NumPart_Total_HighWord");

  snapHeader.close();

}


void readTNGParticle::readStars(int cutID){

  std::string snapFileName = snapNameTemplate + "." + std::to_string(cutID) + ".hdf5";
  snap = H5::H5File(snapFileName, H5F_ACC_RDONLY );

  // Open Header group
  snapStars = snap.openGroup("PartType4");

  gmfMetallicity = readTNGParticle::getSnapVectorValue("GFM_Metallicity");
  gmfInitialMass = readTNGParticle::getSnapVectorValue("GFM_InitialMass");
  gmfStellarFormationTime = readTNGParticle::getSnapVectorValue("GFM_StellarFormationTime");
  masses = readTNGParticle::getSnapVectorValue("Masses");
  readTNGParticle::readCoordinates();

  snapStars.close();

}


void readTNGParticle::readFof(int cutID){

  std::string fofFileName = fofNameTemplate + "." + std::to_string(cutID) + ".hdf5";
  fof = H5::H5File(fofFileName, H5F_ACC_RDONLY );

  H5::Group subhalo = fof.openGroup("Subhalo");

  try {
    H5::Exception::dontPrint();

    H5::DataSet dataset = subhalo.openDataSet("SubhaloLenType");
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t naxes[2];
    dataspace.getSimpleExtentDims(naxes, NULL);

    subhaloLenType = Eigen::MatrixXd::Zero(naxes[1], naxes[0]);
    
    dataset.read(subhaloLenType.data(), H5::PredType::NATIVE_DOUBLE);
    
    starsLenType = std::vector<int>(naxes[0], 0.);
    gasLenType = std::vector<int>(naxes[0], 0.);
    
    
    for (int i = 0; i < naxes[0]; i++){
      starsLenType[i] = subhaloLenType(4,i);
      gasLenType[i] = subhaloLenType(0,i);      
    }
  }
  catch( H5::GroupIException error ){
    starsLenType = std::vector<int>(0);
    gasLenType = std::vector<int>(0);
  }
 

  fof.close();

}


void readTNGParticle::readOffset(){

  offset = H5::H5File(offsetFileName, H5F_ACC_RDONLY );

  H5::Group subhalo = offset.openGroup("Subhalo");

  try {

    H5::Exception::dontPrint();
    
    H5::DataSet dataset = subhalo.openDataSet("SnapByType");
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t naxes[2];
    dataspace.getSimpleExtentDims(naxes, NULL);

    snapByType = Eigen::MatrixXd::Zero(naxes[1], naxes[0]);
    
    dataset.read(snapByType.data(), H5::PredType::NATIVE_DOUBLE);
    
    starsByType = std::vector<int>(naxes[0], 0.);
    gasByType = std::vector<int>(naxes[0], 0.);
    

    for (int i = 0; i < naxes[0]; i++){
      starsByType[i] = snapByType(4,i);
      gasByType[i] = snapByType(0,i);
    }
  }
  catch( H5::GroupIException error ){
    ;
  }

  offset.close();

} 



/*
* GAS functions
*/

std::vector<double> readTNGParticle::getGASSnapVectorValue(std::string keyword){

  H5::DataSet dataset = snapGAS.openDataSet(keyword);
  H5::DataSpace dataspace = dataset.getSpace();

  hsize_t naxes[1];

  dataspace.getSimpleExtentDims(naxes, NULL);

  std::vector<double> tmp(naxes[0], 0.);
  dataset.read(tmp.data(), H5::PredType::NATIVE_DOUBLE);

  return tmp;

}


void readTNGParticle::readGASCoordinates(){

  H5::DataSet dataset = snapGAS.openDataSet("CenterOfMass");
  H5::DataSpace dataspace = dataset.getSpace();
  hsize_t naxes[2];
  dataspace.getSimpleExtentDims(naxes, NULL);

  coordinates_gas = Eigen::MatrixXd::Zero(naxes[1], naxes[0]);

  dataset.read(coordinates_gas.data(), H5::PredType::NATIVE_DOUBLE);

  x_gas = std::vector<double>(naxes[0], 0.);
  y_gas = std::vector<double>(naxes[0], 0.);
  z_gas = std::vector<double>(naxes[0], 0.);

  for (int i = 0; i < naxes[0]; i++){
    x_gas[i] = coordinates_gas(0,i);
    y_gas[i] = coordinates_gas(1,i);
    z_gas[i] = coordinates_gas(2,i);
  }

} 


void readTNGParticle::readGASElements(){

  H5::DataSet dataset = snapGAS.openDataSet("GFM_Metals");
  H5::DataSpace dataspace = dataset.getSpace();
  hsize_t naxes[2];
  dataspace.getSimpleExtentDims(naxes, NULL);

  elements_gas = Eigen::MatrixXd::Zero(naxes[1], naxes[0]); 

  dataset.read(elements_gas.data(), H5::PredType::NATIVE_DOUBLE);

  H_gas = std::vector<double>(naxes[0], 0.);
  He_gas = std::vector<double>(naxes[0], 0.);
  C_gas = std::vector<double>(naxes[0], 0.);

  for (int i = 0; i < naxes[0]; i++){
    H_gas[i] = elements_gas(0,i);
    He_gas[i] = elements_gas(1,i);
    C_gas[i] = elements_gas(2,i);
  }

} 


void readTNGParticle::readGAS(int cutID){

  std::string snapFileName = snapNameTemplate + "." + std::to_string(cutID) + ".hdf5";
  snap = H5::H5File(snapFileName, H5F_ACC_RDONLY );

  // Open Header group
  snapGAS = snap.openGroup("PartType0"); //gas

  gmfMetallicity_gas = readTNGParticle::getGASSnapVectorValue("GFM_Metallicity");
  HIAbundance_gas = readTNGParticle::getGASSnapVectorValue("NeutralHydrogenAbundance");
  masses_gas = readTNGParticle::getGASSnapVectorValue("Masses");
  density_gas = readTNGParticle::getGASSnapVectorValue("Density");
  InternalEnergy_gas = readTNGParticle::getGASSnapVectorValue("InternalEnergy");
  eAbundance_gas = readTNGParticle::getGASSnapVectorValue("ElectronAbundance");
  readTNGParticle::readGASCoordinates();
  readTNGParticle::readGASElements();

  snapGAS.close();

} 


/*
* get functions for Header info
*/

double readTNGParticle::getBoxSize(){

  return boxSize;

} //getBoxSize END

double readTNGParticle::getOmegaZero(){

  return omegaZero;

} //getOmegaZero END

double readTNGParticle::getOmegaLambda(){

  return omegaLambda;

} //getOmegaLambda END

double readTNGParticle::getTime(){

  return time;

} //getTime END

double readTNGParticle::getRedshift(){

  return redshift;

} //getRedshift END

std::vector<double> readTNGParticle::getMassTable(){

  return massTable;

} //getMassTable END

std::vector<int> readTNGParticle::getNumPartTotal(){

  std::vector<int> n(numPartTotal.size(), 0.);

  for (int i = 0; i < numPartTotal.size(); i++){
    n[i] = numPartTotal[i] | (numPartTotal_HW[i] << 32);
  }

  return n;

} //getNumPartTotal END

/*
* get functions for Stars info
*/

std::vector<double> readTNGParticle::getMetallicity(){

  return gmfMetallicity;

} //getMetallicity END

std::vector<double> readTNGParticle::getInitialMass(){

  return gmfInitialMass;

} //getInitialMass END

std::vector<double> readTNGParticle::getStellarFormationTime(){

  return gmfStellarFormationTime;

} //getStellarFormationTime END

std::vector<double> readTNGParticle::getMasses(){

  return masses;

} //getMasses END

std::vector<double> readTNGParticle::getX(){

  return x;

} //getX END

std::vector<double> readTNGParticle::getY(){

  return y;

} //getX END

std::vector<double> readTNGParticle::getZ(){

  return z;

} //getX END

/*
* get functions for fof info
*/

std::vector<int> readTNGParticle::getStarsLenType(){

  return starsLenType;

} //getStarsLenType END

std::vector<int> readTNGParticle::getGasLenType(){

  return gasLenType;

}

/*
* get functions for offset info
*/

std::vector<int> readTNGParticle::getStarsByType(){

  return starsByType;

} //getStarsByType END

std::vector<int> readTNGParticle::getGasByType(){

  return gasByType;

}


/*
* get functions for GAS info
*/

std::vector<double> readTNGParticle::getGASMetallicity(){

  return gmfMetallicity_gas;

} //getMetallicity END

std::vector<double> readTNGParticle::getGASMasses(){

  return masses_gas;

} //getMasses END

std::vector<double> readTNGParticle::getGASDensity(){

  return density_gas;

} //getDensity END

std::vector<double> readTNGParticle::getGASX(){

  return x_gas;

} //getX END

std::vector<double> readTNGParticle::getGASY(){

  return y_gas;

} //getX END

std::vector<double> readTNGParticle::getGASZ(){

  return z_gas;

} //getX END

std::vector<double> readTNGParticle::getGASHIAbundance(){

  return HIAbundance_gas;

} //getHIAbundance END 

std::vector<double> readTNGParticle::getGASH(){

  return H_gas;

} //getGASH END

std::vector<double> readTNGParticle::getGASHe(){

  return He_gas;

} //getGASHe END

std::vector<double> readTNGParticle::getGASC(){

  return C_gas;

  } //getGASC END

std::vector<double> readTNGParticle::getGASIntEnergy(){

  return InternalEnergy_gas;

} //get internal energy END

std::vector<double> readTNGParticle::getGASeAbundance(){

  return eAbundance_gas;

} //get electron abundance END
