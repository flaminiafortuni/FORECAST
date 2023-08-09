# FORECAST
FORECAST is a new flexible and adaptable software package that generates realistic astronomical images and galaxy surveys by forward modeling the output snapshot of any hydrodynamical cosmological simulation.
With customizable options for filters, size of the field of view, and survey parameters, it allows users to tailor the synthetic images to their specific research requirements.

### Table of Contents
1. [In few words](#in-few-words)
2. [Pipeline](#pipeline)
3. [Requirements](#requirements)
4. [Configuration file](#configuration-file)
5. [Input files format](#input-files-format)
6. [Output description](#output-description)
7. [How to use it?](#how-to-use-it)
8. [Contributing](#contributing)
9. [Citation and Acknowledgement](#citation-and-acknowledgment)


    
## In few words
FORECAST constructs a light cone centered on the observer's position exploiting the output snapshots of a chosen simulation and computes the observed fluxes of each simulated stellar element, modeled as a Single Stellar Population (SSP), in any chosen set of pass-band filters, including k-correction, IGM absorption, and dust attenuation. These fluxes are then used to create an image on a grid of pixels, to which observational features such as background noise and PSF blurring can be added with our Python script for the post-processing of the images (or by the user with their own script).


## Pipeline
The FORECAST code, currently available in *beta* version, is structured into four interconnected modules; each module relies on the output of the previous one. The first and third modules require the input file from the snapshots of hydrodynamical simulation to extract relevant properties of the simulated stellar and gas resolution elements.

The architecture of each module is not intrinsically parallel (e.g., it doesn't exploit MPI protocols), but it has been designed to allow the user to independently run it on multiple snapshots to realize multiple light-cone partitions (chunks of the cone) simultaneously.
Each module is dependent on the output of the previous module, and the first and the third module need the input file of the hydrodynamical simulations to extract the properties of the simulated stellar and gas resolution elements.
* [lightcone](lc) The first module handles the construction of the light cone, exploiting the data products of a chosen hydrodynamical simulation. In output, it produces an ASCII file containing the properties of the SSPs included in the field of view (IDs, coordinates, redshift, and physical properties), which is used as input file for the next module. This step might be skipped if a user already has their light cone, as long as the input file for the next module is written in the proper format. 
* [dust-free](df) The second module computes the flux of each SSP within the field of view, applying k-correction. It takes as input the catalog generated in the previous step. It generates an ASCII file containing the properties of the star particles and their dust-free fluxes in the selected filters.
* [dust-corrected](dc) The third module addresses the computation of the dust attenuation caused by dust distributed in and around each galaxy. This attenuation is traced by gas resolution elements included in the hydrodynamical simulation, providing the necessary data to calculate the dust attenuation effects (gas metallicity and neutral hydrogen column density). The module relies on data products from the hydrodynamical simulation to extract gas properties; it also takes as input the catalog generated in the previous step. As output, it generates an ASCII file containing the properties of the star particles and their dust-corrected fluxes in the selected filters.
* [IGM](igm) The final module adds the IGM correction to dust-corrected fluxes, taking as input the catalog generated in the previous step and producing the final output catalog, which includes the physical properties of the stellar particles and their corrected fluxes.
### Post-processing scripts
An independent C++ script handles the arrangement of the fluxes on a grid of pixels with Npix-per-side chosen by the user.   
Two additional independent scripts, written in Python, are made available (i) to build the galaxy catalog, in ASCII format, from the particle catalog given in output by FORECAST; (ii) to post-process the FORECAST images with our noise and PSF pipeline.

## Requirements
FORECAST is written in C and C++; it is supported by independent libraries to make the code more readable and user-friendly. 
It requires the following C/C++ standard libraries: 
* [GSL](https://www.gnu.org/software/gsl/)
* [openBLAS](https://www.openblas.net/)
* [LAPACK](http://www.netlib.org/lapack/)
* [CCfits](https://heasarc.gsfc.nasa.gov/fitsio/CCfits/)
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
* [FFTW](https://www.fftw.org/)
* [Eigen](https://eigen.tuxfamily.org/)
* [Armadillo](http://arma.sourceforge.net/)
* [H5Cpp](http://h5cpp.org/)
* [HDF5 C++](https://www.hdfgroup.org/)
* gcc compiler.      

## Configuration file
FORECAST is adaptable to the choices of the user by selecting a set of input parameters. 
The input parameters can be changed in the *.ini file(s).
In particular, it is possible to choose the hydrodynamical simulation that provides the backbone of the light cone (box with side-length L_box); the highest redshift to be included, z_s, which determines the maximum distance covered by the light cone, D_s; the dimensions of the field of view, L_{fov}; the resolution of the ideal simulated images, setting the number of pixels per image side-length N_{pix}; the filters, by the order they are written in the filters.dat file.


## Input files format
FORECAST works by taking as input the physical properties of stellar particles and gas cells recorded in each snapshot of the chosen simulation. Currently, 

For stellar particles (" _PartType4_ ") 
- " _Coordinates_ " x,y,z-comoving coordinates within the simulated volume, in ckpc
- " _Masses_ " stellar mass, in $M_{\odot}$
- " _InitialMass_ " initial stellar mass, in $M_{\odot}$
- " _GFM_Metallicity_ " metallicity as $M_{Z}/M_{TOT}$
- " _GFM_StellarFormationTime_ " scale factor indicating the exact time when a star is formed, to compute age in yr
- " particle subhalo membership ID " this info is not available; it is reconstructed with a designed algorithm that uses _SubFind_ subhalo catalogs and _Offset_ files, exploiting the specific organization of halos and subhalos within the catalog files.

For gas cells (" _PartType0_ "):
- " _Coordinates_ " x,y,z-comoving coordinates of the geometrical center within the snapshot, in ckpc
- " _Masses_ " gas mass, in $M_{\odot}$; needed to compute Volume
- " _Density_ " gas comoving volume, in $ckpc^{3}$; needed to compute Volume
- " _GFM_Metallicity_ " gas metallicity, in $M_{Z}/M_{TOT}$ (with $M_Z$ the total mass all metal elements)
- " _ElectronAbundance_ " gas cell fractional electron number density in $ckpc^{-2}$ with respect total hydrogen number density, 
- " _InternalEnergy_ " gas cell internal (thermal) energy per unit mass $u$, needed to derive the neutral hydrogen column density
- " cell subhalo membership ID " same as for stellar particles

For the exact descriptions of the fields, please refer to * [IllustrisTNG/DataSpecification]([https://heasarc.gsfc.nasa.gov/fitsio/](https://www.tng-project.org/data/docs/specifications/#parttype0)).

## Output description
The output of the code is the catalog including the physical properties and the true fluxes of the simulated stellar particles. It is used to build the galaxy catalog.
    The catalogs, both the particles and the galaxy ones, have different sizes depending on the number of particles (or galaxies) included in the field of view and typically grow in size as the redshift increases since more structures are included. The total size of output files is 1 TB.
    The output images are recorded on 16-bit floating-point FITS files. Each plane (projection on a bi-dimensional map of fluxes from a volume of the Universe included in the field of view, in a redshift range) occupies 3 GB, while the size of the final stacked image is 5 GB.
    
## How to use it?
- Download the code github.com/FORECAST
- Download the required C/C++ libraries
- Download the hydrodynamical simulation input files
- Change paths and parameters in the input configuration file \*.ini (in each module)
- Compile Makefile
- run bash \*.sh file (one partition per time; refer to planes_list.txt file to choose the partition/snapshot-plane run) 

## Contributing

## Citation and Acknowledgment

When using the code, please cite the release paper [Fortuni et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230519166F/abstract). The paper has been submitted to A&A; we will update this section.

```@ARTICLE{2023arXiv230519166F,
       author = {{Fortuni}, Flaminia and {Merlin}, Emiliano and {Fontana}, Adriano and {Giocoli}, Carlo and {Romelli}, Erik and {Graziani}, Luca and {Santini}, Paola and {Castellano}, Marco and {Charlot}, St{\'e}phane and {Chevallard}, Jacopo},
        title = "{FORECAST: a flexible software to forward model cosmological hydrodynamical simulations mimicking real observations}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Astrophysics of Galaxies, 85},
         year = 2023,
        month = may,
          eid = {arXiv:2305.19166},
        pages = {arXiv:2305.19166},
          doi = {10.48550/arXiv.2305.19166},
archivePrefix = {arXiv},
       eprint = {2305.19166},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv230519166F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

