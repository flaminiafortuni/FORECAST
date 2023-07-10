# FORECAST
*The text is extracted from Fortuni et al. 2023*

FORECAST is a new flexible and adaptable software package that performs forward modeling of the output of any cosmological hydrodynamical simulations to create a wide range of realistic synthetic astronomical images. With customizable options for filters, size of the field of view and survey parameters, it allows users to tailor the synthetic images to their specific requirements.

    
**In few words**
FORECAST constructs a light-cone centered on the observer's position exploiting the output snapshots of a simulation and computes the observed flux of each simulated stellar element, modeled as a Single Stellar Population, in any chosen set of pass-band filters, including k-correction, IGM absorption and dust attenuation. These fluxes are then used to create an image on a grid of pixels, to which observational features such as background noise and PSF blurring can be added. This allows to obtain simulated galaxies with realistic morphologies and star formation histories.


**Pipeline**
The FORECAST code is built in four modules: lc, df, dc, igm (stand for: lightcone construction module; dust-free module; dust-corrected module; igm module). The architecture of each module is not intrinsically parallel (e.g., it doesn't exploit MPI protocols), but it has been designed to allow the user to independently run it on multiple snapshots to realize multiple light-cone partitions (chunks of the cone) simultaneously.
    Each module is dependent from the output of the previous module, and the first and the third module need the input file of the hydrodynamical simulations to extract the properties of the simulated stellar and gas resolution elements.
    (lc) The first module handles the construction of the light-cone, exploiting the data products of a chosen hydrodynamical simulation. In output, it produces an ASCII file containing the properties of the SSPs included in the FoV (IDs, coordinates, redshift and physical properties), which is used as input file for the next module. This step might be skipped if a user already has their light-cone, as long as the input file for the next module is written in the proper format. 
    (df) The second module computes the dust-free flux of each SSPs included in the FoV. It assembles an ASCII file with the properties of the star particles and their dust-free fluxes in chosen filters.
    (dc) The third module addresses the computation of dust-corrected fluxes and it requires the data products of the hydrodynamical simulation to extract the properties of gas resolution elements belonging to the sources included in the FoV. In output it is given the same file produced by the previous module, also including dust-corrected fluxes for stellar particles and the gas mass-weighted mean of the gas metallicity and the neutral hydrogen column density. 
    (igm) The final module adds the IGM correction to dust-corrected fluxes, producing the final output catalogue, which includes the physical properties of the stellar particles and their corrected fluxes, and the mean properties of the gas   
    An independent C++ script handles the arrangement of the fluxes on a grid of pixels with Npix-per-side chosen by the user.   
    Two additional independent scripts, written in python, are made available (i) to build the galaxy catalogue, in ASCII format, from the particle catalogue given in output by FORECAST; (ii) to post-process the FORECAST images with our noise and PSF pipeline.

**Input file format**
----

**Configuration file .ini**
FORECAST is adaptable to the choices of the user by selecting a set of input parameters. 
In particular, it is possible to choose the hydrodynamical simulation that provides the backbone of the light-cone (box with side-length L_box); the highest redshift to be included, z_s, which determines the maximum distance covered by the light-cone, D_s; the dimensions of the field of view, L_{fov}; the resolution of the ideal simulated images, setting the number of pixels per image side-length N_{pix}; the fiters, by the order they are written in the filters.dat file.


**Ouput**
The output of the code is the catalogue including the physical properties and the true fluxes of the simulated stellar particles. It is used to build the galaxy catalogue.
    The catalogues, both the particles and the galaxy ones, have different sizes depending on the number of particles (or galaxies) included in the FoV, and typically grow in size as the redshift increases since more structures are included. The total size of output files is 1 TB.
    The output images are recorded on 16-bit floating-point FITS files. Each plane (projection on a bi-dimensional map of fluxes from a volume of the Universe included in the field of view, in a redshift range) occupies 3 GB, while the size of the final stacked image is 5 GB.
    
 
**Requirements**
FORECAST is written in C and C++; it is supported by independent libraries to make the code more readable and user-friendly. 
It requires the following C/C++ standard libraries: GSL {https://www.gnu.org/software/gsl/}{gsl},  openBLAS {https://www.openblas.net/}, LAPACK  {http://www.netlib.org/lapack/}, CCfits {https://heasarc.gsfc.nasa.gov/fitsio/CCfits/},  CFITSIO {https://heasarc.gsfc.nasa.gov/fitsio/}, FFTW {https://www.fftw.org/},  Eigen {https://eigen.tuxfamily.org/},  Armadillo {http://arma.sourceforge.net/},  H5Cpp {http://h5cpp.org/}, HDF5 C++ {https://www.hdfgroup.org/} and a gcc compiler.   
In the input configuration file of the code, the user chooses the image simulation parameters (e.g., the dimension of the FoV, the filters, the hydrodynamical simulation.  
   The code requires the input files of the chosen hydrodynamical simulation to build the light-cone and produce the final images. The data products of the numerous currently available simulations are organized differently, changing from one to another simulation, and are stored with different formats; as example IllustrisTNG stores a single snapshot in multiple \texttt{.hdf5} files, while in the \textsc{eagle} simulation \citep{eagle15} the snapshots are available for public download via an SQL web interface. Thus, the code requires these data products to be uniformed in a specific format, in order to be easily read and processed.   
      
   We release a beta version of the code, which can be read and improved by the scientific community with a request for access to its source through our website, at {www.astrodeep.eu/FORECAST}.
   
   
**How to use it?**
- Download the code github.com/FORECAST
- Download the required C/C++ libraries
- Download the hydrodynamical simulation input files
- Change paths and parameters in the input configuration file \*.ini (in each module)
- Compile Makefile
- run bash \*.sh file (one partition per time; refer to planes_list.txt file to choose the partition/snapshot-plane run) 
