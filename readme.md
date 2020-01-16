# C-MILOS


## Description 

Implementacion of MILOS in C 


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.


## Requeriments 

### Libraries

The following libraries and tools must be installed in your system: 

- [OpenMPI](https://www.open-mpi.org/) (Minor version tested 1.4-4)
- [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) (Minor version tested 3.3.4.0)
- [FFTW](http://www.fftw.org/)  (Minot version tested 3.3.3)
- [GSL](https://www.gnu.org/software/gsl/) (Minor version tested 1.13-3)
  
There are many differents ways to install them depending of OS what we are using. In our case we have been using Ubuntu 18.04 as OS, and these are the command to install each library, if you are in the same situation. For other OS, it's in your hands install the specific libraries.

OpenMPI: 

```
sudo apt-get update -y 
sudo apt-get install openmpi-bin
```

CFITSIO:

```
sudo apt-get update -y 
sudo apt-get install libcfitsio*
```

FFTW:

```
sudo apt-get update -y 
sudo apt-get install libfftw3-*
```

GSL:

```
sudo apt-get update -y 
sudo apt-get install libgsl*
```

### Files format

#### .per

#### .fits 

* Spectro 

The **fits** files used for pass to the program the spectro image must contain four dimensions: *number_of_wavelengths*X*number_stokes*X*number_rows*X*number_cols* . The order or these parameters may vary, but to identify each one the header of **fits** file must contain the type of each dimension with this correspondence:

  - Number of Wavelenghts: include CTYPE with the value **'WAVE-GRI'**
  - Number of Stokes: include CTYPE with the value **'STOKES  '**
  - Number of Rows: include CTYPE with the value **'HPLN-TAN'**
  - Number of Cols: include CTYPE with the value **'HPLT-TAN'**

An example can be this:

```
CTYPE1  = 'WAVE-GRI'
CTYPE2  = 'STOKES  '
CTYPE3  = 'HPLN-TAN' 
CTYPE4  = 'HPLT-TAN'  
```

* Wavelengths

If the observed spectra of all pixels use the same wavelength grid, the FITS file must contain a single, 2D array with dimension number of wavelength-points×2. The first column must contain the index with which the spectral line is identified according to the atomic parameter file.

#### .grid



#### .mod 



## Instalation


After install the pre-conditions libraries and download this repository you are ready to compile the code. This action will be do using the command make. there are the options that the command accepts:

* Compile and create executable **milos** 
```
make milos
```
* Compile and create executable **milosMPI**
```
make milosMPI
```
* Compile and create both: **milos** and **milosMPI**
```
make 
```
* Clean objects files and executable files. 
```
make clean
```

## Deployment

### milos

The sequential program must be controlled with a configuration file of type **.trol** . Inside the run folder, you can find an example of this type of file [cmilos.trol](run/cmilos.trol). We refer you to the pdf documentation to know in detail how each parameter works. 

The program must be executed by passing the configuration file as a parameter:

```
./milos run/cmilos.trol
```

### milosMPI

For the parallel program the execution is a bit different. In this case the file that we must pass as a parameter to the executable will be of type **.init** . As for the **.trol** file, we leave you an example inside the run folder [cmilos.init](run/cmilos.init)

The program must be executed using the command of **mpirun** or **mpiexec**. a simple use of parallelization on the local machine, would be to specify the number of processes with the *-np* option:

```
mpiexec -np 12 ./milosMPI run/cmilos.init
```
