# C-MILOS


## Description 

This repository contains an implementacion of MILOS in C and will get you a copy of the project up and running on your local machine for development and testing purposes. An extended user manual can be found [here](c-milos_manual.pdf). But in this page you can find a quick overview about how to install the necessary libraries, the types of files used and how to use the programs (sequential and parallel). 


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

It's used to specify one profile. We will use it and input for inversion one pixel and as output for synthesis of one model.
It contains 6 columns:

* The first is the index of the spectral line used in the spectral lines file.
* The second is the offset of wavelenght respect the central wavelenght. 
* Value of I
* Value of Q
* Value of U
* Value of V

This is an example of one line: 

```
1	-0.350000	8.480462e-01	2.081567e-05	-3.810591e-05	-2.589682e-04
```


#### .fits 

* Spectro 

The **fits** files used for pass to the program the spectro image must contain four dimensions: *number_rows*X*number_cols*number_of_wavelengths*X*number_stokes*X* . The order or these parameters cannot change and for identify each one the header of **fits** file must contain the type of each dimension with this correspondence:

  - Number of Rows: include CTYPE with the value **'HPLN-TAN'**
  - Number of Cols: include CTYPE with the value **'HPLT-TAN'**
  - Number of Wavelenghts: include CTYPE with the value **'WAVE-GRI'**
  - Number of Stokes: include CTYPE with the value **'STOKES  '**

An example can be this:

```
CTYPE1  = 'HPLN-TAN' 
CTYPE2  = 'HPLT-TAN' 
CTYPE3  = 'WAVE-GRI'
CTYPE4  = 'STOKES  ' 
```

* Wavelengths

If the observed spectra of all pixels use the same wavelength grid, the FITS file must contain a single, 2D array with dimension number of wavelength-points×2. The first column must contain the index with which the spectral line is identified according to the atomic parameter file.

* Output Models 

For save the output models of invert one image, the program use FITS. The data is saved in FLOAT precision and the dimensiones of image will be: numberOfRows X numberOfCols X 13. The number 13 comes from the eleven parameters of the model, the number of interations used by the algorithm to found the solution in that pixel and the value of Chisqr calculated for the result model of that pixel respect the input profile. Therefor, the order of the third dimension of the file will be: 

  1. eta0 = line-to-continuum absorption coefficient ratio         
  2. B = magnetic field strength       [Gauss]
  3. vlos = line-of-sight velocity     [km/s]         
  4. dopp = Doppler width              [Angstroms]
  5. aa = damping parameter
  6. gm = magnetic field inclination   [deg]
  7. az = magnetic field azimuth       [deg]
  8. S0 = source function constant
  9. S1 = source function gradient
  10. mac = macroturbulent velocity     [km/s]
  11. alpha = filling factor of the magnetic component [0->1]
  12. Number of iterations needed. 
  13. Value of Chisqr. 

#### .grid

This is the file where you can specify the number of line from your file with spectral lines to use and the range of wavelenghts to use. This range will be specify with an initial lambda, a step between each wavelenght and the final lambda of the range. 

Look this example:

```
Line indices            :   Initial lambda     Step     Final lambda
(in this order)                    (mA)          (mA)         (mA) 
-----------------------------------------------------------------------
1                       :        -350            35           665
```
In the file [malla.grid](run/malla.grid) you can find an extended example. 


#### .mod 

These files will be used for three purposes:

  1. For specify the initial model of a synthesis.  
  2. For specify the initial model of a inversion. 
  3. For save output model when we are doing the inversion of a profile stored in a .per file. 

The order of parameters in the file must be always the same. This is an example: 

```
INITIAL_MODEL_ETHA0         :20
INITIAL_MODEL_B             :1100
INITIAL_MODEL_VLOS          :0.2
INITIAL_MODEL_LAMBDADOPP    :0.03
INITIAL_MODEL_AA            :0.05
INITIAL_MODEL_GM            :120
INITIAL_MODEL_AZI           :150
INITIAL_MODEL_S0            :0.35
INITIAL_MODEL_S1            :0.5
INITIAL_MODEL_MAC           :0
INITIAL_MODEL_ALFA          :1
```


## Instalation

In order to deploy the application, it must first be compiled on the target machine. To do this, you must use the command line option 'make' from same directory where the source code is located. So, the first thing is to position ourselves in the C-MILOS. The code is compiled by default using single precision, but if you want you can compile the code using double precision by specifying the following variable 'use_double=yes' when you run make.

The other options are to choose which version of the code to compile (sequential or parallel), or whether to clean up the generated object code and executables. 

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
