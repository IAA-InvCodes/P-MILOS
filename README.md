# P-MILOS

Authors: Manuel Cabrera, Juan P. Cobos, Luis Bellot Rubio (IAA-CSIC).

For questions, please contact Luis Bellot (lbellot@iaa.es)

This development has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 824135 (SOLARNET).

## Description 


This repository contains P-MILOS, a parallel implementacion of the MILOS inversion code using OpenMPI.  P-MILOS is capable of inverting full Stokes spectropolarimetric data cubes in real time, providing the one-component Milne-Eddington model atmosphere that best fit the observations. In this page we provide a quick overview about how to install the required libraries, how to compile the code, the input/output files, and how to execute the code sequentially and in parallel. An extended user manual can be found [here](p-milos_manual.pdf). 


## Requeriments 

### Libraries

The following programs and libraries must be installed in your system to run P-MILOS: 

- [OpenMPI](https://www.open-mpi.org/) (Oldest version tested 1.4-4)
- [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) (Oldest version tested 3.3.4.0)
- [FFTW](http://www.fftw.org/)  (Oldest version tested 3.3.3)
- [GSL](https://www.gnu.org/software/gsl/) (Oldest version tested 1.13-3)
  
There are many different ways to install the libraries depending on your operating system.  On Ubuntu 18.04, these are the commands we used. 

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

### Input/output files

#### Profile files (.per)

ASCII files with extension .per contain a set of Stokes profiles. They have the same format as SIR profile files. They are used as input when inverting one pixel and as output when synthesizing the profiles from a given model atmosphere. 

Profile files have one row per wavelength position and 6 columns containing:

* The index of the spectral line in the atomic parameter file
* The wavelength offset with respect to the central wavelength (in mA) 
* The value of Stokes I
* The value of Stokes Q
* The value of Stokes U
* The value of Stokes V

This is an example of a profile file containing the Stokes parameters of spectral line number 1 in 30 wavelength positions, from -350 to + 665 mA:

```
1	-350.000000	9.836711e-01	6.600326e-04	4.649822e-04	-3.694108e-03
1	-315.000000	9.762496e-01	1.186279e-03	8.329745e-04	-6.497371e-03
1	-280.000000	9.651449e-01	2.305113e-03	1.581940e-03	-1.135160e-02
1	-245.000000	9.443904e-01	5.032997e-03	3.333831e-03	-2.191048e-02
1	-210.000000	9.018359e-01	1.146227e-02	7.265856e-03	-4.544966e-02
1	-175.000000	8.222064e-01	2.265146e-02	1.368713e-02	-8.623441e-02
1	-140.000000	7.066048e-01	3.263217e-02	1.847884e-02	-1.242511e-01
1	-105.000000	5.799722e-01	3.157282e-02	1.560218e-02	-1.238010e-01
1	-70.000000	4.711627e-01	2.068015e-02	6.887295e-03	-8.728459e-02
1	-35.000000	4.014441e-01	9.837587e-03	-1.054865e-03	-4.476189e-02
1	-0.000000	3.727264e-01	4.631597e-03	-4.830483e-03	-9.482273e-03
1	35.000000	3.799767e-01	5.985593e-03	-3.846331e-03	2.321622e-02
1	70.000000	4.249082e-01	1.378049e-02	1.806246e-03	6.157872e-02
1	105.000000	5.119950e-01	2.571598e-02	1.073034e-02	1.046772e-01
1	140.000000	6.316050e-01	3.361904e-02	1.782490e-02	1.297000e-01
1	175.000000	7.571660e-01	2.941514e-02	1.718121e-02	1.115998e-01
1	210.000000	8.600360e-01	1.762856e-02	1.087526e-02	6.786425e-02
1	245.000000	9.230015e-01	8.213106e-03	5.305210e-03	3.366419e-02
1	280.000000	9.545796e-01	3.605521e-03	2.427676e-03	1.654710e-02
1	315.000000	9.701367e-01	1.734934e-03	1.205792e-03	9.068003e-03
1	350.000000	9.786569e-01	9.418463e-04	6.697733e-04	5.552034e-03
1	385.000000	9.838719e-01	5.600237e-04	4.053978e-04	3.671347e-03
1	420.000000	9.873039e-01	3.608273e-04	2.646724e-04	2.540035e-03
1	455.000000	9.896545e-01	2.543154e-04	1.872806e-04	1.792397e-03
1	490.000000	9.912774e-01	1.968866e-04	1.437179e-04	1.263287e-03
1	525.000000	9.923597e-01	1.690016e-04	1.205914e-04	8.552026e-04
1	560.000000	9.929766e-01	1.638180e-04	1.128983e-04	4.949613e-04
1	595.000000	9.930763e-01	1.826333e-04	1.215294e-04	1.047099e-04
1	630.000000	9.923094e-01	2.385952e-04	1.569032e-04	-4.666936e-04
1	665.000000	9.895550e-01	3.734447e-04	2.531845e-04	-1.597078e-03
```


#### Profile files (.fits) 


P-MILOS can be fed with data cubes containing the Stokes profiles observed over the whole field of view. This is the usual way of inverting data from narrow-band filter imagers such as CRISP. 

The data cubes must be written in FITS format, with one cube containing one spectral scan. They must have four dimensions corresponding to the number of rows, the number of columns, the number of observed wavelengths, and the number of Stokes parameters. These dimensions can appear in any order. The exact order is specified in the FITS header by means of the keywords CTYPE1, CTYPE2, CTYPE3 and CTYPE4. P-MILOS follows the SOLARNET standard:

  -  **HPLN-TAN** indicates a spatial coordinate dimension
  - **WAVE-GRI** indicates a wavelength dimension
  - **STOKES  '** indicates the Stokes parameter dimension

The following example corresponds to a data cube with the x spatial coordinate in the first dimension, the y coordinate in the second dimension, wavelength in the third dimension and the polarization in the fourth dimension. 

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
1                       :        -350,            35,           665
```
In the file [malla.grid](run/malla.grid) you can find an extended example. 


#### .mod 

These files will be used for three purposes:

  1. For specify the initial model of a synthesis.  
  2. For specify the initial model of a inversion. 
  3. For save output model when we are doing the inversion of a profile stored in a .per file. 

The order of parameters in the file must be always the same. This is an example: 

```
eta_0:          14
magnetic field: 1200
LOS velocity:   0
Doppler width:  0.07
damping:        0.05
gamma:          130
phi:            25
S_0:            0.25
S_1:            0.75
v_mac:          1
filling factor: 1
```


## Instalation

In order to deploy the application, it must first be compiled on the target machine. To do this, you must use the command line option 'make' from same directory where the source code is located. So, the first thing is to position ourselves in the P-MILOS. The code is compiled by default using single precision, but if you want you can compile the code using double precision by specifying the following variable 'use_double=yes' when you run make.

The other options are to choose which version of the code to compile (sequential or parallel), or whether to clean up the generated object code and executables. 

* Compile and create executable **milos** 
```
make milos
```
* Compile and create executable **pmilos**
```
make pmilos
```
* Compile and create both: **milos** and **pmilos**
```
make 
```
* Clean objects files and executable files. 
```
make clean
```

## Deployment


### milos

The sequential program must be controlled with a configuration file of type **.mtrol** . Inside the run folder, you can find an example of this type of file [pmilos.mtrol](run/pmilos.mtrol). We refer you to the pdf documentation to know in detail how each parameter works. 

The program must be executed by passing the configuration file as a parameter:

```
./milos run/pmilos.mtrol
```

### pmilos

For the parallel program the execution is a bit different. In this case the file that we must pass as a parameter to the executable will be of type **.init** . As for the **.trol** file, we leave you an example inside the run folder [pmilos.minit](run/pmilos.minit)

The program must be executed using the command of **mpirun** or **mpiexec**. a simple use of parallelization on the local machine, would be to specify the number of processes with the *-np* option:

```
mpiexec -np 12 ./pmilos run/pmilos.minit
```
