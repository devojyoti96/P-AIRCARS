Introduction
============
P-AIRCARS (Polarimetry using Automated Imaging Routine for Compact Arrays of the Radio Sun)
*******************************************************************************************

P-AIRCARS is an automated calibration and imaging routine for polarimetric calibration and imaging of solar observations done with Murchison Widefield Array (MWA) https://www.mwatelescope.org. P-AIRCARS has been developed and maintained by solar physics group at National Centre for Radio Astrophysics, Tata Institute of Fundamental Research (NCRA-TIFR), Pune, India https://www.ncra.tifr.res.in.

P-AIRCARS uses CASA (Common Astronomy Software Application) https://casa.nrao.edu for intensity and bandpass self-calibration and polarisation calibration is performed by our own implementation of the full Jones calibration algorithm described by Mitchell et al. 2008 https://doi.org/10.1109/JSTSP.2008.2005327.

Imaging at P-AIRCARS is performed by WSClean https://wsclean.readthedocs.io. If WSClean is not installed P-AIRCARS performs imaging by CASA. When CASA is used for imaging the computation speed is slow.

Basic philosophy of P-AIRCARS is self-calibration and use the instrumental model to perform the precise solar calibration. Details of the algorithm and implementation can be found in the follwoing papers

	1.AIRCARS (Mondal et al. 2019) https://doi.org/10.3847/1538-4357/ab0a01

	2.Kansabanik et al. 2021a, in preparation 

	3.Kansabanik et al. 2021b, in preparation

P-AIRCARS has several modules and scripts. Modules can be useful for any astronomical self-calibration. Users want to use the modules please the the documentation of the Module details.
Instructions for local machines, laptops and work stations can be found in the installation section. Implementation for high performance computing (HPC) environment will be available shortly. A basic tutorial is also included.

For any queires and issues reach us at:

dkansabanik@ncra.tifr.res.in, paircarsnotification@gmail.com 

System requirements
*******************
++++++++++++++++++++++++++++++
Minimum hardware requirements
++++++++++++++++++++++++++++++
CPUs: 4 cores 

RAM : 4 GB 

Installation space : 700 MB

Disk space : 6 times of the data volume

+++++++++++++++++++++++++
Workstation configuration
+++++++++++++++++++++++++
P-AIRCARS has been tested on two types of workstations with following configurations.

1. Intel Xeon Gold 6138 2 GHz 40 CPUs - 256 GB RAM

   4 x 10TB (RAID-0) 40TB , 2 sets of 3 x 10 TB (RAID-5) , 2 x 20 TB 

2. Boston Intel E5-2699 v3 2.30GHz, 72 CPUs

	256 GB (8 x 32) RAM - 2 x 300 GB SSD System

	Data 41 TB (12 x 4 TB SATA with RAID-5, 256KB stripaize)
	
	Data1 81 TB (12 x 8 TB SATA with RAID-5, 512KB stripaize)
	
	960 GFLOPS

++++++++++++++++++++
Software environment
++++++++++++++++++++
P-AIRACRS has been tested in the follwoing linux environments

1. CentOS 7

2. Ubuntu 20.0.4

Details of other software requirements is available in installation section.
