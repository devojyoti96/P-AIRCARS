Introduction
============
P-AIRCARS (Polarimetry using Automated Imaging Routine for Compact Arrays of the Radio Sun)
*******************************************************************************************
..  image:: https://readthedocs.org/projects/p-aircars/badge/
    :alt: Documentation Status
    :scale: 100%
    :target: http://p-aircars.readthedocs.io

P-AIRCARS is an automated calibration and imaging routine for polarimetric calibration and imaging of solar observations done with Murchison WIdefield Array (MWA) https://www.mwatelescope.org. P-AIRCARS has been developed and maintained by solar physics group at National Centre for Radio Astrophysics, Tata Institute of Fundamental Research (NCRA-TIFR), Pune, India https://www.ncra.tifr.res.in.

P-AIRCARS uses CASA (Common Astronomy Software Application) https://casa.nrao.edu for intensity and bandpass self-calibration and polarisation calibration is performed by our own implementation of the full Jones calibration algorithm described by Mitchell et al. 2008 https://doi.org/10.1109/JSTSP.2008.2005327.

Imaging at P-AIRCARS is performed by WSClean https://wsclean.readthedocs.io. If WSClean is not installed, P-AIRCARS performs imaging by CASA. When CASA is used for imaging the computation speed is slow.

Basic philosophy of P-AIRCARS is the self-calibration and use the instrumental model to perform the precise calibration. Details of the algorithm and implementation can be found in the follwoing papers

	1.AIRCARS (Mondal et al. 2020) https://doi.org/10.3847/1538-4357/ab0a01

	2.Kansabanik et al. 2021a, in preparation 

	3.Kansabanik et al. 2021b, in preparation

P-AIRCARS has several modules and scripts. Modules can be useful for any astronomical self-calibration. Users want to use the modules please the the documentation of the Module details.
Instructions for local machines, laptops and work stations can be found in the installation section. Implementation for high performance computing (HPC) environment will be available shortly. A basic tutorial is also included.

For any queires and issues reach us at:

dkansabanik@ncra.tifr.res.in, paircarsnotification@gmail.com 

Software requirements
*********************

P-AIRACRS has been tested in the follwoing linux environments

1. CentOS 7

2. Ubuntu 20.04

Installation
============
P-AIRCARS is a self-contained python package. It is build on python3. Python version greater than 3.6 is required for P-AIRACRS. It is recommended to use python virtual environment to install P-AIRCARS, but it is not a necessary requirement. Installation steps are as follows.

Creating virtual environment
****************************
Install python virtual environment:

>> python3 -m pip install --upgrade pip
 
>> python3 -m pip install virtualenv

Activate virtual environment:

>> python3 -m venv /path/to/virtualenv/paircars

>> source /path/to/virtualenv/paircars

Check the virtual environment path:

>> which python3

It should show "/path/to/virtualenv/paircars/bin/python3"

Obtaining P-AIRCARS
*******************
P-AIRCARS source code can be downloaded from https://github.com/devojyoti96/P-AIRCARS.git. 

P-AIRCARS is not public now. It will be public very soon. If you want to use P-AIRCARS before it becomes public, please reach us at dkansabanik@ncra.tifr.res.in, paircarsnotification@gmail.com

Installing P-AIRCARS
********************
Move to the P-AIRCARS directory.

Run *setup.py*

>> python3 setup.py install

That's all. It will automatically install all required packages and libraries. This installation process is in very early stage. If you find any issues during installation, please reach us at dkansabanik@ncra.tifr.res.in, paircarsnotification@gmail.com. 

.. note::
   If you installed P-AIRCARS in a virtual environment, please make sure you are inside the virtual environment. P-AIRCARS can not be accessed outside the virtual environment in this case.

   >> which python3

   It should show "/path/to/virtualenv/paircars/bin/python3". If you are not in virtual environment, please activate it. 

Documentation
=============
Details documentation can be found add at https://p-aircars.readthedocs.io/
