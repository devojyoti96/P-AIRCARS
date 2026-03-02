P-AIRCARS Documentation
=======================

**P-AIRCARS** is an automated calibration and imaging pipeline designed for solar radio observations using **MWA** radio telescope. It performs end-to-end calibration, flagging, and imaging with a focus on dynamic solar data, supporting both spectral and temporal flexibility in imaging products.

Introduction
------------

Solar radio data presents unique challenges due to the high variability and brightness of the Sun, as well as the need for high time-frequency resolution. The P-AIRCARS pipeline addresses these challenges by:

- Automating the calibration of interferometric data, including flux, phase, and polarization calibrations
- Supporting time-sliced and frequency-sliced imaging workflows
- Pipeline workflow is managed by ``prefect`` leveraging Dask for scalable parallel processing
- Providing hooks for integration with contextual data from other wavelegths for enhanced solar analysis

Features of P-AIRCARS
----------------------

P-AIRCARS serves as a reference pipeline for science-ready processing of **P-AIRCARS solar observations**. It is designed to:

- Calibrate and image non-trivial solar observations
- Enable reproducible science through consistent and documented data reduction steps
- Designed to use on personal computer, single-node workstations, as well as in high-performance computing multi-node cluster
- A free-tier cloud-based remote logger to monitor pipeline over the internet

Software environment
--------------------
P-AIRCARS is tested on operating systems, Ubunut 22 and Ubunut 24, and under Python 3.10 environment. Using P-AIRCARS in other operating systems and python version does not gurantee successful run and limited support for debugging is available. User may look at **Containersed Use** section in those scenarios. 

Sample dataset
---------------

User can download and test entire P-AIRCARS pipeline using the sample dataset available in Zenodo: https://doi.org/10.5281/zenodo.18641232. Do not use this sample dataset for any publication without permission from the developer.

Contents
---------
.. toctree::
   :maxdepth: 2
   
   quickstart
   install
   initial_setup
   tutorial
   output
   kill
   hpc
   cli
   paircars
   dev
   ack
   lic
   


