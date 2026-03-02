Run P-AIRCARS Pipeline
=======================

Basic run
---------
To run P-AIRCARS pipeline, with default settings for full analysis, run the following command from terminal. Work directory needs not to be created before hand, but the path where it will be created should exist.

.. code-block :: bash

    run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> 
    
.. note ::

   User needs to provide full path. Providing short path or only the directory name will cause error.    
        
Advanced run
------------
For advanced run, user is requested to first check the parameters of **run-mwa-paircars**.

.. code-block :: bash
 
    run-mwa-paircars -h
    
.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-paircars -h
   
Multiple options demonstrated below can be combined to have all of them together.
   
1. To view details of measurement set:

.. code-block :: bash

    show-paircars-ms </full/path/to/measurement_set>
    
Runs with advanced calibration paramaters 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. If one do not have calibrator observations:

.. code-block :: bash

   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> 
    
.. note ::

   P-AIRCARS will self-calibrate and perform flux-density and polarisation calibration based on some assumptions. In these cases, absolute values of the quantities should be carefully considered.
   

2. Do calibration with custom calibration parameters. There are two parameters: **cal_uvrange** and **solint** which can be changed. Example, run the following command to perform gain solutions at 10second interval and >200lambda data:

.. code-block :: bash
    
   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --solint "10s" --uvrange ">200lambda" 
    
3. By default for full-polar data, polarization calibration will be performed. To disable it:

.. code-block :: bash 

   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --no_polcal 
    
    
Runs with advanced imaging paramaters 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
1. Run pipeline to image specific time and frequency range. Default is to use entire time and frequency range. Example for imaging two time ranges given in UTC and frequency ranges given in MHz: 

.. code-block :: bash

   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --timerange "2024/06/10/09:00:00~2024/06/10/09:30:00,2024/06/10/10:15:00~2024/06/10/10:45:00" --freqrange 100~150,200~230 
    
2. Run imaging with custom time and frequency resolution. Default is to use each 1.28 MHz coarse bands and 10s of integration. Example run for imaging at 0.5 second time resolution and 160 kHz (0.16 MHz) frequency resolution:

.. code-block :: bash 
    
   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --image_timeres 0.5 --image_freqres 0.16 
    
2. Default is to make only Stokes I images if `do_polcal=False` and Stokes IQUV, if `do_polcal=True`. To run only Stokes I imaging, even if `do_polcal=True`, run:

.. code-block :: bash
    
   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --pol I 
    
Similarly, all other advanced imaging parameters can be used.

Switching off particular pipeline step(s)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
By default, all steps will be performed by pipeline. Even pipeline was run upto certain stages and then stopped, all steps from beginning will be performed to avoid any potential issue in failure in previous runs. If user is certain that previous run was successful upto certain stages, those stages can be switched.

.. caution :: 
    
    User should not modify any file and directory structure in the work directory. Switching off certain parameters will only allow to run the pipeline forward, if the expected output products from those steps are present with appropriate name in appropriate directory. Otherwise, it will fail.

Take a look at the **Advanced pipeline parameters** in the help page of **run-mwa-paircars**. Each parameters are self explanatory. Some examples are given below:

1. To switch off self-calibration:

.. code-block :: bash
    
   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --no_selfcal 
    
2. To stop final imaging:

.. code-block :: bash
    
   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --no_imaging 
   
3. To switch off self-calibration and final imaging

.. code-block :: bash

   run-mwa-paircars </full/path/to/data_directory> </full/path/to/data_metafits> --workdir </full/path/to/work_directory> --outdir </full/path/to/output_product_directory> --cal_datadir </full/path/to/calibrator_data_directory> --cal_metafits </full/path/to/calibrator_metafits> --no_selfcal --no_imaging 


