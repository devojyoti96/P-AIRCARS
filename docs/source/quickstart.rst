Quickstart
==========
P-AIRCARS is distributed on
`PyPI <https://pypi.org/project/paircars/>`__. To use it:

1. Create conda environment with python 3.10 with compaitable C/C++ libraries

   .. code-block:: bash

      conda create -n paircars_env --override-channels -c conda-forge python=3.10 gcc_linux-64=14 gxx_linux-64=14 gfortran_linux-64=14 cmake pkg-config pip
      
      conda activate paircars_env
      
  .. note::
  
     We suggest using **Mamba** for fast conda installtion and environment creation.

2. Install P-AIRCARS in conda environment

   .. code-block:: bash

      pip install paircars

3. Initiate necessary post-installation setup for metadata 

   .. code-block:: bash

      init-paircars-setup --init 
      
   .. note::
   
      By default, the necessary data will be saved in home directory and requires about 20 GB of disk space. We suggest using any other location with larger disk space and specify that by ``--datadir </full/path/to/paircars_datadir>`` in the above command.
      

4. Before running the pipeline, setup your data as following:
    
    * Create a <target_datadir> and put all coarse channel measurement sets of solar scan of a single observation ID (OBSID) inside it.
    
    * Create a <cal_datadir> and put all coarse channel measurement sets for calibrator observation of a single OBSID inside it.
    
5. Run P-AIRCARS pipeline

   .. code-block:: bash

      run-mwa-paircars <full path of target measurement set directory> <full path of target metafits file> --cal_datadir <full path of calibrator measurement set directory> --cal_metafits <full path of calibrator metafits> --workdir <full path of work directory> --outdir <full path of output products directory>

.. note ::

   Always provide the entire direcotry path. Short path or only directory name may cause errors. Keep target measurement sets for a single OBSID and calibrator measurement sets for a single OBSID must be kept in seperate directories. If calibrator is not present, do not provide these information.

That’s all. You started P-AIRCARS pipeline for analysing your MWA solar observation 🎉. Read the ``Directory Structure and Data Products`` section to understand how to find final images.

6. To see all running P-AIRCARS jobs

   .. code-block :: bash
        
      show-paircars-status --show
      
   
7. If P-AIRCARS is running in local machine, view local log of any job using the <jobid>:

   .. code-block :: bash
    
      run-mwa-mwalogger --jobid <jobid>
      
.. note::

   If you are running P-AIRCARS is cluster environment, first checkout **HPC Settings** in the document for viewing P-AIRCARS log remotely using prefect dashboard.
      
      
8. Output products will be saved in : ``<path of output products directory>``.

