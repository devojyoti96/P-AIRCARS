Create and Activate Conda Environment
=====================================
This guideline provide how to activate **conda** environment using **Mamba** before installing and using P-AIRCARS. 

Create conda environment in default conda directory
----------------------------------------------------
This will create the conda environment in the default conda directory where conda is installed. If a custom directory is chosen during installation, environment will be created in that directory.

Create conda environment
~~~~~~~~~~~~~~~~~~~~~~~~

1. **To create a Python 3.10 environment with compaitable C/C++ compilers:**

.. code-block:: bash
   
   mamba create -n paircars_env --override-channels -c conda-forge python=3.10 gcc_linux-64=14 gxx_linux-64=14 gfortran_linux-64=14 cmake pkg-config pip
   
   
Activate and deactivate conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **To activate the environment:**

.. code-block:: bash

   conda activate paircars_env
   
.. note ::

   For the first time, after activating conda environment, type ``python -c "import sys; print('\n'.join(sys.path))"``. This should not show any local paths. If it shows local paths, conda environment is leaking into local python environements, which may cause version conflicts.  

2. **To deactivate:**

.. code-block:: bash

   conda deactivate


