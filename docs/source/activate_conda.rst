Create and Activate Conda Environment
=====================================
This guideline provide how to activate **conda** environment before installing and using P-AIRCARS. 

Create conda environment in default conda directory
----------------------------------------------------
This will create the conda environment in the default conda directory where conda is installed. If a custom directory is chosen during installation, environment will be created in that directory.

Create conda environment
~~~~~~~~~~~~~~~~~~~~~~~~

1. **To create a minimal Python 3.10 environment:**

.. code-block:: bash

   conda create -n paircars_env python=3.10
   
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

Create conda environment in custom conda directory
---------------------------------------------------
To create in a custom directory ``</path/to/env>``, follow the steps below.

.. tip ::

    It is recommended to install in custom path in HPC architechture or your default conda path has limited disk space or does not have global access.

Create conda environment
~~~~~~~~~~~~~~~~~~~~~~~~

1. **To create a minimal Python 3.10 environment:**

.. code-block:: bash

   conda create -p </path/to/env>/paircars_env python=3.10
   
Activate and deactivate conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **To activate the environment:**

.. code-block:: bash

   conda activate </path/to/env>/paircars_env

2. **To deactivate:**

.. code-block:: bash

   conda deactivate

