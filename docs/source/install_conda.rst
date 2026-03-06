Install Conda: Miniconda
========================

This guide provides instructions for installing **Miniconda** into a **custom directory**, without conflicting with local site packages, creating a Python environment, and activating it.

Overview
--------

- **Miniconda**: A minimal Conda installer (recommended for lightweight setups).

.. note::

   Miniconda is recommended for users who want to install only the packages they need.

Install Miniconda in a Custom Directory
---------------------------------------

1. **Download the Installer**

   Choose your OS and download from the official Miniconda page:

   https://docs.conda.io/en/latest/miniconda.html

   For Linux:

   .. code-block:: bash

      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      
2. **Setup environment variables**
    
   Setup these environment variables to avoid conflict with local site packages.
   
   .. code-block:: bash
   
      export PYTHONNOUSERSITE=1
      unset PYTHONPATH 
    

3. **Install into a Custom Directory**

   Replace `/path/to/miniconda3_custom` with your desired location:

   .. code-block:: bash

      bash Miniconda3-latest-Linux-x86_64.sh -b -p /path/to/miniconda3_custom
      
   .. note::
      
      In HPC cluster, it is recommended to set **/path/to/miniconda3_custom** to a location which is accessible by all nodes and with fast disk speed. Read you HPC documentation carefully and check whether **conda** is already installed and available as **module** or not. 
      
4. **Enable the 'conda' Command**
   
   .. important::
      To avoid using the full path (`/path/to/miniconda3_custom/bin/conda`) every time, do the follwing:
   
   .. code-block:: bash

      /path/to/miniconda3_custom/bin/conda init
      source ~/.bashrc    # or ~/.zshrc, depending on your shell

   After this, `conda` will be available as a global command.
   
 
