Install Conda: Mamba
====================
This guide provides instructions for installing **Mamba** into a **custom directory**, without conflicting with local site packages, creating a Python environment, and activating it.

Overview
--------
- **Mamba** is fast envrionment solver. We recommend using this.

  .. note::

     User can use standard **anaconda** or **miniconda** as well. 

Install Mamba in a Custom Directory
---------------------------------------
1. **Download the Installer**

   .. code-block:: bash

      wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh    

2. **Install into a Custom Directory**

   Replace `/path/to/conda3_custom` with your desired location:

   .. code-block:: bash

      bash Miniforge3-Linux-x86_64.sh -b -p /path/to/conda3_custom
      
   .. note::
      
      In HPC cluster, it is recommended to set **/path/to/conda3_custom** to a location which is accessible by all nodes and with fast disk speed. Read you HPC documentation carefully and check whether **conda** is already installed and available as **module** or not. 
      
4. **Enable the 'conda' Command**
   
   .. important::
      To avoid using the full path (`/path/to/conda3_custom/bin/conda`) every time, do the follwing:
   
   .. code-block:: bash

      /path/to/conda3_custom/bin/conda init
      source ~/.bashrc    # or ~/.zshrc, depending on your shell

   After this, `conda` will be available as a global command.
   
 
