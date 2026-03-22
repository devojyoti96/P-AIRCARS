Initial Setup
=============
After installation of **P-AIRCARS**, before running the pipeline, some initial setup is needed. These include downloading some required metadata for the pipeline. **P-AIRCARS** also provides multiple ways to monitor its progress -- local GUI logger (for local machine only), over email, web-based remote logger, and prefect dashboard.  

Download P-AIRCARS metadata
------------------------------
1. To download and save the required P-AIRCARS metadata in appropriate directory, run from command line:

   .. code-block :: bash
    
      init-paircars-setup --init
    
   .. admonition:: Click here to see parameters
      :class: dropdown

      .. program-output:: init-paircars-setup -h
    
2. If data files are present, but needs to updated, run:

   .. code-block :: bash

       init-paircars-setup --init --update
   
.. important::

   ``prefect`` server will be automatically setup during this process. In local environment, if there is any issue in starting ``prefect`` server, ``prefect`` will automatically fall back to its ephemeral mode. For cluster environment, P-AIRCARS can not be run without ``prefect`` server. By default, ``prefect`` server uses the **port 4260** and its associated PostgreSQL database used **port 5260**. If these port are pre-occupied, P-AIRCARS will try to kill those port to avoid port overload. If not successful, the closest free port will be used. The instructions to access ``prefect`` dashboard will be displayed in the terminal with the steps to access it. If user do not want to close the port, use ``--no_kill_port`` parameter during initiation. 
   

Custom P-AIRCARS metadata directory
-----------------------------------
By default, data directory will be at "~/.paircarspipe/paircarspipe_data". It requires 20 GB space. Sometimes home directory may not have sufficient space. In that case, one can setup data directory in a custom location as follows:

.. code-block :: bash

    init-paircars-setup --init --datadir </full/path/to/custom/datadir>
    
.. note::
   
   This is not data directory for observations. We strongly suggest not to change anything in this directory after P-AIRCARS initiation. 
    
Setup e-mail ids
----------------
To receive remote logger Job ID and password, and e-mail notifications for pipeline progress, use can setup their e-mail id(s) in P-AIRCARS. 

.. code-block :: bash

    init-paircars-setup --init --emails <youremail1@email1.id1>,<youremail2@email2.id2> 
    
If you setup a remote logger as described below, you will receive a unique Job ID and password (user provided or auto-generated six-character) to access logs of a particular pipeline run from the remote logger. Without this password, one can not access logs of that particular pipeline run. This added security as well as privacy when multiple user uses the same remote logger link, for example, an institute based remote logger link.   
    
Setup remote logger link
-------------------------
If remote logger is intended to be used, setup the remote link in P-AIRCARS metadata. By default, P-AIRCARS will set an auto-generated 6-character password for you.

.. code-block :: bash
    
    init-paircars-setup --init --remotelink https://<remote-logger-name>.onrender.com
 
Setup custom remote password
----------------------------
If user wants to set a custom password for them

.. code-block :: bash
    
    init-paircars-setup --init --remotelink https://<remote-logger-name>.onrender.com --remote_password <your-custom-password>

Update remote logger link and/or e-mail ids
-------------------------------------------
If user wants to update the already provided remote logger link or e-mail id(s), simply run the above commands with new values. P-AIRCARS will automatically update the database with these new values.


Setup remote logger
-------------------
Before using remote logger this, create your own remote logger on free-tier cloud platform, https://render.com. One can use, same **remotelink** in multiple machines and users. However, free-tier link has some limitations on bandwidth and concurrency. If you want to use **remotelink** for your institution, we suggest to purchase suitable paid version or setup seperate **remotelink** for different users.

.. toctree::
   :maxdepth: 2
   
   create_remotelogger












    
