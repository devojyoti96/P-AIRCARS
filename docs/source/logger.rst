Monitoring P-AIRCARS Pipeline Jobs
==================================
P-AIRCARS pipeline jobs are run in background and in parallel as ``prefect`` tasks managed by the ``prefect`` workflow. P-AIRCARS provides multiple ways to monitor its progress:

1. **Over email**: If e-mail ids are setup, pefect will provide progress report of the P-AIRCARS job over email, which are sufficient to basic monitoring of the pipeline. 

2. **Local GUI logger**: If P-AIRCARS is running in local environment, a local GUI logger is available to monitor logs.

.. note ::
 
   This GUI logger is not available in HPC cluster environment.

3. **Remote logger**: If remote logger is setup, P-AIRCARS logs can be monitored using the web-based remote logger.

4. **Prefect dashboard**: P-AIRCARS uses ``prefect`` server, and all P-AIRCARS jobs on a single machine under a single user can be monitored from a single place using ``prefect`` dashboard.

.. admonition :: Recommendation

   We recommend using ``prefect`` dashboard for monitoring, it ``prefect`` is setup properly.


1. If ``prefect`` server is running, progress of the pipeline can be monitored from ``prefect`` dashboard. To view how to connect to ``prefect`` dashboard:

   .. code-block :: bash
        
      init-paircars-prefect status
    
2. If P-AIRCARS is not running with ``prefect`` server mode in local environment, use local logger GUI.

.. toctree::
   :maxdepth: 2
   
   locallogger

.. warning ::  
  
   If you are working on a headless machine, such as work station, remote server or high-performace computing node, we recommend not to use local GUI logger, as it may encounter issue in deploying the PyQT based GUI window.
   
3. User may use remote logger:
   
.. toctree::
   :maxdepth: 2
   
   remotelogger

