Slurm environment
------------------

Log in into cluster
~~~~~~~~~~~~~~~~~~~
First login into login node of the cluster. Install P-AIRCARS and initiate, such that ``prefect`` server starts in the login node.

.. code-block:: bash

   init-paircars-setup --init --remotelink <remotelink> --emails <emails> --datadir </full/path/to/fast/datadir> 

.. note::

   ``prefect`` server consumes very small amount of CPU and memory. Hence, it is okay to setup it in login node, which may not be calculated for resource utilisation. Otherwise, if it started in compute nodes, it will be counted in resource utilisation and also closed after the maximum timelimit of the compute node. 
   
    
Gather cluster specifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First step is to gather slurm cluster informations. Crucial informations are, partition name, account name (if required by the cluster), node name list, and time limit. Use the following command:

.. code-block:: bash

   sinfo
   
Run P-AIRCARS from login node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
After the setup is done, run basic P-AIRCARS command as follows. Do not forget to include ``--cluster`` and ``--partition`` keywords, and ``--account`` keyword if required by the cluster. 

.. code-block:: bash

   run-mwa-paircars <full path of target measurement set directory> <full path of target metafits file> --cal_datadir <full path of calibrator measurement set directory> --cal_metafits <full path of calibrator metafits> --workdir <full path of work directory> --outdir <full path of output products directory> --cluster --partition <partition_name>
   
   
.. note:: 
   
   Sometimes loading different modules and submodules takes time in cluster environment. Do not loos patience, until P-AIRCARS shows either job submitted successfully or not.
   
Monitoring P-AIRCARS
~~~~~~~~~~~~~~~~~~~~
1. If you setup e-mail ids, you will get notifications of P-AIRCARS progress.

2. If you setup ``remotelogger``, you can monitor the P-AIRCARS progress in ``remotelogger``.

3. To use ``prefect`` dashboard, first create ``ssh`` tunnel to ``prefect`` port. The instructions can be seen as

.. code-block:: bash

   init-paircars-prefect status
   
Then from local machine, connect ``ssh`` tunnel

.. code-block:: bash

   ssh -N -L <prefect_port>:localhost:<prefect_port> <username>@<remote.cluster.name>
   
.. note::
  
   Ensure ``prefect_port`` is not being used in local machine.
   
Launch in web-browser: ``http://localhost:<prefect_port>/dashboard``
