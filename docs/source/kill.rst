Stop P-AIRCARS Job
==================

Since P-AIRCARS runs several seperate parallel processes in background, if user wants to stop the pipeline it is not straightforward. Use the following P-AIRCARS command line tool and ``<jobid>`` to stop a particular P-AIRCARS pipeline job. It will stop all child processes of that particular P-AIRCARS job.

.. code-block :: bash

    kill-paircars-job --jobid <jobid>
    
    


