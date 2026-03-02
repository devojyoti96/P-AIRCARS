Developer Instruction
=====================
1. Developers are instructed to follow the similar format of the scripts.

2. Include new functions in proper modules. If suitable module is not found, include them in ``basic_utils``.

3. After including any function, add a ``pytest`` in appropriate testing scripts and test it. Without proper ``pytest`` new functions will not be included in the main branch and future release.

Instructions for testing
-------------------------
## Steps to test
Install P-AIRCARS in developer mode

.. code-block:: bash
    
   git clone git@github.com:devojyoti96/P-AIRCARS.git
   
   cd P-AIRCARS
    
   pip install .[dev] 


1. Go to test directory
   
   .. code-block:: bash
   
      cd </full/path/to/repository>/paircars/tests

2. Download test data

   .. code-block:: bash
   
      python3 download_test_data.py
   
    
3. Run utils module test:

   .. code-block:: bash
   
      pytest -s -v utils

    
4. Run pipeline module test:
    
   .. code-block:: bash
   
      pytest -s -v pipeline


