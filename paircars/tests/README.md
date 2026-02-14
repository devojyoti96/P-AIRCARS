# Tests for P-AIRCARS

This directory contains the full testing suite for the **P-AIRCARS** pipeline using ``pytest``.

## Steps to test
Install ``pytest`` using ``pip install pytest`` before running the test.


1. Go to test directory
   ```
   text
   cd <repo_path>/paircars/tests

2. Download test data

    ```text
    python3 download_test_data.py
    ```
    
3. Run utils module test:

    ```
    text
    pytest -s -v utils
    ```
    
4. Run pipeline module test:
    
   ```
   text
   pytest -s -v pipeline
   ```
    




