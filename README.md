<p align="center">
  <img src="https://raw.githubusercontent.com/devojyoti96/P-AIRCARS/refs/heads/master/dark_logo.png" alt="P-AIRCARS Logo" width="200"/>
</p>
<p align="center">
  <h1>P-AIRCARS</h1> An automated spectropolarimetric calibration and imaging pipeline designed for solar radio observations using the <strong>Murchision Widefield Array (MWA)</strong> radio telescope. It performs end-to-end calibration, flagging, and imaging with a focus on dynamic solar data, supporting both spectral and temporal flexibility in imaging products.
</p>

## Background

<!-- start elevator-pitch -->

Solar radio data presents unique challenges due to the high variability and brightness of the Sun, as well as the need for high time-frequency resolution. The **P-AIRCARS** pipeline addresses these challenges by:

- Automating the calibration of interferometric data, including flux, phase, and polarization calibrations
- Supporting time-sliced and frequency-sliced imaging workflows
- Leveraging Dask for scalable parallel processing
- Providing hooks for integration with contextual data from other wavelegths for enhanced solar analysis

<!-- end elevator-pitch -->

## Documentation

P-AIRCARS documentation is available at: [p-aircars.readthedocs.io]

[p-aircars.readthedocs.io]: https://p-aircars.readthedocs.io 

## Software environment

**P-AIRCARS** is tested on Ubuntu22, Ubuntu 24, and CentOS 7 with Python 3.10. P-AIRCARS may not work in other operating system and python versions. If user wants to use P-AIRCARS in other environments, limited support is available in debugging or solving the issues. User may look at **Containersed Use** section in the docuement for these scenarios.

## Quickstart

<!-- start quickstart -->

**P-AIRCARS** is distributed on [PyPI]. To use it, install it in isolated conda environment. If **conda** is not installed in your system, see document for Conda installation instructions.

1. Set some environment variable

    ```text
     export PYTHONNOUSERSITE=1
      
     unset PYTHONPATH  
    ```

2. Create conda environment with python 3.10 with compaitable C/C++ libraries

    ```text
    conda create -n paircars_env --override-channels -c conda-forge python=3.10 gcc_linux-64=14 gxx_linux-64=14 gfortran_linux-64=14 cmake pkg-config pip
    
    conda activate paircars_env
    ```
    
    We suggest using **Mamba** for fast conda installtion and environment creation.

3. Install P-AIRCARS in conda environment

   ```text
   pip install paircars
   ```

4. Initiate necessary metadata 

    ```text
    init-paircars-setup --init
    ```
    By default, the necessary data will be saved in home directory and requires about 20 GB of disk space. We suggest using any other location with larger disk space and specify that by ``--datadir </full/path/to/paircars_datadir>`` in the above command.
    
5. Before running the pipeline, setup your data as following:
    
    -- Create a <target_datadir> and put all coarse channel measurement sets of solar scan of a single observation ID (OBSID) inside it.
    
    -- Create a <cal_datadir> and put all coarse channel measurement sets for calibrator observation of a single OBSID inside it.
    
6. Run P-AIRCARS pipeline

    ```text
    run-mwa-paircars <full path of target measurement set directory> --cal_datadir <full path of calibrator measurement set directory> --workdir <full path of work directory> --outdir <full path of output products directory>
    ```    
    
    N.B.: Always provide the entire direcotry path. Short path or only directory name may cause errors. Keep target measurement sets for a single OBSID and calibrator measurement sets for a single OBSID must be kept in seperate directories. If calibrator is not present, do not provide these information.

That's all. You started P-AIRCARS pipeline for analysing your MWA solar observation 🎉.

7. To see all running P-AIRCARS jobs

    ```text
    show-paircars-status --show
    ```
    
8. If P-AIRCARS is running in a local machine, see local log of any job using the <jobid>

   ```text
   run-mwa-mwalogger --jobid <jobid>
   ```
   
   N.B.: If you are running P-AIRCARS is cluster environment, first checkout **HPC Settings** in the document for viewing P-AIRCARS log remotely using prefect dashboard.

   
9. Output products will be saved in : `<path of output products directory>`

[pypi]: https://pypi.org/project/paircars/

<!-- end quickstart -->

## Sample dataset
User can download and test entire P-AIRCARS pipeline using the sample dataset available in Zenodo: https://doi.org/10.5281/zenodo.18641232. Do not use this sample dataset for any publication without permission from the developer.


## Acknowledgements

P-AIRCARS is developed by Devojyoti Kansabanik (NCRA-TIFR, Pune, India and CPAESS-UCAR, Boulder, USA) and an incarnation of [AIRCARS][aircars]. Other contributors are, Surajt Mondal (NCRA-TIFR, Pune, India), Soham Dey (NCRA-TIFR, Pune, India), and Puja Majee (NCRA-TIFR, Pune, India). If you use **P-AIRCARS** for analysing your MWA solar observations, include the following statement in your paper, and cite the following papers:

[aircars]: https://github.com/devojyoti96/AIRCARS 
```text
This MWA solar observations are analysed using P-AIRCARS pipeline. 
```

1. Cite P-AIRCARS software in zenodo: [https://doi.org/10.5281/zenodo.18625477][kansbanikzenodo]

[kansabanikzenodo]: https://doi.org/10.5281/zenodo.18625477

2. [Kansabanik et al., 2025, ApJS, v278:26][kansabanik2025]

[kansabanik2025]: https://doi.org/10.3847/1538-4365/adc443

3. [Kansabanik et al., 2023, ApJS, v264:47][kansabanik2023]

[kansabanik2023]: https://doi.org/10.3847/1538-4365/acac79

4. [Kansabanik et al., 2022, ApJ, v932:110][kansabanik2022a]

[Kansabanik2022a]: https://doi.org/10.3847/1538-4357/ac6758

5. [Kansbanik 2022, Solar Physics, v297:122][kansabanik2022b]

[kansabanik2022b]: https://doi.org/10.1007/s11207-022-02053-x

6. [Mondal et al., 2019, ApJ, v875:97][mondal2019]

[mondal2019]: https://doi.org/10.3847/1538-4357/ab0a01

If you use observations before 2015, include this additonal statement and citation:

```text
Flux calibration of the observations are done using the menthod described in the following paper.
```

7. [Kansabanik et al., 2022, ApJ, v927:17][kansabanik2022c]

[kansabanik2022c]: https://doi.org/10.3847/1538-4357/ac4bba 

P-AIRCARS name is given by Dr. Barnali Das (NCRA-TIFR, Pune, India)

## License

This project is licensed under the MIT License.
