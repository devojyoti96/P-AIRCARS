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

P-AIRCARS documentation is available at: [paircars.readthedocs.io]

[paircars.readthedocs.io]: https://p-aircars.readthedocs.io 

## Quickstart

<!-- start quickstart -->

**P-AIRCARS** is distributed on [PyPI]. To use it:

1. Create conda environment with python 3.10

    ```text
    conda create -n paircars_env python=3.10
    conda activate paircars_env
    ```

2. Install P-AIRCARS in conda environment

   ```text
   pip install paircars
   ```

3. Initiate necessary metadata and prefect server

    ```text
    init-paircars-setup --init --prefect_server
    ```
    
4. Before running the pipeline, setup your data as following:
    
    -- Create a <target_datadir> and put all coarse channel measurement sets of solar scan of a single observation ID (OBSID) inside it.
    
    -- Create a <cal_datadir> and put all coarse channel measurement sets for calibrator observation of a single OBSID inside it.
    
5. Run P-AIRCARS pipeline

    ```text
    run-mwa-paircars <path of target measurement set directory> <path of target metafits file> --cal_datadir <path of calibrator measurement set directory> --cal_metafits <path of calibrator metafits> --workdir <path of work directory> --outdir <path of output products directory>
    ```    
    
    N.B.: Keep target measurement sets for a single OBSID and calibrator measurement sets for a single OBSID must be kept in seperate directories. If calibrator is not present, do not provide these information.

That's all. You started P-AIRCARS pipeline for analysing your MWA solar observation 🎉.

6. To see all running P-AIRCARS jobs

    ```text
    show-paircars-status --show
    ```
    
7. To see prefect dashboard (only work if you started prefect server)

   ```text
   run-mwa-mwalogger
   ```
      
8. If you did not start prefect server, see local log of any job using the <jobid>

   ```text
   run-mwa-mwalogger --jobid <jobid>
   ```
   
9. Output products will be saved in : `<path of output products directory>`

[pypi]: https://pypi.org/project/paircars/

<!-- end quickstart -->

## Sample dataset
User can download and test entire P-AIRCARS pipeline using the sample dataset available in Zenodo: https://doi.org/10.5281/zenodo.18641232. Do not use this sample dataset for any publication without permission from the developer.


## Acknowledgements

P-AIRCARS is developed by Devojyoti Kansabanik (NCRA-TIFR, Pune, India and CPAESS-UCAR, Boulder, USA) and Surajt Mondal (NCRA-TIFR, Pune, India) and an incarnation of [AIRCARS][aircars]. If you use **P-AIRCARS** for analysing your MWA solar observations, include the following statement in your paper, and cite the following papers:

[aircars]: https://github.com/devojyoti96/AIRCARS 
```text
This MWA solar observations are analysed using P-AIRCARS pipeline. 
```

1. Cite P-AIRCARS software in zenodo: https://doi.org/10.5281/zenodo.18625477

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


## License

This project is licensed under the MIT License.
