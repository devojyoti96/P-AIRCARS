Directory Structure and Data Products
=====================================
Once user started a P-AIRCARS pipeline job, P-AIRCARS assign a unique JobID based on current time in millisecond precision in YYYYMMDDHHMMSSmmm format.
Note down the Job ID to view the logger.

The following output will appear in terminal:

.. code-block :: console

    ########################################
    Starting P-AIRCARS Pipeline....
    #########################################

    ###########################
    P-AIRCARS Job ID: <YYYYMMDDHHMMSSmmm>
    Work directory: <workdir>
    Final product directory: <outputdir>
    ###########################

Directory structure
-------------------

All intermediate data products will be saved in ``<workdir>/<target_obsid>_<jobid>``.
    
.. note::
   
   ``jobid`` is added to work directory name to avoid conflicts.
   

All final data products will be saved in ``<outputdir>/<target_obsid>_target``.

All calibrator data products will be saved in ``<outdir>/<calibrator_obsid>_cal``.

.. note :: 

   In local workstations, it is okay to choose the same ``<workdir>`` and ``<outputdir>``. In HPC environment, generally, high-speed disks are used during data-processing, which may have limited storage life-time, and has seperate long-term storage disks. It is recommended to choose ``<workdir>`` path inside the high-speed disk and ``<outputdir>`` inside the long-term storage disk. Otherwise, there may be possiblity that final data-products will be removed after certain time. 

.. note::
   
   Inside ``<workdir>`` and ``<outputdir>``, another directory with target observation ID, ``<targrt_obsid>`` will be created. These will be the ``<workdir>`` and ``<outputdir>`` for that OBSID. This ensures if different OBSID targets are provided same ``<workdir>`` and ``<outputdir>``, there will be no mixup. Some of these follwing directories will be only present when back are kept.

.. admonition:: Click here to see directory structure in work directory
   :class: dropdown
   
   .. mermaid::

       graph LR
           WD["Work directory:<br>workdir"] --> CAL["Calibrator ms:<br>calibrator.ms"]
           WD["Work directory:<br>workdir"] --> CAL["Bandpass calibration tables:<br>calibrator*.bcal"]
           WD["Work directory:<br>workdir"] --> CAL["Crossphase calibration tables:<br>calibrator*.kcrosscal"]
           WD --> SCMS["`Self-cal ms(s):<br>selfcal_*_spw_*.ms`"]
           WD --> INTSCDIR["`Intensity self-cal directories:<br>selfcal*_spw_*_selfcal.int`"]
           WD --> POLSCDIR["`Polarisation self-cal directories:<br>selfcal*_spw_*_selfcal.pol`"]
           WD --> TMS["`Target ms(s):<br>targets*_spw_*.ms`"]
           WD --> LOG["Log directory:<br>logs"]
           LOG --> LOGF["*.log"]
                      
.. admonition:: Click here to see directory structure in output directory
   :class: dropdown
   
   .. mermaid::

       graph LR
           WD["Output directory:<br>{outdir}"] --> WDC["Calibrator output directory:{outdir}/{cal_obsid}_cal"]
           WD["Output directory:<br>{outdir}"] --> WDT["Target output directory:{outdir}/{target_obsid}_target"]
           WDC --> CALTABLE["Caltable directory:<br>caltables"]
           WDC --> DPC["`Diagnostic plots:<br>diagnostic_plots`"]
           WDC --> FS["`Flag summary:<br>flag_summary`"]
           WDC --> FV["`Flag backup:<br>ms_flags`"]
           DPC --> DPPDF["`Diagnostic plots of calibrator ms and caltables in PDF:<br>*.pdf`"]
           CALTABLE --> BCAL["`Bandpass caltable:<br>*.bcal`"]
           CALTABLE --> BCALATT["`Attenuation corrected bandpass caltable:<br>*.bcal.att`"]
           CALTABLE --> KCROSSCAL["`Crossphase caltable:<br>*.kcrosscal`"]
                     
           WDT --> DS["`Dynamic spectra:<br>dynamic_spectra`"]
           WDT --> DS["`Calibrated visibilities:<br>calibrated_ms`"]
           WDT --> DPT["`Diagnostic plots:<br>diagnostic_plots`"]
           WDT --> SELFCALTABLE["Caltable directory:<br>caltables"]
           WDT --> IMG["`Image directory:<br>imagedir_f_*_t_*_w_briggs_*`"]
           IMG --> IMAGE["Fits image:<br>images"]
           IMG --> MODEL["Fits models:<br>models"]
           IMG --> RES["Fits residual:<br>residuals"]
           IMG --> PBCOR["`Primary beam<br>corrected images:<br>pbcor_images`"]
           IMG --> TBIMG["`Brightness temperature images:<br>tb_images`"]
           IMG --> OVRPNG["Overlays of EUV:<br>PNG format:<br>overlays_pngs"]
           
           
Data products
-------------
Pipeline produces calibrated visibilities as well as several imaging products.

Dynamic spectrum
~~~~~~~~~~~~~~~~
Dynamic spectra for all (or the ones selected) target scans are available in ``dynamic_spectra`` directory inside the output directory ``<outputdir>/<target_obsid>``.

Diagnostic plots
~~~~~~~~~~~~~~~~
Diagnostic plots for all measurement sets and calibration tables in pdf format in ``diagnostic_plots`` directory inside the output directory ``<outputdir>/<target_obsid>``.

Flag summary
~~~~~~~~~~~~
Flag summary files for each measurement sets

Measurement set flags
~~~~~~~~~~~~~~~~~~~~~
Flags of the final calibrated measurement sets are saved in ``ms_flags`` directory inside the ``<outputdir>/<target_obsid>``. These flags can be used later to restore and re-image calibrated measurement sets. s

Calibrated visibilities
~~~~~~~~~~~~~~~~~~~~~~~
By default calibrated measurement sets will kept in ``<outdir>/<target_obsid>/calibrated_ms`` directory after final imaging with the naming format ``target_<target_obsid>_ch<coarse_chan>_spw_<chanrange>.ms`. If ``no_calibrated_ms`` parameter is turned on during processing, they will not be kept.

Imaging products 
~~~~~~~~~~~~~~~~
Imaging products are available in ``imagedir_f_<freqres>_t_<timeres>_w_<weight>_<robust>`` directory inside output directory ``<outputdir>/<target_obsid>``. If imaging is performed with different time and frequency resolutions or different weighting schemes, seperate image directories with corresponding parameters will have the corresponding images. 

1. **Image fits in RA/DEC** - Fits images in RA/DEC coordinate are available in ``imagedir_f_<freqres>_t_<timeres>_w_<weight>_<robust>/images`` directory inside work directory. These are not primary beam corrected.

.. note ::
    
   All fits images have some P-AIRCARS specific metadata in the header and some image statistics.
   
.. admonition:: Click here to see details of these metadata
   :class: dropdown
                       
    PIPELINE= 'P-AIRCARS' # Pipeline name
    
    AUTHOR  = 'DevojyotiKansabanik' # Pipeline developer         
                                                           
    MAX     =  ``<maxval>`` # Maximum value on the solar disc       
                                               
    MIN     =  ``<minval>`` # Minimum value on the solar disc      
                                               
    RMS     =  ``<rms>`` # RMS value outside solar disc       
                                            
    SUM     =  ``<sum>`` # Total sum on the solar disc  
                                            
    MEAN    =  ``<mean>`` # Mean value on the solar disc   
                                                  
    MEDIAN  =  ``<median>`` # Median value on the solar disc  
                                                 
    RMSDYN  =  ``<rmsdyn>`` # RMS based dynamic range, ``<maxval/rms>``
                                                    
    MIMADYN =  ``<minmaxdyn>`` # Min-max based dynamic range, ``<maxval/abs(minval)>``  
    
    CALAPP  =   ``TRUE/FALSE`` # Whether calibrator solutions were applied or not
    
    POLSELF =   ``TRUE/FALSE`` # Whether polarisation selfcal is performed or not
    
    LEAKCOR =   ``TRUE/FALSE`` # Whether polarisation leakage correction is made or not
    
    LEAKUNIT =   ``PERCENT`` # Residual leakage unit if LEAKCOR is TRUE
    
    QLEAK    =    ``<qleak>`` # Residual Stokes I to Stokes Q leakage if LEAKCOR is TRUE
    
    ULEAK    =    ``<uleak>`` # Residual Stokes I to Stokes U leakage if LEAKCOR is TRUE 
    
    VLEAK    =    ``<vleak>`` # Residual Stokes I to Stokes V leakage if LEAKCOR is TRUE 
        
.. note::

    If CALAPP is False, flux calibration is performed using the method described in `Kansabanik et al., 2022, ApJ, v927:17 <https://doi.org/10.3847/1538-4357/ac4bba>`_. In this case, flux density and spectrum may not be very accurate. We generally expect observations before middle of 2014 have this issue. These observations may also have unreliable pol-rotation (Stokes Q, U, V mixing).
    
    If POLSELF is False, polarisation self-calibration could not be performed, and hence only image-based leakage correction is done. If this case, very low-level polarised sources should be not be considered as true source. 
 
2. **Primary beam corrected image fits** - Primary beam corrected fits images are available in ``imagedir_f_<freqres>_t_<timeres>_w_<weight>_<robust>/pbcor_images`` directory inside work directory.

3. **Brightness temperature image fits** - Brightness temperature fits images are available in ``imagedir_f_<freqres>_t_<timeres>_w_<weight>_<robust>/tb_images`` directory inside work directory. 

4. **CLEAN model and residual fits** - CLEAN model and residual images are available in ``imagedir_f_<freqres>_t_<timeres>_w_<weight>_<robust>/models`` and ``imagedir_f_<freqres>_t_<timeres>_w_<weight>_<robust>/residuals`` directory inside work directory. These are not saved only if ``keep_backup`` option is switched on.

5. **Radio images in helioprojective coordinates** - Directory names ``hpcs`` inside directories like, ``images, pbcor_images, tb_images`` inside the image directory will have the FITS images in helioprojective coordinates. Images in PNG and PDF formats are also available in ``pngs`` and ``pdfs`` directories inside the parent directories.

.. note ::
  
   Header of helioprojective maps have wavelength information in unit of ``centimeter`` or ``meter``.  

6. **Overlays on EUV images** - Overlays on EUV images are saved in PNG formats in ``imagedir_f_<freqres>_t_<timeres>_w_<weight>_<robust>/overlays_pngs``. By default, images at 30s time interval per corase channel will be overlayed, unless requested to make overlays for all images. 



