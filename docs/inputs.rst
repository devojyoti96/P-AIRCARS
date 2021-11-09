Details of Inputs
=================
In this section the details of the P-AIRCARS inputs are given.

Inputs
******

1. **Data Directory** : Name of the directory where all measurement sets to be imaged are present. If user want to run in *offline* mode, please be sure that all metafits files correcponding to all measurement sets are also present in the *Data Directory*. This is a mandatory input.

2. **Base Directory** : Name of the directory where all calibration and imaging analysis will be done. Either give an existing directory, otherwise P-AIRCARS will create the *Base Directory* in the given path. This is a mandatory input.

3. **Image Directory** : Name of the directory to save final images. If no input is given final images will be saved in a sub-directory named 'final_images' inside the *Base Directory*.

4. **Time range** : Time ranges to be imaged. Time range format is hh0:mm0:ss0.ff0~hh1:mm1:ss1.ff1,hh2:mm2:ss2.ff2~hh3:mm3:ss3.ff3,...

5. **Channel range** : Channel ranges to be imaged. Channel range format is ch0~ch1,ch2~ch3,....

6. **Caltable Names** : Prior caltable names obtained from calibrator observations. Give the full paths separated by commas.

7. **Robustness** : Determine the robustness of the calibration solutions and sef-calibration convergence. It takes three values 0, 1, 2 with increasing robustness. Default value is 1.

8. **Quality** : Determine the image quality. It takes three values 0, 1, 2 with increasing quality of the images. Default value is 1.

9. **Verbose** : If selected, verbose output logs will be saved. Also all the images, measurement sets and calibration tables are also saved. In verbose mode amount of disk space requirements will be more. Turn on only for debuggung purposes. Default it is turned off.

10. **Save Logs** : If selected, logs for self-calibration and imaging rounds will be saved. Default it is turned off.

11. **Interactive** : By default it is turned off. If selected, user can do interactive cleaning and also can change inputs manually during runtime. Only useful for specific debugging.

12. **Decorrelation Correction** : If selected, perform the amplitude decorrelatation correction. By default it is switched on. Details of the correction can be found in Mondal et al. 2020 https://doi.org/10.3847/1538-4357/ab0a01

13. **Reference Antenna** : Select the reference antenna from first 30 core antennas. Default is 1.

14. **Auto-calculate Parameters** : If selected, P-AIRCARS will calculate calibration and imaging parameters automatically from the data. By default, it is switched on. When it is turned off, user can give their own calibration and imaging parameters in advanced inputs. Only experienced users should turn this off.

15. **Use WSClean** :  Use WSClean for imaging. By default it is switched on. If it is turned off, CASA will be used for imaging. Imaging using CASA is significantly slow compared to WSClean. Thus users are instructed to use WSClean. If WSClean is not installed, P-AIRACARS will automatically use CASA for imaging even if Use WSClean is turned on.

16. **CASA Auto-masking** : When CASA is being used for imaging, CASA auto-masking can be used by enabling this option. By default it is not activated.

17. **Notification** : If selected, P-AIRCARS will automatically send the updates about the progress of calibration and imaging to the user specific e-mail id.

18. **CASA Mask** : When using CASA for imaging, user may supply a CASA mask file here. By default it is not activated. It will be activated when Use WSClean is turned off.

19. **CASA Mask String** : CASA soecific mask string when using CASA for imaging. By default it is not activated.  It will be activated when Use WSClean is turned off.

20. **e-mail** : User e-mail address to send notifications. If no e-mail id is given, no notification can be sent.

21. **Bandpass** : Whether perform bandpass self-calibration or not. By default it is turned on.

22. **Polcal** : Whether perform polarisation calibration and imaging. By default it is turned on. When it off, P-AIRCARS will only make Stokes-I images.

23. **Frequency Interval** : Frequency interval to make images. By default it is 160 kHz.

24. **Time Interval** : Time interval to make images. By default it is 0.5 s.

25. **Frequency Width** : Bandwidth to average. By default it is 160 kHz.

26. **Time Width** : Temporal averaging. By default it is 0.5 s.

27. **CPU(%)** : Percentage of total available cpus to be used. By default it is 50\%.

28. **Clear Virtual SCreens** : Clear the virtual screens opened during the P-AIRCARS run.

29. **Save models** : If selected, save final model images. By default it is turned off.

30. **Save residuals** : If selected, save final residual images. By default it is turned off.

31. **XY cutout** : Size of the final image cutout. Default size is 3 degree.

32. **Perform Flagging** : Perform flagging during the calibration or not. By default it is on.

33. **Use aNKflag** : If selected, use aNKflag for flagging. If aNKflag is not installed, do not use.

34. **HPC environment** : Switch on if running is high performance computing (HPC) environment. When it is on fill the settings for HPC environment. 

35. **Input file** : Load the P-AIRCARS inputs from previous saved P-AIRCARS input file. Extension of P-AIRCARS input file is '.paircars'. After loading inputs will be filled automatically.

36. **Restart P-AIRCARS** : Restart P-AIRCARS from the interupted step. It is off by default.

37. **Fresh start P-AIRCARS** : Start P-AIRCARS in fresh. If selected, all previous analysis in the same *Base Directory* will be deleted. By default it is switched on.


Advanced Inputs
***************
Advanced input block will be activated when *Auto-calculate Parameters* is switched off. Use *Advanced Inputs* carefully and suggest to use only by experienced users.

1. **Pixelsize** : Pixel size of the image in arcsec or arcmin, e.g. '20.0arcsec' or '1.0arcmin'.   

2. **Number of pixels** : Total number of pixels in the image. Product of *Pixelsize* and *Number of pixels* determine the size of the field-of-view.

3. **Multiscale scales** : Multiscale scales in units of number of pixels.

4. **UV-taper** : UV taper string in CASA format

5. **Start Sigma** : Starting value of sigma for rms based thresholding in self-calibration.

6. **Sigma Step** : Stepsize to reduce the sigma value with self-calibration iterations.

7. **Minimum Sigma** : Minimum value of sigma to reduce during self-calibration.

8. **Residual Flux Fraction** : Fraction of flux remained in the residual image to stop the self-calibration.

9. **UV range** : UV range in CASA format for self-calibration.

10. **Skip Frequency** :  Frequency step to skip during calibration. Calibration will be performed after *Skip Frequency* steps.

11. **Skip Time** : Time step to skip during calibration. Calibration will be performed after *Skip Time* steps.

12. **Gain Minimum SNR** : Minimum SNR of the gain solutions to be used. Gain solutions less than minimum SNR will be rejected.

13. **DR RMS step** : If rms based dynamic range increases less than this value, self-calibration converged.

14. **DR Negative step** : If negative based dynamic range increases less than this value, self-calibration converged.

15. **Minimum DR** : Minimum dynamic range to save a final image.

16. **Maximum DR** : Maximum dynamic range required. If the imaging dynamic range reached this value, self-calibration stopped.

17. **Minimum Selfcal SNR** : Minimum antenna based SNR required to allow the self-calibration to start.

18. **Extra time averaging** : If large calibration solutions flagged, average over *Extra time* to gain SNR. The value is in second.

19. **Maximum  average time** : Maximum allowed value of the *Extra time*.

High Performance Computing (HPC) settings
*****************************************






