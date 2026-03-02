Advanced CLI
=============

Calibration related CLI
-----------------------
    
1. Flagging of calibrators, use ``run-mwa-flag`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-flag -h  
   
1. Simulate visibilities for calibrators, use ``run-mwa-import-model`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-import-model -h

3. Perform basic calibration, use ``run-mwa-basic-cal`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-basic-cal -h
   
4. Apply basic calibration solutions, use ``run-mwa-apply-basiccal`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-apply-basiccal -h
   
5. Split measurement set for self-calibration or final imaging, use ``run-mwa-split`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-split -h
   
6. Perform self-calibration, use ``run-mwa-selfcal`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-selfcal -h
   
7. Apply self-calibration solutions, use ``run-mwa-apply-selfcal`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-apply-selfcal -h
   
8. Move phasecenter to solar center, use ``run-mwa-movetosun`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-movetosun -h   
   
   
MWA specific CLI
----------------
1. Make frequency interpolated MWA primary beam, use ``run-mwa-beam-interpolate`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-beam-interpolate -h
   
Solar specific CLI
------------------

1. To correct sidereal motion of the Sun, if the Sun is not tracked by the correlator delay center, use ``run-mwa-solar-siderealcor`` . This is useful for observations where the Sun is in sidelobe of the telescope primary beam.

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-solar-siderealcor -h
   
2. Make dynamic spectra of solar scans, use ``run-mwa-makeds`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-makeds -h
 
3. Make overlays on EUV images, use ``run-mwa-euvoverlay`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-euvoverlay -h  
   
   
Imaging related CLI
-------------------
   
1. Perform spectro-polarimetric snapshot imaging, use ``run-mwa-imaging`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-imaging -h
   
2. Perform primary beam correction of MWA primary beam, for a single image, use ``run-mwa-singlepbcor`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-singlepbcor -h
   
3. Perform primary beam corrections of MWA primary beam for all images in a directory, use ``run-mwa-pbcor`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-pbcor -h
   
Ploting related CLI
--------------------
1. To make diagnostic plots of measurement sets, use ``run-mwa-msplot`` .

.. admonition:: Click here to see parameters
   :class: dropdown

   .. program-output:: run-mwa-msplot -h
   




