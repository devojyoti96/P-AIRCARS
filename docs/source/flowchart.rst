P-AIRCARS Flowchart
=====================

The pipeline follows several steps. By default, all steps are executed. If the user desires, they can selectively disable any step. However, the pipeline will only function correctly if its internal logic is maintained.

For example, if the user disables splitting the target scans and no target directory exists beforehand, imaging will not be performed. However, if self-calibration is disabled, the pipeline will skip that step and proceed with applying basic calibrations and generating final images.

.. admonition:: Recommendation
   :class: tip

   It is recommended to go through the flowchart of the pipeline and understand it before modifying the pipeline keys.


