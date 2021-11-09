Robustness and Quality
======================

In this section two major input parameters for P-AIRCARS, *Quality* and *Robustness*, will be discussed.


Quality (Q)
***********
*Quality* parameters determines the quality of the final images. With the increase of *Quality* value, antennas will be added in smaller steps, self-calibration will be continued for smaller value of residual fluxes and final CLEANing of the images will be deeper. *Quality* takes three values; 0, 1, 2. Default value is 1.

Robustness (R)
**************
*Robustness* parameter determines the robustness of the calibration and imaging. Increasing the robustness parameter will use strong convergence criteria for self-calibration, uses higher value of minimum SNR for gain solutions, performs calibration at smaller time and frequency intervals. *Robustness* takes three values; 0, 1, 2. Default value is 1.

The calibration and imaging parameters are determined based on the combination of *Quality* and *Robustness*. The calibration and imaging parameters for each of the combinations with increasing significance are listed below. With the increasing significance, the calibration and imaging quality will be better, but the computational time will also increase. Default value of 1 for *Quality* and *Robustness* is an optimization between final image quality and computational time.

Combinations
************
  
**Q = 0, R = 0**

residual_frac	=	0.03

start_sigma		=	9.0

min_sigma		=	6

gain_minsnr		=	3.0		

DR_delta_rms	=	25

DR_delta_neg	=	20

min_DR			=	20

max_DR			=	100

skip_time		=	480 s

skip_freq		=	2560 kHz

**Q = 0, R = 1**

residual_frac	=	0.03

start_sigma		=	9.0

min_sigma		=	7

gain_minsnr		=	3.0		

DR_delta_rms	=	22

DR_delta_neg	=	18

min_DR			=	25

max_DR			=	500

skip_time		=	240 s

skip_freq		=	2560 kHz

**Q = 0, R = 2**

residual_frac	=	0.03

start_sigma		=	9.0

min_sigma		=	8

gain_minsnr		=	3.0		

DR_delta_rms	=	20

DR_delta_neg	=	15

min_DR			=	30

max_DR			=	1000

skip_time		=	120 s

skip_freq		=	2560 kHz


**Q = 1, R = 0**

residual_frac	=	0.015

start_sigma		=	10.0

min_sigma		=	7

gain_minsnr		=	4		

DR_delta_rms	=	20

DR_delta_neg	=	15

min_DR			=	30

max_DR			=	1000

skip_time		=	120 s

skip_freq		=	1280 kHz

**Q = 1, R = 1 (Default)**

residual_frac	=	0.015

start_sigma		=	10.0

min_sigma		=	8

gain_minsnr		=	4		

DR_delta_rms	=	18

DR_delta_neg	=	12

min_DR			=	35

max_DR			=	5000

skip_time		=	60 s

skip_freq		=	1280 kHz

**Q = 1, R = 2**

residual_frac	=	0.015

start_sigma		=	10.0

min_sigma		=	9

gain_minsnr		=	4		

DR_delta_rms	=	15

DR_delta_neg	=	10

min_DR			=	40

max_DR			=	10000

skip_time		=	30 s

skip_freq		=	1280 kHz

**Q = 2, R = 0**

residual_frac	=	0.01

start_sigma		=	11.0

min_sigma		=	8

gain_minsnr		=	4.5		

DR_delta_rms	=	18

DR_delta_neg	=	12

min_DR			=	40

max_DR			=	10000

skip_time		=	15 s

skip_freq		=	640 kHz

**Q = 2, R = 1**

residual_frac   =   0.01

start_sigma		=	11.0

min_sigma		=	9

gain_minsnr		=	4.5		

DR_delta_rms	=	15

DR_delta_neg	=	10

min_DR			=	45

max_DR			=	50000

skip_time		=	12 s

skip_freq		=	640 kHz

**Q = 2, R = 2**

residual_frac   =  0.01

start_sigma		=	11.0

min_sigma		=	10

gain_minsnr		=	4.5		

DR_delta_rms	=	12

DR_delta_neg	=	8

min_DR			=	50

max_DR			=	100000

skip_time		=	10 s

skip_freq		=	640 kHz




















































