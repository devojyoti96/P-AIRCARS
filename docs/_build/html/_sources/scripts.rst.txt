paircars scripts
================

control_paircars
----------------
Usage:  PAIRCARS master controller for each day calibration

Options: -h, \--help            		show this help message and exit

		\--basedir=Directory path, 	Name of base directory for a given day

run_paircars
------------
Usage:  Perform self calibration of a single time and frequency slice

Options: -h, \--help            	show this help message and exit

		\--msname=Measurement Set,	Name of measurement set of a single time and frequency slice

		\--metafits=Metafits file,	Name of metafits file of the observation

		\--basedir=Directory path,	Name of the base directory

		\--workdir=Directory path,	Name of the working directory

		\--ref_freq_avg=Float,  	Frequency averaging for reference ms

		\--ref_time_avg=Float, 		Time averaging for reference ms

		\--ref_time_freq=Boolean,	Reference measurement set or not

		\--do_bandpass=Boolean,		Perform bandpass calibration or not

		\--do_polcal=Boolean,  	 	Perform polarisation calibration or not

		\--cal_attenuation=Float,	Attenuation in dB for calibrator observation

		\--caltables=String, comma separated,	Previous calibration tables

		\--scratch=Boolean,    		Start from scratch or not for reference time frequency slice

		\--wsclean=Boolean,    		Use WSClean for imaging or not

		\--cpu_frac=Float,     		Fraction of cpu to use

run_intensity_selfcal
---------------------
Usage:  Perform intensity self calibration of a single time and frequency slice

Options: -h, \--help            		show this help message and exit

		\--msname=Measurement Set,		Name of measurement set of a single time anf frequen slice

		\--metafits=Metafits file,		Name of metafits file of the observation

		\--workdir=Directory path,		Name of the working directory
	
		\--dopoint=Boolean,    			Want to try with point source model

		\--verbose=Boolean,    			Verbose mode

		\--interactive=Boolean,			Interactive mode

		\--fresh=Boolean,      			Start fresh self calibration loop

		\--reduce\_flags=Boolean,		Try to reduce flag solutions if it is more than 5%

		\--scratch=Boolean,    			Start from scratch or not for reference time frequency slice

		\--caltables=String, comma separated, Previous caltables

		\--wsclean=Boolean,    			Use WSClean for imaging or not

run_bandpass_selfcal
--------------------
Usage:  Perform bandpass self calibration

Options:  -h, \--help           		show this help message and exit

		  \--msname=Measurement Set,	Name of measurement set of a single time and frequency slice

		  \--metafits=Metafits file,	Name of metafits file of the observation

		  \--workdir=Directory path,	Name of the working directory

		  \--verbose=Boolean,    		Verbose mode

		  \--interactive=Boolean,		Interactive mode

		  \--fresh=Boolean,       		Start fresh self calibration loop

		  \--caltables=String, comma separated,	Previous caltables

		  \--wsclean=Boolean,     		Use WSClean for imaging or not

run_pol_selfcal
---------------
Usage:  Perform polarisation self calibration of a single time and frequency slice

Options:  -h, \--help					show this help message and exit

		  \--msname=Measurement Set,	Name of measurement set of a single time anf frequency slice

		  \--metafits=Metafits file,	Name of metafits file

		  \--workdir=Directory path,	Name of the working directory
	
		  \--verbose=Boolean,     		Verbose mode

		  \--interactive=Boolean,		Interactive mode

		  \--fresh=Boolean,       		Start fresh self calibration loop
		
		  \--gaincal=Boolean,     		Perform gaincal using leakage corrected model (Only do when no calibrator observation is present)

		  \--caltables=String, comma separated,	Previous caltables

		  \--wsclean=Boolean,     		Use WSClean for imaging or not








