
1. System requirements:
	OS: the code was tested on Windows 11 Pro, RAM>=32 GB
	MATLAB: MATLAB R2024a
	Required MATLAB package:
		Aerospace Toolbox v24.1
		Curve Fitting Toolbox v24.1
		Financial Toolbox v24.1
		Mapping Toolbox v24.1
		Optimization Toolbox v24.1
		Parallel Computing Toolbox v24.1
		Statistics and Machine Learning Toolbox v24.1

2. Instruction on the use of the code and how to run the demo:
	2.1 Open MATLAB R2024a and set the working directory to the path where "GWKarst" is stored;
	2.2 Open the "RUN_main_script.m" file and run the code;
	2.3 The model will be run and the progress will be shown in the command window.

	Expected output: the result of the demo model simulation will be stored in the struct variable "result". It includes the timeseries of streamflow simulation at the predefined locations and the timeseries of head simulations of the karst grids (here is the head for the karst matrix).

	Expected run time: it took about 2 mins for the demo case on the computer with 13th Gen Intel(R) Core(TM) i7-13700   2.10 GHz.
	
