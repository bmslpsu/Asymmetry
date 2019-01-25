# Asymmetry

Analysis functions for Asymmetry project

## Intructions
1. Download folders "Experiment_Asymmetry" & "Experiment_Asymmetry_Control" from box.com
2. Edit main function's root folders to your folder location "...\Experiment_Asymmetry" and "...\Experiment_Asymmetry_Control"
3. Add all functions to path and run main functions

## File Structure
* In the folder "Experiment_Asymmetry" DAQ files are stored in the root, video files in "Vid", and angle files in "Angle".
* In "Experiment_Asymmetry_Control" there are folders for high constrast, low contrast, & interpolated motion. Each of these folders includes subfolders labeled by spatial frequency (0 = Random). There are no video or angle files associated with these experiments.

## Analysis Functions
### Main

* ```Analyze_Asymmetry.m```

	*Analysis function for 45 deg/s random spatial frequency >>> WBA plots for CW & CCW motion (video and DAQ)*
	
* ```Analyze_Asymmetry_Control.m```

	*Analysis function for 30-150 deg/s and all spatial frequencys >>> WBA plots for CW & CCW motion (DAQ only)*
	
* ```MakeSpatFreqFig.m```

	*Runs ```Analyze_Asymmetry_Control.m``` and compiles results into figure*
	
### Supplementary

* ```PlotPatch.m```

	*Plots means with STD patch*

* ```hampel.m```

	*Hampel filter*

* ```scatplot.m```

	*Scatter plot with densities*
	
## Experiment Functions

* ```Experiment_Asymmetry.m```

	* Runs asymmetry experiment using controller V3 with LED panels, NiDAQ, Basler acA640-120gm camera *
