# 22q_longitudinal
Analysis of age-related trajectories for thalamocortical functional network connectivity in 22qDel and typical controls. 

## Overview
*[2q_longitudinal_gamm_final.md](https://github.com/charles-schleifer/22q_tcc_longitudinal/blob/main/22q_longitudinal_gamm_final.md): R markdown knitted output with code and plots for main analysis. 
* analysis steps include:
  * preprocessing and thalamocortical functional connectivity calculation (see: https://github.com/charles-schleifer/22q_hoffman)
  * data harmonization (longComBat)
  * generalized additive mixed models (GAMMs)
  * plotting connectivity vs age curves
  * rate of change
  * group differences between curves
  * demographics table 

## Dependencies
* requires wb_command for ciftiTools functions that read/plot MRI data. 
  * download: https://www.humanconnectome.org/software/get-connectome-workbench
  * script expects the workbench directory to be `/Applications/workbench/bin_macosx64` (either download to this location, or edit the path in the script)
* to read fMRI results from the hoffman2 server, the script expects that the server `hoffman2.idre.ucla.edu:/u/project/cbearden/data` is mounted to your local machine at the path `~/Desktop/hoffman_mount` using an application such as SSHFS (mac download: https://osxfuse.github.io/)
  * requires first-level MRI results to be already computed on server (see https://github.com/charles-schleifer/22q_hoffman)


 
