# 2017-Cu-TGS-Void-Swelling
Data repository for 2017 Scripta paper on detecting void swelling using transent grating spectroscopy

%Read me file for data and analysis associated with ``Detecting ion irradiation-induced void swelling in pure copper using transient grating spectroscopy" - Dennett et. al. 

Include in this repository are the following
- Raw, unprocessed SEM, TEM, and STEM images of lamella lifted out from irradiated samples
- Binary processed STEM images used for the determination of areal swelling versus depth
- Raw TGS data as collected
- Iterations of SRIM profiles used to determine the ion implantation energy to best match TGS depth
- Matlab processing scripts used to determine the dominant frequency of each TGS measurement, these frequencies can later be used to determine the acoustic wave speed from the calibrated grating spacing.
	- The function param_extract.m is used to find these dominant frequencies, instructions on its use can be found in the comments of the function
