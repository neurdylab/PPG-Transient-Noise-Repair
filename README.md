<!-- ABOUT THE PROJECT -->
## About The Project

This project, coded in MATLAB, aims to detect and interpolate over transient noises in PPG waveforms. The project loads a directory containing PPG waveform files, preprocesses them, then uses the collected information to reinterpolate over noisy regions.

The noise detection is performed using information including peak amplitude, inter-beat intervals, etc., and the interpolation over noisy region is based on information of neighboring beats.

Note that the noise interpolation is only performed on the bandpass-filtered, not the original, version of the waveform.


<!-- RUNNING THE CODE -->
## Running The Code

* After downloading the repository, open run_fix.m in MATLAB.
* In the first line of code in run_fix.m, replace the '...' in brackets with the absolute directory directly containing the PPG files.
* The preprocessed files and corrected files will be contained in two separate folders in the directory.
* If the sampling rate is not 2000Hz, replace the fs_phys parameter in preproc_phsio.m