Preprocessing
=============

``fMRwhy`` includes pipelines for basic anatomical and functional preprocessing.

As a first step, the T1-weighted anatomical image is coregistered to the template functional image (defined in the settings file)
using SPM12's coregister/estimate functionality, which maximizes normalised mutual information to generate a 12 degree-of-freedom transformation matrix.
Before resampling to the functional resolution, this coregistered T1-weighted image is segmented using tissue probability maps and SPM12's unified segmentation algorithm (Ashburner and Friston, 2005).
This yields subject-specific probability maps for gray matter, white matter, CSF, soft tissue, bone and air in the subject functional space.
All of these probability maps are then resampled (using SPM12's coregister/write) to the subject functional resolution.
Masks are then generated for grey matter, white matter, CSF, and the whole brain (which is a combination - logical OR after thresholding - of the previous three masks).
These are overlaid on the coregistered and resampled T1w image, to allow visual inspection of segmentation and registration quality in the QC report.

Anatomical regions of interest were then taken from the cytoarchitecture-based atlases in the SPM Anatomy Toolbox (Eickhoff et al., 2005).
For the motor cortex, regions 4a and 4p were used. For the amygdala, regions LB, IF, SF, MF, VTM, and CM were used. For the fusiform gyrus, regions FG1, FG2, FG3, and FG4 were used.
Regions of interest were transformed from MNI152 space to the subject functional space using SPM12 normalise/write, as well as the inverse transformation field that was saved as part of the segmentation procedure mentioned above.
The regions of interest for this study include the left motor cortex (for the motor processing tasks), the bilateral amygdala (for the emotion processing tasks) and the fusiform gyrus (for the emotionProcessing task).
These ROIs are overlaid on the coregistered and resampled T1-weighted image, to allow visual inspection of normalisation quality.

Functional data were preprocessed, starting with estimating realignment parameters for each functional time series using SPM12's realign/estimate,
which performs a 6 degree-of-freedom rigid body transformation that minimizes the sum of squared differences between each volume and the template volume.
Realignment parameters were estimated for the second-echo time series of each run.
Then, slice timing correction was done with SPM12, which corrects for differences in image acquisition time between slices.
Each echo time series of all functional runs were slice time corrected.
3D volume realignment followed, which applied spatial transformation matrices derived from the previously estimated realignment parameters to all echo time series of all functional runs.
Both raw time series and slice time corrected time series were realigned.
Lastly, all echo time series of all functional runs were spatially smoothed using a Gaussian kernel filter with FWHM = 7 mm (i.e. double the voxel size).
Smoothing was performed on raw, slice time corrected and realigned time series data.

Next, several signal time series were calculated or extracted for use as possible GLM regressors in functional task analysis, or for quality control.
From the realignment parameters (3 translation and 3 rotation parameters per volume), a Volterra expansion yielded derivatives, squares and squares of derivatives (Friston et al., 1996).
Framewise displacement (FD, Power et al., 2012) was also calculated from the realignment parameters, and volumes were marked as outliers based on different thresholds of, respectively, 0.2 mm and 0.5 mm.
RETROICOR regressors (Glover et al., 2000) were generated from the cardiac and respiratory signals using the TAPAS PhysIO toolbox, which yielded 6 cardiac regressors, 8 respiratory regressors, 4 interaction regressors,
and additionally a cardiac rate regressor (CR; the cardiac rate time series convolved with the cardiac response function; Chang et al., 2009) and a respiratory volume per time regressor (RVT; respiratory volume per time convolved with the respiratory response function; Birn et al. 2006; Birn et al., 2008).
From the slice time corrected and realigned time series data (of all functional runs), signals were extracted per voxel and spatially averaged within the previously generated tissue masks to yield tissue compartment signals for gray matter, white matter, cerebrospinal fluid (CSF) and the whole brain. 

The last set of preprocessing steps included calculation of image quality metrics and visualizations, using the BIDS-compatible fmrwhy_bids_workflowQC pipeline from the fMRwhy toolbox.
Operations on functional time series data were all done on detrended (linear and quadratic trends) realigned data, except where otherwise specified.
Temporal signal-to-noise ratio (tSNR) maps were calculated for all runs by dividing the voxel-wise time series mean by the voxel-wise standard deviation of the time series.
Tissue compartment averages were then extracted from these tSNR maps. Percentage difference maps (from the time series mean) were calculated per volume for use in carpet plots (or gray plots).