Core Features
=============


Focus on functional MRI
-----------------------

``fMRwhy`` was created with a core focus on functional MRI and basic structural (i.e. T1-weighted) preprocessing applications.

BIDS-compatibility
------------------

``fMRwhy`` strives to stay compatible with the `Brain Imaging Data Structure`_.
This includes understanding the structure of a BIDS dataset, such as the number and names of subjects, sessions, tasks, runs, and more.
This compatibility allows automated pipelines to be run for any fMRI-based BIDS dataset.
Additionally, ``fMRwhy`` will output derivative data in a BIDS-compatible fashion.

Visual fMRI quality control
---------------------------

This is currently the core functionality. :func:`fmrwhy_workflow_qc` is an automated, BIDS-compatible quality checking and reporting pipeline that generates subject-specific fMRI quality reports. `An example is available here`_.
It requires a settings file to be prepopulated by the user based on the data and the user's preferences for processing steps. It can run on a full BIDS dataset with T1w and BOLD data, and will automatically derive the structure of the data in order to process all tasks, sessions and runs.
Detailed usage information is available in the :ref:`quality_reporting` section, and the function's :ref:`api_ref`.

Multi-echo fMRI preprocessing
-----------------------------

Multi-echo fMRI has known benefits for improving fMRI signal recovery, increasing the signal-to-noise ratio, and separating BOLD and non-BOLD fluctuations.
``fMRwhy`` aims to make multi-echo processing methods accessible to researchers using SPM12 and Matlab.
For this, much inspiration is drawn from the `tedana community`_.

Accessible and extensible SPM12 batch processing
------------------------------------------------

In an attempt to help users move from clicking through a GUI to building reproducible SPM12 analysis scripts,
``fMRwhy`` makes SPM12 batch processes and other functions available as extensible functions. This includes standard SPM12 functionality such as realignment, reslicing, smoothing, and more.

Processing utilities
--------------------

Many helper functions are also available to assist with preprocessing, image calculation and visualization scripts. To find out more, please see the :ref:`api_ref`.


What ``fMRwhy`` is NOT
----------------------

Importantly, while ``fMRwhy`` aims to make these features and future developments available to a wider community of MATLAB users under an open license, it is not (nor is it intended to be) a fully-fledged fMRI analysis package.
It would not suffice nor should it be used as a replacement tool for the likes of SPM12, FSL, fmriprep, AFNI, or their constituent functionality.
It is better positioned as an auxiliary tool for MATLAB and SPM12 users who work with BIDS-formatted fMRI datasets and who want to:

1. Generate automatic, user-friendly quality reports
2. Build reproducible SPM12-based processing pipelines
3. Preprocess multi-echo fMRI data
4. Build/test real-time analysis pipelines
5. Use auxiliary tools for fMRI visualization and image manipulation


.. _Brain Imaging Data Structure: https://bids.neuroimaging.io/
.. _An example is available here: https://jsheunis.github.io/fmrwhy_sample_QCreport.html
.. _tedana community: https://tedana.readthedocs.io/