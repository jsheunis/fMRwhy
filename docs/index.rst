.. fMRwhy documentation master file, created by
   sphinx-quickstart on Mon Nov 16 22:51:54 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

fMRwhy: BIDS-compatible fMRI analysis with SPM12
================================================

.. image:: http://img.shields.io/badge/license-LGPL--2.1-blue.svg
   :target: https://opensource.org/licenses/LGPL-2.1
   :alt: License

.. image:: https://readthedocs.org/projects/fmrwhy/badge/?version=latest
   :target: https://fmrwhy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


About
-----

.. image:: ../assets/fmrwhy_logo_2020.png
  :target: https://fmrwhy.readthedocs.io/

``fMRwhy`` is a Matlab- and SPM12-based toolbox with a variety of helper functions and BIDS-compatible workflows to assist in your fMRI quality checking, preprocessing and analysis journey.

With ``fMRwhy`` you are provided with the tools and shown *how* to calculate interesting quality metrics, *how* to visualize outcomes, *how* to analyse your data with batch scripts, and *how* to build a BIDS compatible analysis pipeline, all to flexible levels of automation.
``fMRwhy`` does not help with the *why* questions, which are arguably the most important ones that need to be considered right at the start of your fMRI research journey.

``fMRwhy`` currently has the following features:

Focus on functional MRI
```````````````````````

``fMRwhy`` was created with a core focus on functional MRI and basic structural (i.e. T1-weighted) preprocessing applications.

BIDS-compatibility
``````````````````

``fMRwhy`` strives to stay compatible with the `Brain Imaging Data Structure`_.
This includes understanding the structure of a BIDS dataset, such as the number and names of subjects, sessions, tasks, runs, and more.
This compatibility allows automated pipelines to be run for any fMRI-based BIDS dataset.
Additionally, ``fMRwhy`` will output derivative data in a BIDS-compatible fashion.

Visual fMRI quality control
````````````````````````````

This is currently the core functionality. ``fmrwhy_bids_workflowQC`` is an automated, BIDS-compatible quality checking and reporting pipeline that generates subject-specific fMRI quality reports. An example is available `here`_.
It requires a settings file to be prepopulated by the user based on the data and the user's preferences for processing steps. It can run on a full BIDS dataset with T1w and BOLD data, and will automatically derive the structure of the data in order to process all tasks, sessions and runs.
Detailed usage information is available here, and the function's API is available here.

Multi-echo fMRI preprocessing
`````````````````````````````

Multi-echo fMRI has known benefits for improving fMRI signal recovery, increasing the signal-to-noise ratio, and separating BOLD and non-BOLD fluctuations.
``fMRwhy`` aims to make multi-echo processing methods accessible to researchers using SPM12 and Matlab.
For this, much inspiration is drawn from the `tedana community`_.

Accessible and extensible SPM12 batch processing
````````````````````````````````````````````````

In an attempt to help users move from clicking through a GUI to building reproducible SPM12 analysis scripts,
``fMRwhy`` makes SPM12 batch processes and other functions available as extensible functions. This includes standard SPM12 functionality such as realignment, reslicing, smoothing, and more.

Processing utilities
`````````````````````

Many helper functions are also available to assist with preprocessing, image calculation and visualization scripts. To find out more, please see the `API`_.


.. _Brain Imaging Data Structure: https://bids.neuroimaging.io/
.. _here: https://jsheunis.github.io/fmrwhy_sample_QCreport.html
.. _tedana community: https://tedana.readthedocs.io/
.. _API: https://fmrwhy.readthedocs.io/en/latest/#api

License Information
-------------------

``fMRwhy`` is licensed under GNU Lesser General Public License version 2.1.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quality
   multiecho
   resources
   contributing
   codeofconduct
   api






Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
