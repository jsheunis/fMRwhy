.. fMRwhy documentation master file, created by
   sphinx-quickstart on Mon Nov 16 22:51:54 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

fMRwhy: BIDS-compatible fMRI analysis with SPM12
================================================

``fMRwhy`` is a MATLAB- and SPM12-based toolbox with a variety of helper functions and BIDS-compatible workflows to assist in your fMRI quality checking, preprocessing and analysis journey.

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


With ``fMRwhy`` you are provided with the tools and shown *how* to calculate interesting quality metrics, *how* to visualize outcomes, *how* to analyse your data with batch scripts, and *how* to build a BIDS compatible analysis pipeline, all to flexible levels of automation.
``fMRwhy`` does not help with the *why* questions, which are arguably the most important ones that need to be considered right at the start of your fMRI research journey.

``fMRwhy`` currently has the following `features`_:

* Focus on functional MRI
* BIDS-compatibility
* Visual fMRI quality control
* Multi-echo fMRI preprocessing
* Accessible and extensible SPM12 batch processing
* Processing utilities

.. _features: https://fmrwhy.readthedocs.io/features.html

Background
----------

``fMRwhy`` started as a loose collection of MATLAB- and SPM12-based scripts created to reproduce useful quality control measures and interesting visualizations as reported in journal articles,
such as calculating the temporal signal-to-noise ratio of an fMRI time series or creating a carpet plot from fMRI data.
The core idea behind this was: if the results of our research are supported by useful quality control methods, why not implement these in a reproducible and extensible way so that the greater community can benefit from and contribute to it?
This evolved over time into a modular set of SPM12 batch process wrapper functions that simplified the process of creating reproducible fMRI preprocessing and quality reporting scripts.
With the goal of allowing automated quality reporting workflows, BIDS-compatibility was added using ``bids-matlab`` as a dependency. 


License Information
-------------------

``fMRwhy`` is licensed under GNU Lesser General Public License version 2.1.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   features
   installation
   quality
   multiecho
   resources
   contributing
   api






Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
