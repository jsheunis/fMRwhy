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

.. image:: https://github.com/jsheunis/fMRwhy/blob/master/assets/fmrwhy_logo_2020.png
  :target: https://fmrwhy.readthedocs.io/

``fMRwhy`` is a Matlab- and SPM12-based toolbox with a variety of helper functions and BIDS-compatible workflows to assist in your fMRI quality checking, preprocessing and analysis journey.

With ``fMRwhy`` you are provided with the tools and shown *how* to calculate interesting quality metrics, *how* to visualize outcomes, *how* to analyse your data with batch scripts, and *how* to build a BIDS compatible analysis pipeline, all to flexible levels of automation.

``fMRwhy`` does not help with the *why* questions, which are arguably the most important ones that need to be considered right at the start of your fMRI research journey.

``fMRwhy`` currently supports:

* fmrwhy_bids_worfklowQC:
  an automated, BIDS-compatible quality checking and reporting pipeline. It requires a settings file to be prepopulated by the user based on the data and the user's preferences for processing steps. It can run on a full BIDS dataset with T1w and BOLD data, and will automatically derive the structure of the data in order to process all tasks, sessions and runs. It then generates a quality control/checking report in HMTL format for each individual functional run. *(further description to be populated)*
* Multi-echo fMRI preprocessing:
  *(further description to be populated)*

Several helper functions are also available *(documentation pending)*


License Information
-------------------

``fMRwhy`` is licensed under GNU Lesser General Public License version 2.1.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   api






Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
