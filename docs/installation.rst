Installation
============

To run ``fMRwhy`` on your local machine, you will first need to have Matlab 2016b or a more recent version installed on your system.

To install ``fMRwhy``, clone the GitHub repository to your machine:

.. code-block:: bash
  git clone https://github.com/jsheunis/fMRwhy.git

Or download and extract the zipped code base from the green Code button on the GitHub repository page.

To add ``fMRwhy`` to your Matlab path, run the following from the Matlab command window:

.. code-block:: bash
  addpath(genpath('path/to/your/fMRwhy/directory'))

where you will have to replace ``path/to/your/fMRwhy/directory`` with the actual path on your system.
This will add ``fMRwhy`` and all subdirectories to the Matlab path.

Dependencies
------------

``fMRwhy`` also has several software dependencies, which you are required to install in order to use all of its features.
When running the main quality checking workflow ``fmrwhy_bids_worfklowQC``,
these dependencies will be checked and you will be notified if they are not installed or not on your Matlab path.

``fMRwhy`` requires installation of the following toolboxes:

* `SPM12`_ 
* `bids-matlab`_
* `dicm2nii`_
* `TAPAS (PhysIO)`_
* `RainCloudPlots`_

After installing each toolbox (for which the process should be very similar to the one for ``fMRwhy``),
please remember to add each directory to the Matlab path.

.. _SPM12: https://github.com/spm/spm12/releases/tag/r7771
.. _bids-matlab: https://github.com/bids-standard/bids-matlab
.. _dicm2nii: https://github.com/jsheunis/dicm2nii/releases/tag/v0.2
.. _TAPAS (PhysIO): https://github.com/translationalneuromodeling/tapas/releases/tag/v4
.. _RainCloudPlots: https://github.com/RainCloudPlots/RainCloudPlots/releases/tag/v1.1

