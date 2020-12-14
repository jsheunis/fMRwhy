Installation
============

To run ``fMRwhy`` on your local machine, you will first need to have Matlab 2016b or a more recent version installed on your system.

Installation requires a number of steps:

Step 1 - Download ``fMRwhy``
----------------------------

If you are using git, clone the ``fMRwhy`` GitHub repository to your machine:

.. code-block:: bash

  git clone https://github.com/jsheunis/fMRwhy.git

Otherwise, download and extract the zipped code base from the green Code button on the GitHub repository page.


Step 2 - Download dependencies
------------------------------

``fMRwhy`` also has several software dependencies, which you are required to install in order to use all of its features.
These dependencies can be checked by running ``fmrwhy_util_checkDependencies`` from the Matlab command window (after adding ``fMRwhy`` to your Matlab path, see below),
and you will be notified if a dependency is not installed or if it is not on your Matlab path.

For full functionality, ``fMRwhy`` requires installation of the following toolboxes:

* `SPM12`_ 
* `bids-matlab`_
* `dicm2nii`_
* `TAPAS (PhysIO)`_
* `RainCloudPlots`_

Step 3 - Add tools to your Matlab path
--------------------------------------

To add ``fMRwhy`` to your Matlab path, run the following from the Matlab command window:

.. code-block:: bash

  addpath(genpath('path/to/your/fMRwhy/directory'))

where you will have to replace ``path/to/your/fMRwhy/directory`` with the actual path on your system.

This will add ``fMRwhy`` and all subdirectories to the Matlab path.

After installing each dependent toolbox (for which the process will likely be similar to the one for ``fMRwhy``),
please remember to add each toolbox directory to the Matlab path as well.


.. _SPM12: https://github.com/spm/spm12/releases/tag/r7771
.. _bids-matlab: https://github.com/bids-standard/bids-matlab
.. _dicm2nii: https://github.com/jsheunis/dicm2nii/releases/tag/v0.2
.. _TAPAS (PhysIO): https://github.com/translationalneuromodeling/tapas/releases/tag/v4
.. _RainCloudPlots: https://github.com/RainCloudPlots/RainCloudPlots/releases/tag/v1.1

