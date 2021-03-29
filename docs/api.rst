.. _api_ref:

.. currentmodule:: fmrwhy

API
===

.. Note::  This API is currently under development and documentation is still being updated.

+-------------------------+---------------------------------------------+
| Module                  | Description                                 |
+=========================+=============================================+
| :mod:`fmrwhy.batch`     | Wrapper functions for SPM12 batch processes |
+-------------------------+---------------------------------------------+
| :mod:`fmrwhy.bidsutils` | Utilities for BIDS-compatibility            |
+-------------------------+---------------------------------------------+
| :mod:`fmrwhy.multiecho` | Multi-echo fMRI analysis functions          |
+-------------------------+---------------------------------------------+
| :mod:`fmrwhy.qc`        | Quality control functions                   |
+-------------------------+---------------------------------------------+
| :mod:`fmrwhy.realtime`  | Real-time fMRI analysis functions           |
+-------------------------+---------------------------------------------+
| :mod:`fmrwhy.settings`  | Defaults and workflow-specific settings     |
+-------------------------+---------------------------------------------+
| :mod:`fmrwhy.utils`     | Common utility functions                    |
+-------------------------+---------------------------------------------+
| :mod:`fmrwhy.workflows` | BIDS-compatible analysis pipelines          |
+-------------------------+---------------------------------------------+


|

.. _api_batch_ref:

.. automodule:: fmrwhy.batch

.. currentmodule:: fmrwhy.batch

:mod:`fmrwhy.batch`: Wrapper functions for SPM12 batch processes
----------------------------------------------------------------



|

.. _api_bidsutils_ref:

.. automodule:: fmrwhy.bidsutils

.. currentmodule:: fmrwhy.bidsutils

:mod:`fmrwhy.bidsutils`: Utilities for BIDS-compatibility
---------------------------------------------------------

|

.. autofunction:: fmrwhy.bidsutils.fmrwhy_bids_setupQcDerivDirs
.. autofunction:: fmrwhy.bidsutils.fmrwhy_bids_setupQcSubDirs
.. autofunction:: fmrwhy.bidsutils.fmrwhy_bids_getAnatDerivs
.. autofunction:: fmrwhy.bidsutils.fmrwhy_bids_constructFilename



|

.. _api_multiecho_ref:

.. automodule:: fmrwhy.multiecho

.. currentmodule:: fmrwhy.qmultiechoc

:mod:`fmrwhy.multiecho`: Multi-echo fMRI analysis functions
-----------------------------------------------------------


|

.. _api_qc_ref:

.. automodule:: fmrwhy.qc

.. currentmodule:: fmrwhy.qc

:mod:`fmrwhy.qc`: Quality control functions
-------------------------------------------



|

.. _api_realtime_ref:

.. automodule:: fmrwhy.realtime

.. currentmodule:: fmrwhy.realtime

:mod:`fmrwhy.realtime`: Real-time fMRI analysis functions
---------------------------------------------------------



|

.. _api_settins_ref:

.. automodule:: fmrwhy.settings

.. currentmodule:: fmrwhy.settings

:mod:`fmrwhy.settings`: Defaults and workflow-specific settings
---------------------------------------------------------------

|

.. autofunction:: fmrwhy.settings.fmrwhy_defaults
.. autofunction:: fmrwhy.settings.fmrwhy_settings_validate



|

.. _api_utils_ref:

.. automodule:: fmrwhy.utils

.. currentmodule:: fmrwhy.utils

:mod:`fmrwhy.utils`: Common utility functions
--------------------------------------------------

|

.. autofunction:: fmrwhy.utils.fmrwhy_util_checkDependencies
.. autofunction:: fmrwhy.utils.fmrwhy_util_saveNiftiFrom4D


|

.. _api_workflows_ref:

.. automodule:: fmrwhy.workflows

.. currentmodule:: fmrwhy.workflows

:mod:`fmrwhy.workflows`: BIDS-compatible analysis pipelines
-----------------------------------------------------------

|

.. autofunction:: fmrwhy.workflows.fmrwhy_workflow_qc
.. autofunction:: fmrwhy.workflows.fmrwhy_workflow_qcSubReport
