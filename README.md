# fMRwhy - Facilitating BIDS-compatible fMRI analysis with SPM12 and Matlab
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-1-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

***CURRENTLY UNDER DEVELOPMENT***

---

## What is `fMRwhy`?

This is a Matlab- and SPM12-based toolbox with a variety of helper functions and BIDS-compatible workflows to assist in your fMRI quality checking, preprocessing and analysis journey.

With `fMRwhy` you are shown *how* to calculate interesting quality metrics, *how* to visualize outcomes,
*how* to analyse your data with batch scripts, and *how* to build a BIDS compatible analysis pipeline,
all to flexible levels of automation. 

It cannot help with the *why* questions, which are arguably the most important ones that need to be addressed right at the start of your fMRI research journey.

## Functionality

`fMRwhy` currently supports:

-  `fmrwhy_bids_worfklowQC`: an automated, BIDS-compatible quality checking and reporting pipeline. It requires a settings file to be prepopulated by the user based on the data and the user's preferences for processing steps. It can run on a full BIDS dataset with T1w and BOLD data, and will automatically derive the structure of the data in order to process all tasks, sessions and runs. It then generates a quality control/checking report in HMTL format for each individual functional run. (*further description to be populated*)

- Multi-echo fMRI preprocessing (*further description to be populated*)

Several helper functions are also available (*for which documentation is pending*)

## Installation

To run `fMRwhy` on your local machine, you will need Matlab 2016b (or a more recent version) installed on your system.

To install `fMRwhy`, clone this GitHub repository to your machine:

```
git clone https://github.com/jsheunis/fMRwhy.git
```

Or download and extract the zipped code base from the green `Code` button on this page.

`fMRwhy` also has several software dependencies, which you are required to install. When running `fmrwhy_bids_worfklowQC`, these dependencies will be checked and you will be notified if they are not installed or not on your Matlab path.

To add `fMRwhy` to your Matlab path, run the following from the Matlab command window:

```
addpath(genpath('path/to/your/fMRwhy/directory'))
```

This will add `fMRwhy` and all subdirectories to the Matlab path.

## Dependencies

`fMRwhy` requires installation of the following toolboxes:

- [SPM12](https://github.com/spm/spm12/releases/tag/r7771)
- [bids-matlab](https://github.com/bids-standard/bids-matlab)
- [dicm2nii](https://github.com/jsheunis/dicm2nii/releases/tag/v0.2)
- [TAPAS (PhysIO)](https://github.com/translationalneuromodeling/tapas/releases/tag/v4.0.0)
- [RainCloudPlots](https://github.com/RainCloudPlots/RainCloudPlots/releases/tag/v1.1)

After installing each toolbox (for which the process should be very similar to the one for `fMRwhy`), please remember to add each directory to the Matlab path.

## Usage

In order to run `fmrwhy_bids_worfklowQC` on a BIDS-validated dataset, please follow these steps:

1. Create a `scripts` directory in a location of your choice (and with a name of your choice) in which to save `fMRwhy`-related scripts.
2. Create a copy of the settings file `fMRwhy/workflows/fmrwhy_settings_template.m` and put it in our `scripts` directory. You can rename it to make it more unique and recognisable for your analysis.
3. Update your new settings m-file with information derived from your BIDS dataset and based on your preferences for processing steps. The settings file provides guidance on the required changes, which includes (but is not limited to) aspects like:
    - The BIDS dataset directory location
    - The list of subjects for which you want to run the workflow
    - The template task/session/run/volume for realignment steps
    - Image dimensions
    - Requirements related to regions of interest
    - Inclusion/exclusion of physiological signal processing
    - Inclusion/exclusion of a lis of confounds
4. Run the workflow from the Matlab command window:

```
fmrwhy_bids_worfklowQC('path/to/your/new/settings/file')
```

The fMRwhy `fmrwhy_bids_worfklowQC` will then do the following:
- It creates a `derivatives` directory inside your BIDS dataset directory
- Inside the `derivatives` directory, it creates `fMRwhy`-related output directories: `fmrwhy-preproc` and `fmrwhy-qc`
- For each subject being processed, it copies the raw data from the BIDS directory to the `fmrwhy-preproc` directory.
- It then proceeds with basic preprocessing steps required for the quality control pipeline, the outputs of which are stored in the `fmrwhy-preproc` directory.
- It follows with the quality control pipeline, the outputs of which are stored in the `fmrwhy-qc` directory.
- It ends with the reporting step, the outputs of which are stored in the `fmrwhy-qc` directory. Inside the subject's directory, a report directory named in the format `report_yyyymmddhhmmss` will be created, which contains an html file `sub-XXX_desc-QCreport_yyyymmddhhmmss.html` that can be viewed in your browser.

## Contributing

Contributions to `fMRwhy` are very welcome and encouraged. Ultimately, the goal is for this toolbox to be used and maintained by a community of fMRI researchers.

To provide feedback, report errors, ask questions or suggest improvements, please create a GitHub [issue](https://github.com/jsheunis/fMRwhy/issues)

If you have written code to solve an issue or add a feature/improvement, please fork the repository and submit a
[pull request](https://github.com/jsheunis/fMRwhy/pulls) with the updates.

## Background
This toobox is a culmination of scripts and functions from
[here](https://github.com/jsheunis/matlab-spm-scripts-jsh),
[here](https://github.com/jsheunis/Neu3CA-RT),
[here](https://github.com/jsheunis/fMRI-Quality-Checker), [here](https://github.com/jsheunis/rtme-fMRI-ISMRMb-2019) and other undocumented sources.


## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/jesperr17"><img src="https://avatars1.githubusercontent.com/u/54264865?v=4" width="100px;" alt=""/><br /><sub><b>Jesper Pilmeyer</b></sub></a><br /><a href="https://github.com/jsheunis/fMRwhy/issues?q=author%3Ajesperr17" title="Bug reports">🐛</a> <a href="#ideas-jesperr17" title="Ideas, Planning, & Feedback">🤔</a> <a href="#userTesting-jesperr17" title="User Testing">📓</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!