_[Note: this repository is not actively maintained, although notifications are monitored]_

# fMRwhy - Facilitating BIDS-compatible fMRI analysis with SPM12 and Matlab

[![All Contributors](https://img.shields.io/badge/all_contributors-3-orange.svg?style=flat-square)](#contributors-)
[![Documentation Status](https://readthedocs.org/projects/fmrwhy/badge/?version=latest)](https://fmrwhy.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/license-LGPL--2.1-blue)](https://opensource.org/licenses/LGPL-2.1)

---

## What is `fMRwhy`?

<img src="assets/fmrwhy_logo_2020.png" alt=""/>

<br>

This is a Matlab- and SPM12-based toolbox with a variety of helper functions and BIDS-compatible workflows to assist in your fMRI quality checking, preprocessing and analysis journey.

With `fMRwhy` you are shown *how* to calculate interesting quality metrics, *how* to visualize outcomes,
*how* to analyse your data with batch scripts, and *how* to build a BIDS compatible analysis pipeline,
all to flexible levels of automation. 

It cannot help with the *why* questions, which are arguably the most important ones that need to be addressed right at the start of your fMRI research journey.

## Functionality

`fMRwhy` currently supports:

-  `fmrwhy_workflow_qc`: an automated, BIDS-compatible quality checking and reporting pipeline. It requires a settings file to be prepopulated by the user based on the data and the user's preferences for processing steps. It can run on a full BIDS dataset with T1w and BOLD data, and will automatically derive the structure of the data in order to process all tasks, sessions and runs. It then generates a quality control/checking report in HMTL format for each individual functional run. (*further description to be populated*)

- Multi-echo fMRI preprocessing (*further description to be populated*)

Several helper functions are also available (*for which documentation is pending*)

## Installation

To run `fMRwhy` on your local machine, you will need MATLAB 2016b or a more recent version installed on your system. Installing `fMRwhy` requires downloading the toolbox along with a list of required dependencies:

- [SPM12](https://github.com/spm/spm12/releases/tag/r7771)
- [bids-matlab](https://github.com/bids-standard/bids-matlab)
- [dicm2nii](https://github.com/jsheunis/dicm2nii/releases/tag/v0.2)
- [TAPAS (PhysIO)](https://github.com/translationalneuromodeling/tapas/releases/tag/v4.0.0)
- [RainCloudPlots](https://github.com/RainCloudPlots/RainCloudPlots/releases/tag/v1.1)

### To install `fMRwhy` using git

First clone this GitHub repository to your machine once located in your preferred directory:

```bash
git clone https://github.com/jsheunis/fMRwhy.git
```

and then install all dependencies (for reference, see the list below):

```bash
cd fMRwhy
git submodule update --init
```

This could take several minutes, because the dependencies are all downloaded and saved to your local `fMRwhy/dependencies` directory.

Lastly, open MATLAB, add theccode directory (`fMRwhy/fmrwhy`) to the MATLAB path, and run `fmrwhy_init.m`:

```matlab
addpath('your/path/to/fMRwhy/fmrwhy')
fmrwhy_init()
```

This adds `fMRwhy` and all required dependencies to the MATLAB path.

### To install `fMRwhy` manually

Download and extract the zipped code base from the green `Code` button on this page. Once extracted, rename the folder to `fMRwhy` (exact case match and remove any extra text, e.g. `-master`).

Next, following the same process, download+extract (+rename where applicable) all dependencies listed above.

Then, importantly, move all dependencies into the `fMRwhy/dependencies` directory. The toolbox will search this location to see if all required dependencies are available.

Lastly, open MATLAB, add theccode directory (`fMRwhy/fmrwhy`) to the MATLAB path, and run `fmrwhy_init.m`:

```matlab
addpath('your/path/to/fMRwhy/fmrwhy')
fmrwhy_init()
```
This adds `fMRwhy` and all required dependencies to the MATLAB path.

## Usage

In order to run `fmrwhy_workflow_qc` on a BIDS-validated dataset, please follow these steps:

1. Create a `scripts` directory in a location of your choice (and with a name of your choice) in which to save `fMRwhy`-related scripts.
2. Create a copy of the settings file `fMRwhy/settings/fmrwhy_settings_template.m` and put it in our `scripts` directory. You can rename it to make it more unique and recognisable for your analysis.
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
fmrwhy_workflow_qc('path/to/your/new/settings/file')
```

The `fmrwhy_workflow_qc` pipeline will then do the following:
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


### Code style guide and quality

We use the [MISS_HIT linter](https://github.com/florianschanda/miss_hit/) to
automatically detect programming and stylistic errors, to enforce styling and other rules, and to check for code quality.

The linter is a Python package, which means that we use Python to check on MATLAB :smile:. This also implies that you need a local Python environment set up in order to make use of the linter. Once your Python environment is ready, MISS_HIT can be installed with:

```
pip3 install -r requirements.txt
```

where the `requirements.txt` file is located in the `fMRwhy` root directory.

The rules followed by the MISS_HIT are already specified in the
[MISS_HIT configuration file](./miss_hit.cfg).

To check the code style of the whole repository, you can can simply type the following in the root directory:

```
mh_style .
```

Some styling issues can be automatically fixed by using the `--fix` flag. You
might be required to rerun this command several times if multiple issues are flagged.

```
mh_style . --fix
```

Code quality can be checked with:

```
mh_metric .
```

To see only the issues that "break" the code quality rules (also specified in the
configuration file) type:

```
mh_metric . --ci
```

Note that you do not have to have this linter installed locally in order for your contributions to be checked. The code style and quality are also checked during the
[continuous integration](.github/workflows/miss_hit.yml) process, which runs automatically on all pushes to the `master` branch and on pull requests to all branches.

For more information about MISS_HIT, please refer to its
[documentation](https://florianschanda.github.io/miss_hit/).

### Pre-commit

Analogous to the continuous integration process running on GitHub, you can also automatically run the MISS_HIT linter on your local commits. We use [pre-commit hook](https://pre-commit.com/) for this purpose, because it can reformat files as you commit them locally.

As with MISS_HIT, pre-commit is a Python package and needs to be installed in your local Python development environment. As with MISS_HIT, you can install pre-commit by using our `requirements.txt` file 
```bash
pip install -r requirements.txt
```

Then, install the hook:
```bash
pre-commit install
```

You're done. Now, `mh_style --fix` will run every time you commit a local change.

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
    <td align="center"><a href="https://github.com/jesperr17"><img src="https://avatars1.githubusercontent.com/u/54264865?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jesper Pilmeyer</b></sub></a><br /><a href="https://github.com/jsheunis/fMRwhy/issues?q=author%3Ajesperr17" title="Bug reports">🐛</a> <a href="#ideas-jesperr17" title="Ideas, Planning, & Feedback">🤔</a> <a href="#userTesting-jesperr17" title="User Testing">📓</a></td>
    <td align="center"><a href="https://remi-gau.github.io/"><img src="https://avatars3.githubusercontent.com/u/6961185?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Remi Gau</b></sub></a><br /><a href="#ideas-Remi-Gau" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/jsheunis/fMRwhy/commits?author=Remi-Gau" title="Code">💻</a> <a href="https://github.com/jsheunis/fMRwhy/commits?author=Remi-Gau" title="Tests">⚠️</a> <a href="#question-Remi-Gau" title="Answering Questions">💬</a> <a href="#userTesting-Remi-Gau" title="User Testing">📓</a></td>
    <td align="center"><a href="https://jsheunis.github.io/"><img src="https://avatars0.githubusercontent.com/u/10141237?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Stephan Heunis</b></sub></a><br /><a href="https://github.com/jsheunis/fMRwhy/issues?q=author%3Ajsheunis" title="Bug reports">🐛</a> <a href="https://github.com/jsheunis/fMRwhy/commits?author=jsheunis" title="Code">💻</a> <a href="https://github.com/jsheunis/fMRwhy/commits?author=jsheunis" title="Documentation">📖</a> <a href="#ideas-jsheunis" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-jsheunis" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-jsheunis" title="Maintenance">🚧</a> <a href="#projectManagement-jsheunis" title="Project Management">📆</a> <a href="#question-jsheunis" title="Answering Questions">💬</a> <a href="#userTesting-jsheunis" title="User Testing">📓</a></td>
    <td align="center"><a href="https://luiseudave.netlify.app/"><img src="https://avatars.githubusercontent.com/u/38909206?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Luis Eudave</b></sub></a><br /><a href="https://github.com/jsheunis/fMRwhy/issues?q=author%3Anegatoscope" title="Bug reports">🐛</a> <a href="https://github.com/jsheunis/fMRwhy/commits?author=negatoscope" title="Code">💻</a> <a href="https://github.com/jsheunis/fMRwhy/commits?author=negatoscope" title="Documentation">📖</a> <a href="#ideas-negatoscope" title="Ideas, Planning, & Feedback">🤔</a> <a href="#question-negatoscope" title="Answering Questions">💬</a> <a href="#userTesting-negatoscope" title="User Testing">📓</a></td>
    <td align="center"><a href="https://huijbers.github.io/"><img src="https://avatars.githubusercontent.com/u/17592262?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Willem Huijbers</b></sub></a><br /><a href="https://github.com/jsheunis/fMRwhy/issues?q=author%3Ahuijbers" title="Bug reports">🐛</a> <a href="#ideas-huijbers" title="Ideas, Planning, & Feedback">🤔</a> <a href="#question-huijbers" title="Answering Questions">💬</a> <a href="#userTesting-huijbers" title="User Testing">📓</a></td>
    <td align="center"><a href="https://github.com/lydiatgit"><img src="https://avatars.githubusercontent.com/u/453067?v=4?s=100" width="100px;" alt=""/><br /><sub><b>lydiatgit</b></sub></a><br /><a href="https://github.com/jsheunis/fMRwhy/issues?q=author%3Alydiatgit" title="Bug reports">🐛</a> <a href="#ideas-lydiatgit" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/jsheunis/fMRwhy/commits?author=lydiatgit" title="Code">💻</a> <a href="#question-lydiatgit" title="Answering Questions">💬</a> <a href="#userTesting-lydiatgit" title="User Testing">📓</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
