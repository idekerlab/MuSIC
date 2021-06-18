# The Multi-Scale Integrated Cell (MuSIC)

MuSIC is a hierarchical map of human cell architecture created from integrating immunofluorescence images in the Human Protein Atlas with affinity purification experiments from the BioPlex resource. Integration involves configuring each approach to produce a general measure of protein distance, then calibrating the two measures using machine learning.

This repository contains the source code and data for reproducing the results in [Qin et al. 2020](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1), as well as step-by-step introduction for the MuSIC tool.

The generated MuSIC map and more details are available at [https://nrnb.org/music/](https://nrnb.org/music/).

## Installation

Here are the steps for setting up the enviroment for running MuSIC tool:

- Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), [Anaconda](https://www.anaconda.com/products/individual#Downloads) or other virtual environment manager.

- For Ubuntu (or similar Linux environments) users, run the following command to install the required libraries:
```bash
sudo apt install build-essential python-dev libxml2 libxml2-dev zlib1g-dev libigraph0-dev libmpc-dev
```
- First create a Python virtual environment, for example, in Conda this would be done with:
```
conda create -n music python=3.6.2
```
- Then, we need to install the Data-Driven Ontology Toolkit (DDOT) by following the instructions at: https://github.com/michaelkyu/ddot#install-the-ddot-python-package

- Finally, we can clone the MuSIC git repository and install the required Python packages:
```
git clone https://github.com/idekerlab/MuSIC
cd MuSIC
source activate music # activate the virtual environment
pip install -r ./requirements.txt
```

This code was tested on an Ubuntu 20.04 system.

## Usage

For detailed usage of the MuSIC tool, we provide a [Step-by-step Guide for creating MuSIC maps](./MuSIC-guide.md).

## Reproducing the results
If you want to reproduce the result described in [Qin et al. 2020](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1), please start a Jupyter notebook server and run this notebook: [Reproduce-MuSIC-paper-results.ipynb](./Reproduce-MuSIC-paper-results.ipynb).

## License

[MIT](./LISENCE.txt)

## Citation

If you find this code useful in your research, please consider citing:

**[Qin et al., “Mapping cell structure across scales by fusing protein images and interactions”](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1)**

```
@article {Qin2020.06.21.163709,
	author = {Qin, Yue and Winsnes, Casper F. and Huttlin, Edward L. and Zheng, Fan and Ouyang, Wei and Park, Jisoo and Pitea, Adriana and Kreisberg, Jason F. and Gygi, Steven P. and Harper, J. Wade and Ma, Jianzhu and Lundberg, Emma and Ideker, Trey},
	title = {Mapping cell structure across scales by fusing protein images and interactions},
	year = {2020},
	doi = {10.1101/2020.06.21.163709},
	URL = {https://www.biorxiv.org/content/early/2020/06/22/2020.06.21.163709}
}
```
