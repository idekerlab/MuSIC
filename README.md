# The Multi-Scale Integrated Cell (MuSIC)

MuSIC is a hierarchical map of human cell architecture created from integrating immunofluorescence images in the Human Protein Atlas with affinity purification experiments from the BioPlex resource. Integration involves configuring each approach to produce a general measure of protein distance, then calibrating the two measures using machine learning.

Comprehensive information for MuSIC available at **https://nrnb.org/music/**.

Additional resources on GitHub:

**[A Step-By-Step Guide to Building a MuSIC Map](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map)**

**[Accompanying Jupyter Notebook to Build MuSIC v1.0](https://github.com/idekerlab/MuSIC/blob/master/Step-by-step%20guide%20to%20build%20MuSIC%20v1.ipynb?)**

**[Accompanying bash script to Build MuSIC v1.0 (example_buid_music_v1.sh)](https://github.com/idekerlab/MuSIC/blob/master/example_buid_music_v1.sh)**

Please cite **[Qin et al., “Mapping cell structure across scales by fusing protein images and interactions”](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1)**.

## Requirements
- [Anaconda](https://www.anaconda.com/products/individual#Downloads) (optional but highly recommended)
- build-essential 
- python-dev 
- libxml2 
- libxml2-dev 
- zlib1g-dev 
- libigraph0-dev 
- libmpc-dev

## Installation
1. Create an Anaconda virtual environment. This is optional but highly recommended. Takes ~10 minutes.
```
conda create -n music python=3.6.2 anaconda
source activate music
```

2. Download MuSIC and install dependencies.

```
git clone https://github.com/idekerlab/MuSIC.git
cd MuSIC
pip install -r ./installation/requirements.txt
```

3. Install [CliXO v1.0](https://github.com/fanzheng10/CliXO-1.0) and [DDOT](https://github.com/michaelkyu/ddot) by running the following command line.

```
./installation/install.sh
```

## Toy example
The toy script will run through MuSIC pipeline using 100 proteins with random embeddins.
```
./toy_example.sh
```

The toy hierarchy is stored in two output files with details available [here](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map#output-file-outprefixlouvainddot).
```
# Output: hierarchical relationship among systems and genes
head ./Examples/toy_output/toy.louvain.ddot

# Output: specific protein assignment for each identified system
head ./Examples/toy_output/toy.louvain.termStats
```
