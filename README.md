# The Multi-Scale Integrated Cell (MuSIC)

MuSIC is a hierarchical map of human cell architecture created from integrating immunofluorescence images in the Human Protein Atlas with affinity purification experiments from the BioPlex resource. Integration involves configuring each approach to produce a general measure of protein distance, then calibrating the two measures using machine learning.

Web-based exploration of comprehensive information for MuSIC is available at **https://nrnb.org/music/**.

If your research utilizes the MuSIC hierarchy or a customized hierarchy constructed using the MuSIC pipeline, please cite **[Qin et al., “Mapping cell structure across scales by fusing protein images and interactions”](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1)**.


# Set up an environment for MuSIC

0. Requirements
- [Anaconda](https://www.anaconda.com/products/individual#Downloads) (optional but highly recommended)
- APT packages including build-essential, python-dev, libxml2, libxml2-dev, zlib1g-dev, libigraph0-dev, libmpc-dev


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

3. Install hierarchy building softwares, [CliXO v1.0](https://github.com/fanzheng10/CliXO-1.0) and [DDOT](https://github.com/michaelkyu/ddot), by running the following command line.

```
./installation/install.sh
```


# MuSIC pipeline execution

1. Test of the pipeline: run the following bash script to execute MuSIC pipeline for a toy example including 100 proteins with random embeddings.
```
./toy_example.sh
```
  i. The bash script runs a series of python scripts to build a hierarchy for the queried proteins.
    ```
    a) calibrate_pairwise_distance.py: generate gold-standard protein-protein proximity values
    b) random_forest_samples.py, run_random_forest.py, random_forest_output.py: build random forest to predict protein-protein proximity from data embeddings
    c) community_detection.py: analyze proximity data to identify protein communities at progressive resolutions
    ```
    
  ii. The resulting hierarchy is stored in two output files. Details about the file format can be found [here](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map#output-file-outprefixlouvainddot).
    ```
    # Output: hierarchical relationship among systems and genes
    head ./Examples/toy_output/toy.louvain.ddot

    # Output: specific protein assignment for each identified system
    head ./Examples/toy_output/toy.louvain.termStats
    ```

2. To run the MuSIC pipeline for user-specified input (proteins), follows steps detailed in the following document:
[A Step-By-Step Guide to Building a MuSIC Map](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map)**

    Command line example of running the MuSIC pipeline is given in a bash script file:
[Accompanying Bash Script to Build MuSIC v1.0 (example_buid_music_v1.sh)](https://github.com/idekerlab/MuSIC/blob/master/example_buid_music_v1.sh)**


3. To run the MuSIC pipeline on jupyter notebook, consider using the following jupyter notebook as a starting point:
[Accompanying Jupyter Notebook to Build MuSIC v1.0](https://github.com/idekerlab/MuSIC/blob/master/Step-by-step%20guide%20to%20build%20MuSIC%20v1.ipynb?)**



