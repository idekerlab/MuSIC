# The Multi-Scale Integrated Cell (MuSIC)

MuSIC is a hierarchical map of human cell architecture created from integrating immunofluorescence images in the Human Protein Atlas with affinity purification experiments from the BioPlex resource. Integration involves configuring each approach to produce a general measure of protein distance, then calibrating the two measures using machine learning.

Web-based exploration of comprehensive information for MuSIC is available at: **https://nrnb.org/music/**.

**[A Step-By-Step Guide to Building a MuSIC Map](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map)**
- **[Accompanying Jupyter Notebook to Build MuSIC v1.0](https://github.com/idekerlab/MuSIC/blob/master/Step-by-step%20guide%20to%20build%20MuSIC%20v1.ipynb?)**
- **[Accompanying Bash Script to Build MuSIC v1.0 (example_buid_music_v1.sh)](https://github.com/idekerlab/MuSIC/blob/master/example_buid_music_v1.sh)**

If you find MuSIC helpful for your research, please cite **[Qin et al., “Mapping cell structure across scales by fusing protein images and interactions”](https://www.nature.com/articles/s41586-021-04115-9)**.


## Set up an environment for MuSIC

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



## MuSIC pipeline execution

1. Test of the pipeline: run the [toy_example.sh](https://github.com/idekerlab/MuSIC/blob/master/toy_example.sh) bash script to execute MuSIC pipeline for a toy example including 100 proteins with random embeddings.
```
./toy_example.sh
```
   * If getting errors like "no modules named tqdm" or "no modules named dill", try reactivate environment with the following command lines:
       * ```
         # reactivate environment
         conda deactivate
         source activate music
         
         # run toy example script again
         ./toy_example.sh
         ```
   * toy_example.sh runs a series of python scripts to infer a hierarchy for the query proteins using random embeddings. The resulting hierarchy is stored in two output files. Details about the file format can be found [here](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map#output-file-2).
        ```
        # Output: hierarchical relationship among systems and genes
        head ./Examples/toy_output/toy.louvain.ddot

             # column 1: the parent system
             # column 2: the child system or gene
             # column 3: property of child in the second column (default indicates column 2 is a system, gene indicates column 2 is a gene)

        # Output: specific protein assignment for each identified system
        head ./Examples/toy_output/toy.louvain.termStats

             # column 1: unique identifier for each system
             # column 2 (Number_of_proteins): total number of proteins belonging to the system
             # column 3 (Proteins): comma separated list of proteins belonging to the system
             # column 4 (median_recal_nm): median of predicted distance, in nm, among all pairs of proteins in the system
             # column 5 (Estimated_size_in_nm): predicted size, in nm, of the system
        ```     

2. To run the MuSIC pipeline for user-specified input (proteins), follow the steps detailed in the following document:
   **[A Step-By-Step Guide to Building a MuSIC Map](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map)**

   Command lines for building MuSIC v1 map as presented in the paper are provided both in a bash script file
   (**[Accompanying Bash Script to Build MuSIC v1.0 (example_buid_music_v1.sh)](https://github.com/idekerlab/MuSIC/blob/master/example_buid_music_v1.sh)**)
   and as a jupyter notebook
   (**[Accompanying Jupyter Notebook to Build MuSIC v1.0](https://github.com/idekerlab/MuSIC/blob/master/Step-by-step%20guide%20to%20build%20MuSIC%20v1.ipynb?)**).
   
   * The CliXO version used in the original study is provided in CliXO_MuSIC.zip. For downloading the latest CliXO version, please follow installation instructions detailed previously. In comparison to CliXO used in the original study, the latest CliXO is faster but could yield suboptimal clusters.




