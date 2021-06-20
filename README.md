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

   * toy_example.sh runs a series of python scripts to infer a hierarchy for the query proteins using random embeddings.
   
        ```
        # Step 1: Generate gold-standard protein-protein proximity values
            python calibrate_pairwise_distance.py --protein_file ./Examples/toy/toy_proteins.txt --outprefix ./Examples/toy_output/toy

        # Step 2: Build random forest to predict protein-protein proximity from data embeddings
            python random_forest_samples.py --outprefix ./Examples/toy_output/toy --protein_file ./Examples/toy/toy_proteins.txt --emd_files ./Examples/toy/toy_IF_image_embedding.csv ./Examples/toy/toy_APMS_embedding.csv --emd_label IF_emd APMS_emd --n_samples 1000

            # run random forest for 5 folds
            for ((fold = 1; fold <= 5; fold++));
            do
                python run_random_forest.py --outprefix ./Examples/toy_output/toy --fold $fold --emd_label IF_emd APMS_emd;
            done

            python random_forest_output.py --outprefix ./Examples/toy_output/toy

        # Step 3: Analyze proximity data to identify protein communities at progressive resolutions
            python community_detection.py --outprefix ./Examples/toy_output/toy --clixo_i ./Examples/toy_output/toy_predicted_proximity.txt --predict_nm_size --keep_all_files
        ```

   * The resulting hierarchy is stored in two output files. Details about the file format can be found [here](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map#output-file-outprefixlouvainddot).

      ```
      # Output: hierarchical relationship among systems and genes
      head ./Examples/toy_output/toy.louvain.ddot

         # column 1: the parent system
         # column 2: the child system or gene
         # column 3: property of child in the second column
            default: child is a system 
            gene: child is a gene

      # Output: specific protein assignment for each identified system
      head ./Examples/toy_output/toy.louvain.termStats

         # column 2 (Number_of_proteins): total number of proteins belonging to the system
         # column 3 (Proteins): comma separated list of proteins belonging to the system
         # column 4 (median_recal_nm): median of predicted distance, in nm, among all pairs of proteins in the system
         # column 5 (Estimated_size_in_nm): predicted size, in nm, of the system
      ```

2. To run the MuSIC pipeline for user-specified input (proteins), follows steps detailed in the following document:
   [A Step-By-Step Guide to Building a MuSIC Map](https://github.com/idekerlab/MuSIC/wiki/A-Step-By-Step-Guide-to-Building-a-MuSIC-Map)**

   Command line example of running the MuSIC pipeline is given in a bash script file: 
   [Accompanying Bash Script to Build MuSIC v1.0 (example_buid_music_v1.sh)](https://github.com/idekerlab/MuSIC/blob/master/example_buid_music_v1.sh)**


3. To run the MuSIC pipeline on jupyter notebook, consider using the following jupyter notebook as a starting point:
   [Accompanying Jupyter Notebook to Build MuSIC v1.0](https://github.com/idekerlab/MuSIC/blob/master/Step-by-step%20guide%20to%20build%20MuSIC%20v1.ipynb?)**



