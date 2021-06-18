# Step-by-step Guide for creating MuSIC maps

This page describes a step-by-step guide to create the MuSIC map, including how to install the necessary packages and other dependencies; the specification of input data; and detailed documentation for each step and its parameters. An [accompanying Jupyter notebook](https://github.com/idekerlab/MuSIC/blob/master/Step-by-step%20guide%20to%20build%20MuSIC%20v1.ipynb) instantiates these steps with the exact series of code executions that were used to generate MuSIC v1 presented in the paper. Notably, these steps can be used with other inputs and parameter settings to build other multi-scale cell models.

Before following this guide, please make sure you setup the Python environment according to the [instructions for installation](https://github.com/idekerlab/MuSIC/tree/master#installation).

## Table of contents
[Input Data Preparation](#input-data-preparation)

[Step 1: Generate gold-standard protein-protein proximity values](#step-1-generate-gold-standard-protein-protein-proximity-values)

[Step 2: Build random forest to predict protein-protein proximity from data embeddings](#step-2-build-random-forest-to-predict-protein-protein-proximity-from-data-embeddings)

[Step 3: Analyze proximity data to identify protein communities at progressive resolutions](#step-3-analyze-proximity-data-to-identify-protein-communities-at-progressive-resolutions)

[Additional Resources](#additional-resources)

## Input Data Preparation
In MuSIC v1 study, we demonstrated the MuSIC pipeline with embeddings from immunofluorescence images and protein physical association data. However, the application of MuSIC pipeline is not limited to only these two data modalities. The key is to describe each individual protein with respect to the specific measurement platform. Of note, we recommend embedding each protein in the same number of dimensions for different data types.

To run MuSIC, you need to first generate the two types of protein embeddings from immunofluorescence(IF) images and protein-protein interactions(PPI). 

As described in the following sections, we use IF images from [HPA](https://www.proteinatlas.org/humanproteome/cell) and PPI data from [BioPlex](https://bioplex.hms.harvard.edu/).

#### HPA immunofluorescence image embedding
We here provide the 1024-dimension embeddings for the 1,451 images used in MuSIC v1 (`/Examples/IF_image_embedding.csv`). Please refer to https://github.com/CellProfiling/densenet for image embedding code and contact [Casper Winsnes](mailto:casper.winsnes@scilifelab.se?subject=[MuSIC%20GitHub]%20Image%20embeddings) for related questions.

#### BioPlex protein embedding
We here provide the 1024-dimension embeddings for the 661 proteins used in MuSIC v1 (`/Examples/APMS_embedding.MuSIC.csv`), as well as embeddings for all 10,961 BioPlex v2 proteins generated in MuSIC v1 study (https://www.dropbox.com/s/zb1i0vzcsntlcp3/APMS_embedding.BioPlex_v2.csv?dl=0).



#### Customized embedding file format
Besides the provided embeddings, you can also use your own custom embedding files generated from other data sources.

For customized embedding files, please format file to match the style in the example embedding files:
- First column: unique ID for each embedding
- Second column: gene name
- Following columns: each column is one entry in the embedding vector
- Comma separate all columns

***Note: The MuSIC pipeline can handle any number of data types (e.g. IF, APMS, etc.) and any dimension of protein embeddings (i.e. length of feature vector), but all proteins need to have same number of dimension within each individual data type, and we recommend keeping consistent number of dimensions among different types of data.***


## Step 1: Generate gold-standard protein-protein proximity values
<img src="https://github.com/idekerlab/MuSIC/blob/master/Figures/GitHub_calibration.png" width="1000">

As a means of calibrating distance in the embeddings to physical distance in cells, we sampled the literature to assemble a reference set of ten subcellular components with known physical sizes, from protein complexes of <20 nm to organelles >1 µm in diameter. The size of each of these ten components strongly correlated with its number of protein species documented in the Gene Ontology (GO), suggesting a general approximate conversion from the number of protein species to diameter, in nanometers, of a cellular component (**Calibration Function**).

### Usage
```
python calibrate_pairwise_distance.py --protein_file /path/to/file/contain/proteins/to/analyze 
			   	      --outprefix /path/to/output/folder/filePrefix 
```
#### Required arguments for calibrate_pairwise_distance.py:
`--protein_file` Path to the file containing list of proteins to analyze. E.g., /Examples/MuSIC_proteins.txt

`--outprefix` Full path to the folder where results will be saved in with unique file identifier.

#### Optional arguments for community_detection.py:
`--C_file` Path to precomputed matrix containing the size (number of proteins) of smallest Gene Ontology component shared by the gene pair. Computing protein pairwise C matrix can be time consuming when analyzing a large number (a few thousands) of proteins, so we also provide pre-computed protein pairwise C matrix for all 18,395 proteins in Gene Ontology (https://www.dropbox.com/s/w9lxpnw39g64zs8/all_protein_min_GO_size.txt?dl=0).

### Output file: $outprefix.calibrated_distance.csv

This file stores the calibrated distance and proximity for every pair of input proteins that have annotation in Gene Ontology (GO). Below are specific annotation for each column:
- **geneA:** name of gene A
- **geneB:** name of gene B
- **C:** the number of proteins in the smallest GO cellular component to which both are annotated
- **D:** protein-protein distance calibrated from GO
- **log10D:** log10(D)
- **P:** protein-protein proximity, -log10(D)



## Step 2: Build random forest to predict protein-protein proximity from data embeddings
![Calibration](https://github.com/idekerlab/MuSIC/blob/master/Figures/GitHub_RandomForest.png)
Using the **Calibration Function**, we label every protein pair with a curated physical distance. With these curated distances as training labels, we will teach random forest regressors to predict the pairwise distance of any protein pair directly from its features embedded from different data modalities. Depending on number of features and samples, training a random forest regressor could take a long time and need a large amount of computational resource. For example, each random forest regressor in the original MuSIC study was trained with ~1M samples consisted of 2060 input features, requiring ~1 day and >100 Gb memory with 24 threads. 

#### Note
1) Original MuSIC study used only immunofluorescence images and protein physical association data. However, the MuSIC pipeline design is generalizable to any number of data modalities and here we allow user input any number of data modalities to integrate.
2) We recommend keeping number of input features to random forest model the same for each data modality.
3) The amount of training time and computational resource can be decreased by using less number of samples and embedding data modalities into smaller number of dimensions.

### Usage
```
python random_forest_samples.py --outprefix /path/to/output/folder/filePrefix
			        --protein_file /path/to/file/contain/proteins/to/analyze
			        --emd_files /path/to/embedding/file1 /path/to/embedding/file2 ...
			        --emd_label emd1 emd2 ...			    
```
To facilitate concurrent training of multiple random forest regressors with high performance computing, we here provide stand alone script for training and predicting random forest model.

```
python run_random_forest.py --outprefix /path/to/output/folder/filePrefix
			    --fold 1
			    --emd_label emd1 emd2 ...
```
After all random forest models finished training, the final protein-protein proximity is determined by the average of all availabel predictions.
```
python random_forest_output.py --outprefix /path/to/output/folder/filePrefix
```
The <outprefix>_predicted_proximity.txt file is used as input for the pan-resolution community detection.

#### Required arguments for random_forest_samples.py:
`--outprefix` Full path to the folder where results will be saved in with unique file identifier. Note that this needs to be the same as previous calibration step.

`--protein_file` Path to the file containing list of proteins to analyze. E.g., /Examples/MuSIC_proteins.txt

`--emd_files` Path to each embedding file generated from different data modalities.

`--emd_label` Label for each embedding file. Enter in the order of `--emd_files`. E.g., IF_emd, APMS_emd

#### Optional arguments for random_forest_samples.py:
`--num_set` Number of training sets for each embedding file. Enter in the order of `--emd_files` (default: auto).

`--n_samples` Maximum number of samples to train/test random forest regressor in each fold of k-fold cross validation (default: 1000000).

`--k` Specify k for k-fold cross validation (default: 5).

#### Required arguments for run_random_forest.py:
`--outprefix` Full path to the folder where results will be saved in with unique file identifier. Note that this needs to be the same as previous calibration step.

`--fold` Specify which fold of k-fold cross validation to train.

`--emd_label` Label for each embedding file. Enter in the order of `--emd_files`. E.g., IF_emd, APMS_emd

#### Optional arguments for run_random_forest.py:
`--train_set` For each embedding data, specify which training set to use. In order with `emd_label`. Default to 1 for each emd_label.

`--n_estimators` The number of trees in the forest (default: 1000).

`--max_depth` The maximum depth of the tree (default: 30).

`--n_jobs` The number of jobs to run in parallel for training random forest regressor (default: 8).

### Output file: $outprefix_predicted_proximity.txt

This file stores the predicted protein-protein proximity for all pairs of the given protein.
- First column: name of gene A
- Second column: name of gene B
- Third column: predicted proximity between gene A and gene B



## Step 3: Analyze proximity data to identify protein communities at progressive resolutions

<img src="https://github.com/idekerlab/MuSIC/blob/master/Figures/GitHub_CommunityDetection.png" width="1000">

Protein communities were identified at multiple resolutions, starting with those that form at the smallest protein-protein distances then progressively relaxing the distance threshold (multi-scale community detection). Communities at smaller distances were contained, in full or in part, inside larger communities as the threshold was relaxed, yielding a structural hierarchy.

Here, we require an approach to network community detection that has two key properties. First is the ability to identify protein communities that form at multiple stringencies of analysis, resulting in a spectrum of community sizes (multi-scale). Second is the ability to identify communities overlapping at multiple extents (i.e., which share common members), including communities that nest partially or wholly within multiple others (pleiotropy and multi-localization). These properties reflect the biological reality that cells are multi-scale structures with many subcomponents, each of which may be involved in multiple higher-order or pleiotropic processes. Note that maximum clique finding is NP-hard. Although utilizing a heuristic approach, CliXO can still take a long time to finish. For MuSIC, pan-resolution community detecion took around 7 hours.

### Usage
Before running the script for community detection, you need to download the CliXO v1.0 dataset from [here](https://github.com/fanzheng10/CliXO-1.0/releases) and place the unzipped folder to somewhere. Below you will need to specify the path to the CliXO folder.

Now you can run the community detection script:
```
python community_detection.py --outprefix /path/to/output/folder/filePrefix 
                              --path_to_clixo /path/to/CliXO/folder
                              --clixo_i /path/to/clixo/inputFile
                              --path_to_alignOntology /path/to/alignOntology/folder
			      --predict_nm_size                         
bash <outprefix>.sh
```
Hierarchy is written in file `$outprefix.louvain.ddot` with specific protein assignment for each system available in file `$outprefix.louvain.termStats`.

#### Required arguments for community_detection.py:
`--outprefix` Full path to the folder where results will be saved in with unique file identifier.

`--path_to_clixo` Full path to CliXO folder.

`--clixo_i` Path to input similarity network for CliXO: a TSV file with three columns. The first two columns should be two strings for the node names (using numbers may cause problem); and the third column should be a value for the edge weight.

`--path_to_alignOntology` Full path to alignOntology folder.

#### Optional arguments for community_detection.py:
`--clixo_a` CliXO -a flag: for the step size of hierarchy construction; usually, a smaller value will create "deeper" hierarchies with more levels from leaves to the root (default: 0.1).

`--clixo_b` CliXO -b flag: for merging overlapping communities. Two existing communities will be merged if their similarity is above the threshold defined by this value. Usually a higher value will create a hierarchy with more smaller communities", which looks "broader" (default: 0.5). 

`--clixo_m` CliXO -m flag: modularity cutoff. Communities lower than this threshold will be removed from the output. Increasing the value will reduce the number of communities (default: 0).

`--clixo_z` CliXO -z flag: modularity cutoff. Communities lower than this threshold will be removed from the output. Increasing the value will reduce the number of communities (default: 0).

`--clixo_s` CliXO -s flag: a cutoff of similarity score, if set, the program will terminate when it reaches this point, and stop looking for more terms from scores lower than this threshold (default: 0).

`--minSystemSize` Minimum number of proteins requiring each system to have (default: 4).

`--ci_thre` Threshold for containment index. Additional hierarchical parent-child containment relations will be assigned between pairs of systems having a containment index above this threshold (default: 0.75). 

`--ji_thre` Threshold for Jaccard index. System having Jaccard index above this threshold with its parent system will be removed (default: 0.9). 

`--niter` Number of iterations Louvain clustering will run to select partition with the best modularity (default: 1000).

`--min_diff` Minimum difference in number of proteins required for every parent-child pair (default: 1).

`--keep_all_files` When this flag is provided, all intermediate output files will be kept.

`--n_samples` Maximum number of samples to use for fitting linear model (default: 1000000).

`--q_step` Step for scanning best quantile to use (default: 0.1).

`--predict_nm_size` When this flag is provided, all systems will have an estimated size in nm. Note that this calcualtion requries <outprefix>_avgPred_ytrue.csv generated from the random_forest_output.py script in the previous step.

### Output file: $outprefix.louvain.ddot

This file stores the hierarchical relationship among systems and genes:
- First column: the parent system 
- Second column: the child system or gene
- Third column: property of child in the second column
    - default: child is a system
    - gene: child is a gene

### Output file: $outprefix.louvain.termStats

This file stores the specific protein assignment for each identified system. First column is the unique identifier for each system. Other columns are annotated below:
- **Number_of_proteins:** total number of proteins belonging to the system
- **Proteins:** comma separated list of proteins belonging to the system
- **median_recal_nm:** median of predicted distance, in nm, among all pairs of proteins in the system
- **Estimated_size_in_nm:** predicted size, in nm, of the system
	
## Additional Resources
MuSIC v1 was visualized and explored with Cytoscape (session file available [here](https://www.ndexbio.org/?#/network/7fc70ab6-9fb1-11ea-aaef-0ac135e8bacf?accesskey=68afa0480a4859906b5d221619ee95679da96059680557f65c3dd9f1842e4930)). Ideker lab has also developed a hierarchical visualization webapp, [HiView](http://hiview.ucsd.edu/), for visualizing hierarchical models and the underlying interaction data that support them.
