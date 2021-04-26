# Multi-Scale Integrated Cell (MuSIC)

![Overview](https://github.com/idekerlab/MuSIC/blob/master/Figures/GitHub_overview.png)

#### Comprehensive information for MuSIC available at https://nrnb.org/music/

The eukaryotic cell is a multi-scale structure with modular organization across at least four orders of magnitude. Two central approaches for mapping this structure – protein fluorescent imaging and protein biophysical association – each generate extensive datasets, but of distinct qualities and resolutions that are typically treated separately. Here, we integrate immunofluorescence images in the Human Protein Atlas with affinity purification experiments from the BioPlex resource to create a unified hierarchical map of eukaryotic cell architecture. Integration involves configuring each approach to produce a general measure of protein distance, then calibrating the two measures using machine learning. The evolving map is called the Multi-Scale Integrated Cell (MuSIC).

Please cite **[Qin et al., “Mapping cell structure across scales by fusing protein images and interactions”](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1)**.



## Dependency
Anaconda users please install relevant packages with the following command lines:
```
conda env create -f environment.yml
source activate music
```
To perform pan-resolution community detection as in MuSIC, please install [CliXO v1.0](https://github.com/fanzheng10/CliXO-1.0) and [alignOntology](https://github.com/mhk7/alignOntology). 

We recommend using high performance computing to run the MuSIC pipeline on new datasets (TODO: add details).




## Data embeddings
In MuSIC v1 study, we demonstrated the MuSIC pipeline with embeddings from immunofluorescence images and protein physical association data. However, the application of MuSIC pipeline is not limited to only these two data modalities. The key is to describe each individual protein with respect to the specific measurement platform. Of note, we recommend embedding each protein in the same number of dimensions for different data types.

#### HPA immunofluorescence image embedding
We here provide the 1024-dimension embeddings for the 1,451 images used in MuSIC v1 (/Examples/IF_image_embedding.csv). Please refer to https://github.com/CellProfiling/densenet for image embedding code.

#### BioPlex Protein Embedding
We here provide the 1024-dimension embeddings for the 661 proteins used in MuSIC v1 (/Examples/APMS_embedding.MuSIC.csv), as well as embeddings for all 10,961 BioPlex v2 proteins generated in MuSIC v1 study (https://www.dropbox.com/s/zb1i0vzcsntlcp3/APMS_embedding.BioPlex_v2.csv?dl=0).




## Calibrate protein-protein distance and proximity from Gene Ontology
![Calibration](https://github.com/idekerlab/MuSIC/blob/master/Figures/GitHub_calibration.png)

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





## Random forest prediction of protein distances
![Calibration](https://github.com/idekerlab/MuSIC/blob/master/Figures/GitHub_RandomForest.png)
Using the **Calibration Function**, we label every protein pair with a curated physical distance. With these curated distances as training labels, we will teach random forest regressors to predict the pairwise distance of any protein pair directly from its features embedded from different data modalities.

### Usage
```
python random_forest_samples.py --outprefix /path/to/output/folder/filePrefix
			        --protein_file /path/to/file/contain/proteins/to/analyze
			        --emd_files /path/to/embedding/file1 /path/to/embedding/file2 ...
			        --emd_label emd1 emd2 ...			    
```


#### Required arguments for random_forest_samples.py:
`--outprefix` Full path to the folder where results will be saved in with unique file identifier. Note that this needs to be the same as previous calibration step.

`--protein_file` Path to the file containing list of proteins to analyze. E.g., /Examples/MuSIC_proteins.txt

`--emd_files` Path to each embedding file generated from different data modalities.

`--emd_label` Label for each embedding file. Enter in the order of `--emd_files`. E.g., IF_emd, APMS_emd

#### Optional arguments for random_forest_samples.py:
`--num_set` Number of training sets for each embedding file. Enter in the order of `--emd_files` (default: auto).

`--n_samples` Maximum number of samples to train/test random forest regressor in each fold of k-fold cross validation (default: 1000000).

`--k` Specify k for k-fold cross validation (default: 5).


## Pan-resolution community detection
![Calibration](https://github.com/idekerlab/MuSIC/blob/master/Figures/GitHub_CommunityDetection.png)

### Usage
```
python community_detection.py --outprefix /path/to/output/folder/filePrefix 
                              --path_to_clixo /path/to/CliXO/folder
                              --clixo_i /path/to/clixo/inputFile
                              --path_to_alignOntology /path/to/alignOntology/folder
                              
bash <outprefix>.sh
```
Hierarchy is written in file `<outprefix>.louvain.ddot` with specific protein assignment for each system available in file `<outprefix>.louvain.termStats`.

*Note: maximum clique finding is NP-hard, although utilized a heuristic approach, CliXO can still take a long time to finish.*

#### Required arguments for community_detection.py:
`--outprefix` Full path to the folder where results will be saved in with unique file identifier.

`--path_to_clixo` Full path to CliXO folder.

`--clixo_i` Path to input similarity network for CliXO: a TSV file with three columns. The first two columns should be two strings for the node names (using numbers may cause problem); and the third column should be a value for the edge weight.

`--path_to_alignOntology` Full path to alignOntology folder.

#### Optional arguments for community_detection.py:
`--clixo_a` CliXO -a flag: for the step size of hierarchy construction; usually, a smaller value will create "deeper" hierarchies with more levels from leaves to the root. (default: 0.1)

`--clixo_b` CliXO -b flag: for merging overlapping communities. Two existing communities will be merged if their similarity is above the threshold defined by this value. Usually a higher value will create a hierarchy with more smaller communities", which looks "broader". (default: 0.5)

`--clixo_m` CliXO -m flag: modularity cutoff. Communities lower than this threshold will be removed from the output. Increasing the value will reduce the number of communities. (default: 0)

`--clixo_z` CliXO -z flag: modularity cutoff. Communities lower than this threshold will be removed from the output. Increasing the value will reduce the number of communities. (default: 0)

`--clixo_s` CliXO -s flag: a cutoff of similarity score, if set, the program will terminate when it reaches this point, and stop looking for more terms from scores lower than this threshold. (default: 0)

`--minSystemSize` Minimum number of proteins requiring each system to have. (default: 4)

`--ci_thre` Threshold for containment index. Additional hierarchical parent-child containment relations will be assigned between pairs of systems having a containment index above this threshold. (default: 0.75)

`--ji_thre` Threshold for Jaccard index. System having Jaccard index above this threshold with its parent system will be removed. (default: 0.9)

`--niter` Number of iterations Louvain clustering will run to select partition with the best modularity. (default: 1000)

`--min_diff` Minimum difference in number of proteins required for every parent-child pair. (default: 1)

`--keep_all_files` When this flag is provided, all intermediate output files will be kept.

## Calibrate to physical distance

This script calibrates the hierarchy to physical distance.

### Usage 
```
python quantile_regression.py --q 0.2 --balance --outdir /path/to/output/folder --rf_pred_path /path/to/file/avg_allyPred.npy --sdir /path/to/output/folder/train_test_data --hierarchy_path /path/to/hierarchy.louvain.termStats --predicted_sim_path /path/to/predicted_resnik_sim.ddot
 ```

### Required arguments for quantile_regression.py

`--q` quantile for quantile regression model

`--outdir` Full path to the folder where results will be saved in with unique file identifier.

`--rf_pred_path` Full path to file with prefix avg_allyPred_modelfname_prefix.npy saved in pred_RF.py

`--sdir` Full path to the folder training and testing data was saved in format_random_forest.py 

`--hierarchy_path` Full path to file with hierarchy saved in community_detection.py with file extension louvain.termStats

`--predicted_sim_path` Full path to random forest predicted similarity matrix predicted_resnik_sim.ddot

### Optional arguments for quantile_regression.py

`--balance` 

