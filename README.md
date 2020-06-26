# Multi-Scale Integrated Cell (MuSIC)
This repo contains code for multi-scale community detection.

## Requirements
1. [CliXO v1.0](https://github.com/fanzheng10/CliXO-1.0)
2. [louvain-igraph v0.6.1](https://github.com/vtraag/louvain-igraph)
3. [DDOT](https://github.com/michaelkyu/ddot)
4. [alignOntology](https://github.com/mhk7/alignOntology)
5. pandas
6. numpy
7. argparse
8. networkx ≥v2.3
9. igraph

## Usage
```
python community_detection.py --outprefix /path/to/output/folder/filePrefix 
                              --path_to_clixo /path/to/CliXO/folder
                              --clixo_i /path/to/clixo/inputFile
                              --path_to_alignOntology /path/to/alignOntology/folder
                              
bash <outprefix>.sh
```
Hierarchy is written in file `<outprefix>.louvain.ddot` with specific protein assignment for each system available in file `<outprefix>.louvain.termStats`.

*Note: maximum clique finding is NP-hard, although utilized a heuristic approach, CliXO can still take a long time to finish.*

#### Required parameters for community_detection.py:
`--outprefix` Full path to the folder where results will be saved in with unique file identifier. E.g. /path/to/output/folder/filePrefix 

`--path_to_clixo` Full path to CliXO folder. E.g. /path/to/CliXO/folder

`--clixo_i` Path to input similarity network for CliXO: a TSV file with three columns. The first two columns should be two strings for the node names (using numbers may cause problem); and the third column should be a value for the edge weight.

`--path_to_alignOntology` Full path to alignOntology folder.

#### Optional parameters for community_detection.py:
`--clixo_a` CliXO -a flag: for the step size of hierarchy construction; usually, a smaller value will create "deeper" hierarchies with more levels from leaves to the root. (default: 0.1)

`--clixo_b` CliXO -b flag: for merging overlapping communities. Two existing communities will be merged if their similarity is above the threshold defined by this value. Usually a higher value will create a hierarchy with more smaller communities", which looks "broader". (default: 0.5)

`--clixo_m` CliXO -m flag: modularity cutoff. Communities lower than this threshold will be removed from the output. Increasing the value will reduce the number of communities. (default: 0)

`--clixo_z` CliXO -z flag: modularity cutoff. Communities lower than this threshold will be removed from the output. Increasing the value will reduce the number of communities. (default: 0)

`--clixo_s` CliXO -s flag: a cutoff of similarity score, if set, the program will terminate when it reaches this point, and stop looking for more terms from scores lower than this threshold. (default: 0)

`--minSystemSize` Minimum number of proteins requiring each system to have. (default: 4)

`--ci_thre` Threshold for Containment Index. (default: 0.75)

`--ji_thre` Threshold for Jaccard Index. (default: 0.9)

`--niter` Number of iterations Louvain clustering will run to select partition with the best modularity. (default: 1000)

`--min_diff` Minimum difference in number of proteins for every parent-child pair. (default: 1)

`--keep_all_files` When this flag is provided, all intermediate output files will be kept.


## Webpage
https://nrnb.org/music/

## Citing MuSIC
Please cite **[Qin et al., “Mapping cell structure across scales by fusing protein images and interactions”](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1)**.
