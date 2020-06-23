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
1. `python community_detection.py`
2. `bash \<outprefix\>.sh`
3. Hierarchy is written in file \<outprefix\>.louvain.ddot with specific protein assignment for each system available in file \<outprefix\>.louvain.termStats

*Note: maximum clique finding is NP-hard, although utilized a heuristic approach, CliXO can still take a long time to finish.*

## Webpage
https://nrnb.org/music/

## Citing MuSIC
Please cite **[Qin et al., “Mapping cell structure across scales by fusing protein images and interactions”](https://www.biorxiv.org/cgi/content/short/2020.06.21.163709v1)**.
