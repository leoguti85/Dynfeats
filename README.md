# Dynfeats

This repository provides information about the implementation of the method proposed in the paper:
Dynamics Based Features for Graph Classification

This is a feature based method for creating feature vector of graphs for classification purposes. More details can be found in the website

USAGE:

1) Loading dataset. ex:
load('datasets/mutag.mat')

2) Generating features. This script will generate all features. ex:
dynamicFeatures(mutag,lmutag,1,3)

3) Selecting features acording with the following codification:

                1 - H = identity
                2 - pagerank
                3 - second eigenvector
                4 - node labels encoded in H
                5 - number of nodes
                6 - number of edges
                7 - Nodes degrees
                8 - Node betweeness
                9 - Local clustering coefficient
                10 - Closeness centrality        
                11 - Degree centrality
                12 - Assortativity
                13 - Number of triangles
                14 - Global clustering coefficient

ex. feats = load_dynfeat([1 4])

4) Training classifier. ex:
result=runntimes(feats, lmutag,10,1,0)

Note: Do not forget to change your libsvm path into runIndependent_dynF.m file.

