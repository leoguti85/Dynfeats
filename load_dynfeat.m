% Load dynamics feature vectors from a dataset
% Author: Leonardo Gutierrez Gomez - leonardo.gutierrez@uclouvain.be
% Copyright 2017 - Leonardo Gutierrez Gomez
% Note: This script works only with the output of dynamcFeatures.m 
% Input:                                  
% 		 feats: array of features to be loaded according with the following coding:
%                1 - H = identity
%                2 - pagerank
%                3 - second eigenvector
%                4 - node labels encoded in H
%                5 - number of nodes
%                6 - number of edges
%                7 - Nodes degrees
%                8 - Node betweeness
%                9 - Local clustering coefficient
%                10 - Closeness centrality        
%                11 - Degree centrality
%                12 - Assortativity
%                13 - Number of triangles
%                14 - Global clustering coefficient

% Output:  res - a Nxp matrix in which each row corresponds to a p
%                dimensional feature vector
%
% Example: Loading features H=identity and node labels
% load_dynfeat([1 4])
%

function [res] = load_dynfeat(feats)
    %path = strcat('/path/to/dynamicFeatures/',dataset,'/dynamic_feats.mat');
    data = load('features/dynamic_feats.mat');  
    
    data_feats = data.data_feats
    res1=[]; res2=res1; res3=res1;res4=res1;res5=res1;res6=res1;res7=res1;
    res8=res1;res9=res1;res10=res1;res11=res1;res12=res1;res13=res1;res14=res1;
    
    for opt=feats
        
        switch opt
            case 1 
                res1 = data_feats.d.f1;  % H=id           
            case 2
                  res2 = data_feats.d.f2 ; % pagerank          
            case 3
                  res3 = data_feats.d.f3;  % second eigenvector          
            case 4
                  res4 = data_feats.d.f4;  % node labels         
            case 5
                  res5 = data_feats.d.f5;  % Num nodes          
            case 6
                  res6 = data_feats.d.f6;  % num edges                        
            case 7
                  res7 = data_feats.d.f7;  % Degrees          
            case 8
                  res8 = data_feats.d.f8;  % Betweeness      
            case 9
                  res9 = data_feats.d.f9;  % Clustering coefficient            
            case 10
                  res10 = data_feats.d.f10;  % Closeness centrality                  
            case 11
                  res11 = data_feats.d.f11;  % degree centrality                  
            case 12
                  res12 = data_feats.d.f12;  % assortativity                
            case 13
                  res13 = data_feats.d.f13;  % num triangules  
            case 14
                  res14 = data_feats.d.f14;  % global clustering coefficient                        
                        
        end
            
        
    end
    res = [res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res13,res14];

end
