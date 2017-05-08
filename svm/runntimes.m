function result = runntimes(K,lk,n,dynFeats, normal)
% Copyright 2012 Nino Shervashidze, Karsten Borgwardt
% edited by Leonardo Gutierrez, leonardo.gutierrez@uclouvain.be
% Input: K  - m x m kernel matrix 
%             if dynFeats = 1, K is a m x k matrix with m examples of dimension k 
%        lk - m x 1 array of class labels
%        n - number of times we want to run svm
%        dynFeats - boolean: 1 if K is a feature vector matrix, 0 if K is a
%        kernel matrix
%        normal -  boolean: 1 it will normalize feature vectors before
%        training the svm, 0 if not.
% Output: result - a structure with fields accuracy, mean, std and
%                  mean, std and se are the mean, the standard
%                  deviation and the standard error of accuracy
%
% Example: runntimes(feats,lmutag,10,1,0)
%
rand('seed', 666)
accuracy = zeros(n,1);

for i = 1:n
 
 if dynFeats
    [junk1, accuracy(i), junk2] = runIndependent_dynF(K, lk,normal)
    
 else
    [junk1, accuracy(i), junk2] = runIndependent(K, lk)
    
 end    

end

result.mean = mean(accuracy)
result.std = std(accuracy)
result.se = result.std / sqrt(n)
result.accuracy = accuracy

result
