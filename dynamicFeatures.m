% Compute dynamics feature vectors for a set of graphs
% Author: Leonardo Gutierrez Gomez - leonardo.gutierrez@uclouvain.be
% Copyright 2017 - Leonardo Gutierrez Gomez
% Input: Graphs - a 1xN array of graphs
% 		  Graphs(i).am is the adjacency matrix of the i'th graph, 
% 		  Graphs(i).al is the adjacency list of the i'th graph, 
%         Graphs(i).nl.values is a column vector of node
%                   labels for the i'th graph. When graphs are unlabeled, 
%                   this field is not used.
%         Graphs(i) may have other fields, but they will not be used here.
%         labels - a 1xN array of graph's labels
%         nodel -  a boolean: 1 if we want to use original node labels, 0 otherwise
%         l - lag value for stability computation. Recommended lag=3
% Output: data_feats - a NxK matrix in which each row corresponds to a K
%               dimensional feature vector
%
% Example: dynamicFeatures(mutag,lmutag,1,3)
%

function [data_feats] = dynamicFeatures(Graphs,labels,nodel,l)

py.importlib.import_module('networkx');

total = size(Graphs,2);
if nodel==1
  keys = node_labels(Graphs);
  vals = 1:size(keys,2);
  nodelabs = containers.Map(keys,vals);
end
for i=1:total
    
    N = size(Graphs(i).am,1);
    A = Graphs(i).am;
    if issparse(A)
        A = full(A);
    end
    
    one = ones(N,1);
        
    G = create_graph(A);
    d = A*one;
    D = diag(d);
    m = double(G.number_of_edges());  
    
    M = inv(D)*A;
    pi = d'/(2*m);
    display(i)
   
    lag = l;        
    H = eye(N);
    [Cov, Aut] = covariance(M,pi,H,lag);  

    [V D]=eigs(M');

    pi2 = V(:,2);
    h2 =  pi2/norm(pi2,1); % second eigvector
    
    res1 = stability(Aut);
    res2 = attribute_covariance(Cov,pi',lag);
    res3 = attribute_covariance(Cov,h2,lag);
    
    if nodel==1
        H_labels = compute_H(Graphs(i).nl.values);       
        [Cov, Aut] = covariance(M,pi,H_labels,lag);
        res4 = stability(Aut);
        r4(i,:) = res4; % Node labels
    end
    
    
    res7 = attribute_covariance(Cov,d,lag);
    
    tmp2 = py.networkx.betweenness_centrality(G);
    res8 = attribute_covariance(Cov,to_matlab(tmp2),lag);
    
    tmp3 = py.networkx.clustering(G);
    res9 = attribute_covariance(Cov,to_matlab(tmp3),lag);
    
    tmp4 = py.networkx.closeness_centrality(G);
    res10 = attribute_covariance(Cov,to_matlab(tmp4),lag);
    
    tmp5 = py.networkx.degree_centrality(G);
    res11 = attribute_covariance(Cov,to_matlab(tmp5),lag);
    
    res12 = py.networkx.degree_pearson_correlation_coefficient(G); % assortativity
    tmp5 = py.networkx.triangles(G);
    res13 = sum(to_matlab(tmp5))/3; % num de triangles
    res14 = py.networkx.transitivity(G); % global clustering coefficient
    
    res15 = res3*res7';
    res16 = res3*res8' ;
          
    r1(i,:) = res1; % H=identity
    r2(i,:) = res2; % pagerank
    r3(i,:) = res3; % second eigenvector    
    r5(i,:) = N;    % Num nodes
    r6(i,:) = m;    % num edges
    r7(i,:) = res7; % degree
    r8(i,:) = res8; % betweeness
    r9(i,:) = res9; % clustering coeff
    r10(i,:) = res10; % closeness centrality
    r11(i,:) = res11; % degree centrality
    r12(i,:) = res12; % asortativity
    r13(i,:) = res13; % num triangles
    r14(i,:) = res14; % global clustering coefficient   

end;
if nodel==1
    feats = struct('f1',r1,'f2',r2,'f3',r3,'f4',r4,'f5',r5,'f6',r6,'f7',r7,'f8',r8,'f9',r9,'f10',r10,'f11',r11,'f12',r12,'f13',r13,'f14',r14);
else
    feats = struct('f1',r1,'f2',r2,'f3',r3,'f5',r5,'f6',r6,'f7',r7,'f8',r8,'f9',r9,'f10',r10,'f11',r11,'f12',r12,'f13',r13,'f14',r14);
end

data_feats =struct('d', feats, 'l',labels);

save('features/dynamic_feats.mat','data_feats');



end

function [G] = create_graph(A)
    G = py.networkx.Graph;
    for i=1:size(A,1)
        for j=1:size(A,1)
            if A(i,j)>=1
                G.add_edge(i,j);
            end
        end
    end

end

function [res] = to_matlab(pylist)
    tmp = cell(pylist.values());
    res = cell2mat(tmp);
    res = double(res');
    
end
function [Cov, Aut] = covariance(M,pi, H,lag)

    PI = diag(pi);
    
    for t = 0:lag
        C = PI*M^t - pi'*pi;
        R = H'*C*H;
        Cov{t+1} = C;
        Aut{t+1} = R;
    end    

end

function [stab] = stability(Aut)

    for j=1:size(Aut,2)
       stab(j) = trace(Aut{j}); 
    end
end


function [u] = attribute_covariance(Cov,v,lag)
    u = zeros(1,lag+1);
    for t = 0:lag
        u(t+1) = v'*Cov{t+1}*v;
    end    

end


function [H] = compute_H(attribs) 
    num_nodes = size(attribs,1);
    keys = unique(attribs);
    num_classes = size(keys,1);
    vals = 1:num_classes;
    map = containers.Map(keys,vals);
    
    H = zeros(num_nodes,num_classes);
    for i=1:num_nodes      
        H(i,map(attribs(i))) = 1;
    end;    
end

function [nodelabs] = node_labels(Graphs)
    nodelabs = []
    for i=1:size(Graphs,2)
       vals = Graphs(i).nl.values;
       res = unique(vals,'stable')';
       nodelabs = [nodelabs res];
       
    end
    nodelabs = unique(nodelabs,'first');
end
