function [result,mean_accuracy,std_accuracy] = runIndependent2(K,lk,normal)
% Copyright 2012 Nino Shervashidze, Karsten Borgwardt
% Edited by Leonardo Gutierrez, leonardo.gutierrez@uclouvain.be
% K = matrix of n feature vectors of dimension k (n*k)
% lk = vector of labels (n*1)
% normal = boolean: 1 if normalize data before train svm, 0 if not

addpath('/path/libsvm-3.22/matlab/');

% randomly permute kernel matrix and labels
r = randperm(size(K,1));
K = K(r,:);

lk = lk(r);
lk=lk';
lkoriginal = lk; 
Koriginal = K;
num_examples = size(K,1)

%% stratified cross-validation
%sum(sum(Koriginal));
%neworder = stratifiedsplit(lk)
%for i = 1:size(neworder,2) 
%m = size(neworder(i).old,1);
%r = randperm(m);
%newlk(neworder(i).new) =  lk(neworder(i).old(r));
%Knew([neworder(i).new]',[neworder(i).new]') = K([neworder(i).old(r)]',[neworder(i).old(r)]');  
%end
%
%sum(sum(Knew)) - sum(sum(Koriginal))
%dbstop
%lk = newlk'
%K = Knew;
%size(lk);
%size(K);
%dbstop 
%Koriginal = K;


% bring kernel matrix into libsvm format
p80 = ceil(size(K,1) * 0.8);
p90 = ceil(size(K,1) * 0.9);

% specify range of c-values
cvalues = (10 .^ [-7:2:7]) / size(K,1);
%cvalues = [10^(-3) 10^(-2) 10^(-1) 10^(0) 10^1 10^2 10^3]

cv = 10;
fs = size(K,1) - p90;

% cross-validation loop
for k = 1:cv
K = Koriginal;
lk = lkoriginal;

K = K([k*fs+1:size(K,1),1:(k-1)*fs,(k-1)*fs+1:k*fs],:);  
lk = lk([k*fs+1:size(K,1),1:(k-1)*fs,(k-1)*fs+1:k*fs]); 

% preprocessing

if normal~=1
    Ktr_std = K(1:p80,:);
    Ktv_std = K(p80+1:p90,:);
else
    [Ktr_std,mu,sigma] = zscore_normalize_data(K(1:p80,:),mean(K(1:p80,:)),std(K(1:p80,:)));
    Ktv_std = zscore_normalize_data(K(p80+1:p90,:),mu,sigma);
end
%K = makepos(K);
%K1 = [(1:size(K,1))', normalizekm(K)];

%if any(strcmp('optimal',options))
imresult=[];
for i = 1:size(cvalues,2)
    % train on 80%, predict on 10% (from 81% to 90%)    
  model = svmtrain(lk(1:p80)', Ktr_std, strcat(['-t 0 -c ' num2str(cvalues(i))]));
  [predict_label, accuracy, dec_values] = svmpredict(lk(p80+1:p90)',Ktv_std, model);
  accuracy80 = accuracy;
  imresult(i)= accuracy(1);
end
  
  if normal~=1
      Ktr_std = K(1:p90,:);
      Kte_std = K(p90+1:num_examples,:);
  else
      [Ktr_std, mu, sigma] = zscore_normalize_data(K(1:p90,:),mean(K(1:p90,:)),std(K(1:p90,:)));
      Kte_std = zscore_normalize_data(K(p90+1:num_examples,:),mu,sigma);
  end
  % determine optimal c
  [junk,optimalc]= max(fliplr(imresult));
  optimalc = size(cvalues,2)+1 - optimalc; 
  % train on 90% with optimal c, predict on 10% (from 91% to 100%)
  model = svmtrain(lk(1:p90)', Ktr_std,strcat(['-t 0 -c ' num2str(cvalues(optimalc))]) );
  [predict_label, accuracy, dec_values] = svmpredict(lk(p90+1:num_examples)', Kte_std, model);
  accuracy90 = accuracy;
  result(k)=accuracy(1);


end  
mean_accuracy =  mean(result) ;
std_accuracy = std(result);


end

%
%% cross-validation
%if any(strcmp('cv',options))
%options = strcat(['-t 4 -v ' num2str(cv) ' -c ' num2str(cvalues(i))])
%result(i) = svmtrain(lk, K1, options); %', num2str(cv)));
%end
%
%end

function [res] = preproc(data)
    
    nfeat = size(data,2);
    res = zeros(size(data));


    for i=1:nfeat
      xmin = min(data(:,i));
      xmax = max(data(:,i));
      res(:,i) = (data(:,i) - xmin)/(xmax - xmin);
    end

end

function [X] = standardize(varargin)
switch nargin
    case 1
        mean_X = mean(varargin{1}); %// Find mean of each column
        std_X = std(varargin{1}); %// Find std. dev. of each column

        X = bsxfun(@minus, varargin{1}, mean_X); %// Subtract each column by its respective mean
        X = bsxfun(@rdivide, X, std_X); %// Take each column and divide by its respective std dev.

    case 3
        mean_X = varargin{2};
        std_X = varargin{3};

        %// Same code as above
        X = bsxfun(@minus, varargin{1}, mean_X);
        X = bsxfun(@rdivide, X, std_X);
end
end

function result = makepos(K)
pd = 0;
addc = 10e-7;
while (pd ==  0)
  
  addc = addc * 10
  try
    if (isinf(addc) == 1)
      pd = 1;
    else 
      chol(normalizekm(K + eye(size(K,1),size(K,1)) * addc));
      pd = 1;
    end
  catch
    
  end
  
end
if (isinf(addc)==0)
  result = K + eye(size(K,1),size(K,1)) * addc;
else
  result = eye(size(K,1));
end
end
