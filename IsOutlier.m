function [isoutlier,Mdl,quartiles,f1,f2] = IsOutlier(x,y,varargin)

[ntrees] = myparse(varargin,'ntrees',200);

tau = [0.25 0.5 0.75];
Tbl = table(x,y);
Mdl = TreeBagger(ntrees,Tbl,'y','Method','regression');

quartiles = quantilePredict(Mdl,x,'Quantile',tau);

iqr = quartiles(:,3) - quartiles(:,1);
k = 1.5;
f1 = quartiles(:,1) - k*iqr;
f2 = quartiles(:,3) + k*iqr;

isoutlier = y < f1 | y > f2;