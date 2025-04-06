function [newX, Z] = learnVector2(LabelMat, dim)
[U,S,V]= svd(LabelMat,'econ');
ds= diag(S);
rnk= min(dim, length(ds));

U= U(:,1:rnk);
S= diag(ds(1:rnk));
V= V(:,1:rnk);

X= U*sqrt(S);
Z= V*sqrt(S);

X= mapminmax(X',0,1)';
Z= mapminmax(Z',0,1)';
newX= X;
end