function Z_cell= prox_largestK(X_cell, l2, lambdaG)
Z_cell= cell(size(X_cell));
for i=1:length(X_cell)
    Xi= X_cell{i};
    Yi= Xi - (lambdaG/l2);
    Yi(Yi<0)=0;
    Z_cell{i}= Yi;
end
end