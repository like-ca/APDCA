function S_cell= updateS_all(Rcell, Gcell, instIdx)
S_cell= cell(length(instIdx),1);
for rr=1:length(instIdx)
    [i,j]= deal(instIdx{rr}(1), instIdx{rr}(2));
    Gi= Gcell{i}; 
    Gj= Gcell{j}; 
    Rij= Rcell{rr};

    A= Gi'*Gi + 1e-12*eye(size(Gi,2));
    B= Gj'*Gj + 1e-12*eye(size(Gj,2));
    S_ij= pinv(A)*(Gi'*Rij*Gj)* pinv(B);
    S_ij(isnan(S_ij))=0;
    S_cell{rr}= S_ij;
end
end
