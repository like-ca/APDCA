function gradS= gradS_f(Rcell, S_cell, G_cell_like, Ws, thetaCell, param)
instIdx= param.instanseIdx;
gradS= cell(size(S_cell));
for rr=1:length(S_cell)
    [i,j]= deal(instIdx{rr}(1), instIdx{rr}(2));
    w_ij= Ws(rr);
    E_ij= G_cell_like{i}* S_cell{rr} * G_cell_like{j}' - Rcell{rr};

    gradS{rr}= 2*w_ij*( G_cell_like{i}'* E_ij * G_cell_like{j} );
end
end

