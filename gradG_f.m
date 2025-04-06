function gradG= gradG_f(Rcell, S_cell_like, G_cell, Ws, thetaCell, param)
nTypes= param.nTypes;
instIdx= param.instanseIdx;
gradG= cell(nTypes,1);
for i=1:nTypes
    gradG{i}= zeros(size(G_cell{i}));
end

for rr=1:length(S_cell_like)
    if isempty(Rcell{rr}), continue; end
    w_ij= Ws(rr);
    [i,j]= deal(instIdx{rr}(1), instIdx{rr}(2));
    E_ij= G_cell{i}* S_cell_like{rr}* G_cell{j}' - Rcell{rr};

    grad_i= 2*w_ij* E_ij * ( G_cell{j}* S_cell_like{rr}' );
    grad_j= 2*w_ij* E_ij'*( G_cell{i}* S_cell_like{rr} );

    gradG{i}= gradG{i} + grad_i;
    gradG{j}= gradG{j} + grad_j;
end

lambda_ = param.lambda;
for i=1:nTypes
    if ~isempty(thetaCell{i})
        gradG{i}= gradG{i} + 2*lambda_*( thetaCell{i}* G_cell{i} );
    end
end
end