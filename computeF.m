function val= computeF(Rcell, G_cell, S_cell, Ws, thetaCell, param)
instIdx= param.instanseIdx;
lambda_ = param.lambda;
lambdaG = param.lambdaG;
k      = param.k_top;

val=0;

for rr=1:length(S_cell)
    if isempty(Rcell{rr}), continue; end
    w_ij= Ws(rr);
    [i,j]= deal(instIdx{rr}(1), instIdx{rr}(2));
    resid= Rcell{rr} - G_cell{i}* S_cell{rr}* G_cell{j}';
    val= val + w_ij* sum(resid(:).^2);
end

for i=1:length(G_cell)
    if ~isempty(thetaCell{i})
        val= val + lambda_* sum(sum( (G_cell{i}'*thetaCell{i}).* G_cell{i}' ));
    end
end

for i=1:length(G_cell)
    x= G_cell{i}(:);
    sum_all= sum(x);
    [xs, ~]= sort(x,'descend');
    topSum= sum(xs(1: min(k, end)));
    val= val + lambdaG*( sum_all - topSum );
end
end
