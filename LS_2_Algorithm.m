function [S_new, V_new] = LS_2_Algorithm(Rcell, S_cell, G_cell, Ws, thetaCell, param, c_iter)
lmin1=1e-3; lmax1=1e3;
lmin2=1e-3; lmax2=1e3;
eta= param.eta;
delta= param.delta;
MAX_LS=20;

l1= 1; 
l2= 1;

S_candidate= S_cell;
V_candidate= G_cell;

for count=1:MAX_LS
    gradS_val= gradS_f(Rcell, S_candidate, G_cell, Ws, thetaCell, param);
    S_candidate= update_S_cell(S_candidate, gradS_val, l1);

    gradG_val= gradG_f(Rcell, S_candidate, G_cell, Ws, thetaCell, param);
    W_sub= subgrad_g(G_cell, param.k_top);

    tmp = cell(size(G_cell));
    for i=1:numel(G_cell)
        tmp{i} = G_cell{i} - (1/l2)*gradG_val{i};
    end
    G_tilde = cell(size(G_cell));
    for i=1:numel(G_cell)
        G_tilde{i} = tmp{i} + (1/l2)*W_sub{i};
    end

    V_candidate= prox_largestK(G_tilde, l2, param.lambdaG);
    F_val= computeF(Rcell, V_candidate, S_candidate, Ws, thetaCell, param);
    
    dS = 0;
    for ii=1:length(S_candidate)
        D = S_candidate{ii} - S_cell{ii};
        dS = dS + sum(D(:).^2);
    end
    dG = 0;
    for ii=1:length(V_candidate)
        D = V_candidate{ii} - G_cell{ii};
        dG = dG + sum(D(:).^2);
    end

    if F_val <= c_iter - delta * ( dS + dG )
        break;
    end
    l1= min(max(l1*eta, lmin1), lmax1);
    l2= min(max(l2*eta, lmin2), lmax2);
end
S_new= S_candidate;
V_new= V_candidate;
end

