function [S_new, Z_new] = LS_1_Algorithm(Rcell, S_cell, Y_cell, Ws, thetaCell, param, c_iter)
lmin1=1e-3; lmax1=1e3;
lmin2=1e-3; lmax2=1e3;
eta= param.eta;
delta= param.delta;
MAX_LS=20;

l1= 1; 
l2= 1;

S_candidate= S_cell;
Z_candidate= Y_cell;

for count=1:MAX_LS
    gradS_val= gradS_f(Rcell, S_candidate, Y_cell, Ws, thetaCell, param);
    S_candidate= update_S_cell(S_candidate, gradS_val, l1);
    
    gradG_val= gradG_f(Rcell, S_candidate, Y_cell, Ws, thetaCell, param);
    W_sub= subgrad_g(Y_cell, param.k_top);
    tmp = cell(size(Y_cell));
    for i=1:numel(Y_cell)
        tmp{i} = Y_cell{i} - (1/l2)*gradG_val{i};
    end
    
    % 再做: [上述结果] + (1/l2)*W_sub
    G_tilde = cell(size(Y_cell));
    for i=1:numel(Y_cell)
        G_tilde{i} = tmp{i} + (1/l2)*W_sub{i};
    end

    Z_candidate= prox_largestK(G_tilde, l2, param.lambdaG);
    F_val= computeF(Rcell, Z_candidate, S_candidate, Ws, thetaCell, param);

    dS = 0;
    for ii=1:length(S_candidate)
        D = S_candidate{ii} - S_cell{ii};
        dS = dS + sum(D(:).^2);
    end

    dG = 0;
    for ii=1:length(Z_candidate)
        D = Z_candidate{ii} - Y_cell{ii};
        dG = dG + sum(D(:).^2);
    end

    if F_val <= c_iter - delta*(dS + dG)
        break;
    end
    l1= min(max(l1*eta, lmin1), lmax1);
    l2= min(max(l2*eta, lmin2), lmax2);
end
S_new= S_candidate;
Z_new= Z_candidate;
end