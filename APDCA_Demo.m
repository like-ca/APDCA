function [S_cell_final, G_cell_final] = APDCA_Demo(Rcell, thetaCell, param)
instIdx = param.instanseIdx;
nTypes  = param.nTypes;
maxIter = param.maxIter;
deltaVal= param.delta;
tauVal  = param.tau;

S_cell= updateS_all(Rcell, param.GcellInit, instIdx);

mus= zeros(length(instIdx),1);
for rr=1:length(instIdx)
    [i,j]= deal(instIdx{rr}(1), instIdx{rr}(2));
    E_ij= Rcell{rr}- param.GcellInit{i}* S_cell{rr} * param.GcellInit{j}';
    mus(rr)= sum(E_ij(:).^2);
end
[~, Ws] = getOptimalWeights(mus, param.alpha);

thetaOld=0; 
thetaCur=1;
q_val=1;

F_val= computeF(Rcell, param.GcellInit, S_cell, Ws, thetaCell, param);
c_val= F_val;

Z_cell= param.GcellInit;
G_old= param.GcellInit;
Gcell= param.GcellInit;

for iter=1:maxIter
    fprintf('APDCA iter %d/%d\n', iter, maxIter);
    Y_cell= cell(nTypes,1);
    for i=1:nTypes
        Y_cell{i} = Gcell{i}  + (thetaOld/thetaCur)*(Z_cell{i}- Gcell{i}) + ((thetaOld-1)/thetaCur)*( Gcell{i} - G_old{i});
    end

    [S_next, Z_next] = LS_1_Algorithm(Rcell, S_cell, Y_cell, Ws, thetaCell, param, c_val);
    F_candidate= computeF(Rcell, Z_next, S_next, Ws, thetaCell, param);

    dS = 0;
    for ii=1:length(S_next)
        D = S_next{ii} - S_cell{ii};
        dS = dS + sum(D(:).^2);
    end
    dG = 0;
    for ii=1:length(Z_next)
        D = Z_next{ii} - Y_cell{ii};
        dG = dG + sum(D(:).^2);
    end

    if F_candidate <= c_val - deltaVal*(dS + dG)
        G_new= Z_next;
        S_new= S_next;
        finalF= F_candidate;
    else
        [S_next2, V_next] = LS_2_Algorithm(Rcell, S_cell, Gcell, Ws, thetaCell, param, c_val);
        F_candidate2= computeF(Rcell, V_next, S_next2, Ws, thetaCell, param);
        if F_candidate <= F_candidate2
            G_new= Z_next;
            S_new= S_next;
            finalF= F_candidate;
        else
            G_new= V_next;
            S_new= S_next2;
            finalF= F_candidate2;
        end
    end

    mus2= zeros(length(instIdx),1);
    for rr=1:length(instIdx)
        [i,j]= deal(instIdx{rr}(1), instIdx{rr}(2));
        E_ij= Rcell{rr}- G_new{i}* S_new{rr} * G_new{j}';
        mus2(rr)= sum(E_ij(:).^2);
    end
    [~, Ws_next] = getOptimalWeights(mus2, param.alpha);

    thetaNew= ( sqrt(4*thetaCur^2 +1)+1 )/2;
    q_next= tauVal*q_val +1;
    c_next= (tauVal*q_val*c_val + finalF)/ q_next;

    G_old= Gcell;
    Gcell= G_new;
    S_cell= S_new;
    Ws= Ws_next;

    thetaOld= thetaCur;
    thetaCur= thetaNew;
    q_val= q_next;
    c_val= c_next;
    F_val= finalF;
    
end
S_cell_final = S_cell;
G_cell_final = Gcell;
if sum(sum(G_cell_final{1})) == 0
    G_cell_final{1} = param.GcellInit{1};
end
end
