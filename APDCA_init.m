clear all; clc;
funspath=[pwd,filesep,'data',filesep];
addpath(funspath);

%% Load Data
load('lncRNAMiA.mat');
R12 = lncMI;       
[nlnc, nmi] = size(lncMI);
load('lncRNAGene.mat');
R13 = LGasso;      
load('lncRNAGOs.mat');
R14 = [lncBPs lncCCs lncMFs];
load('LncDOs.mat');
load('lncDisease2.mat');
LncDO = lncDisease + LncCancer;
LncDO(LncDO>1)=1;
R15 = LncDO;       
load('MiDOs.mat');
R25 = miDOs;       
load('miRNAGene.mat');
R23 = MGasso;      
load('GeneDisease.mat');
R35 = GDasso;      
[npro, ndi] = size(GDasso);
load('HumanGOAs.mat');
R34 = [bpLabels ccLabels mfLabels];
nGO = size(R34,2);
load('GeneDrug.mat');
R63 = GDrgasso;
nDrug = size(R63,1);

%% Filter disease
fun_stat = sum(R15,1);
sel_do_idx = find(fun_stat>0);
R15 = R15(:, sel_do_idx);
R25 = R25(:, sel_do_idx);
R35 = R35(:, sel_do_idx);
ndi = length(sel_do_idx); 

Rcell = {R12,R13,R14,R15,R23,R25,R34,R35,R63};

instanseIdx = {[1, 2],[1, 3],[1, 4],[1, 5],[2, 3],[2, 5],[3, 4],[3, 5],[6, 3]};

k1=50; k2=110; k3=50; k4=70; k5=170; k6=50;
R1 = [R12,R13,R14];      
R5_combined = [R25',R35']; 

[G1,~] = learnVector2(R1,k1);
[G2,~] = learnVector2(R25,k2);
[G3,~] = learnVector2(R35,k3);
[G4,~] = learnVector2(R14',k4);
[G5,~] = learnVector2(R5_combined,k5);
[G6,~] = learnVector2(R63,k6);

Gcell_init = {G1,G2,G3,G4,G5,G6};  

load('HumanPPI.mat');
load('DrugInter.mat');
thetaCell=cell(6,1);
thetaCell{3} = -PPI;
thetaCell{6} = -DDI;

params.maxIter    = 20;
params.delta      = 1e-5;
params.tau        = 0.4;
params.alpha      = 1e-5;
params.lambda     = 1e-2;
params.lambdaG    = 1e-5;
params.nTypes     = 6;
params.instanseIdx = instanseIdx;
params.k_top      = 100;
params.deltaLine  = 1e-3;
params.eta        = 1.01;
params.sigma      = 1e-4;
params.GcellInit  = Gcell_init;

tic;
[S_final, G_final] = APDCA_Demo(Rcell, thetaCell, params);
new_F = G_final{1} * S_final{4} * G_final{5}';
toc;
