%% Loading models for METACONE analysis
% Models need to be loaded and adapted. Make sure that Cobra Toolbox
% has been already initialized.
% Move script to the folder with the models.

clc; clear
% ======================================================
% Change the path to the working directory.
cd /Users/alvin/Doctorado/TESIS/Article2/
% initCobraToolbox(false)
% changeCobraSolver('gurobi','ALL');
% ======================================================

%% Loading 

load('./models/toymodel1_bio.mat')
load('./models/bmodel.mat')
load('./models/iML1515.mat')
load('./models/Bacteroides_dorei_DSM_17855.mat'); BdD = model;
load('./models/Clostridium_symbiosum_WAL_14673.mat'); CsW = model;
load('./models/yeast-GEM.mat'); yeast8 = model; clear model

%% Models modification
% We will add inluin uptake and or consumption  to both models


% Xylan consumption
BdD = addReaction(BdD,'EX_xylan(e)', 'xylan[e] <=> ');
BdD = addReaction(BdD, 'XYLAN_DEGe', '527.0 h2o[e] + xylan[e] -> 528.0 xyl_D[e]');
CsW = addReaction(CsW,'XYLI1' , 'xyl_D[c]  <=> xylu_D[c]');

% Inulin consumption
BdD = addReaction(BdD, 'INULINASE', 'reactionFormula', '29.0 h2o[c] + inulin[c] -> 29.0 fru[c] + glc_D[c]');
BdD = addReaction(BdD, 'INULINabc', 'reactionFormula', 'atp[c] + h2o[c] + inulin[e] -> adp[c] + h[c] + inulin[c] + pi[c]');
BdD = addReaction(BdD, 'EX_inulin(e)', 'reactionFormula', 'inulin[e] <=> ');

CsW = addReaction(CsW, 'INULINASE', 'reactionFormula', '29.0 h2o[c] + inulin[c] -> 29.0 fru[c] + glc_D[c]');
CsW = addReaction(CsW, 'INULINabc', 'reactionFormula', 'atp[c] + h2o[c] + inulin[e] -> adp[c] + h[c] + inulin[c] + pi[c]');
CsW = addReaction(CsW, 'EX_inulin(e)', 'reactionFormula', 'inulin[e] <=> ');

%% We save the models as a cell array

Models = {toymodel1_bio, bmodel, BdD, CsW, iML1515, yeast8}';
Names = {'toym', 'ecore_bio', 'dore', 'symb', 'ecoli_full', 'yeast8'}';

clear ans BdD bmodel CsW gRates iML1515 toymodel1_bio yeast8
