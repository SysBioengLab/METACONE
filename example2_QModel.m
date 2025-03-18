%% 
clc; 
initCobraToolbox(false)
changeCobrasolver('gurobi','ALL');

%% LOADING MODELS
% Run this section from the main repository folder

% load('./models/toymodel1_bio.mat')
% load('./models/bmodel.mat')
load('./models/PD.mat'); PD = model;
load('./models/LS.mat'); LS = model;
% load('./models/yeast-GEM.mat'); yeast8 = model; clear model

%% SETTING PARAMETERS
% metaCone can be run with several different parameters, and so, QModel
% should receive the same arguments, or 
% Make sure you call the function according to your needs.

Alpha         = 0.2; % fraction related to the maximum growth according to FBA
Modality      = 'fast';
Exchanges     = findRxnIDs(bmodel,bmodel.rxns(contains(bmodel.rxns,'EX_')));
biomassIDX    = 0;
keepAll       = true;
eTol          = 1e-7;
vTol          = 1e-8;

%% RUNNING QModel

%% SHOWING THE RESULTS

%% ANALYZING SOLUTIONS
