%% 
clc; 
initCobraToolbox(false)
changeCobrasolver('gurobi','ALL');

%% LOADING MODELS
% Run this section from the main repository folder

load('./models/toymodel1_bio.mat')
load('./models/bmodel.mat')
load('./models/iML1515.mat')
load('./models/PD.mat'); PD = model;
load('./models/LS.mat'); LS = model;
load('./models/yeast-GEM.mat'); yeast8 = model; clear model

%% SETTING PARAMETERS
% metaCone can be run with several different parameters.
% Make sure you call the function according to your needs.
