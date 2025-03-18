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

Alpha         = 0.1; % fraction related to the maximum growth according to FBA
Modality      = 'fast';
eTol          = 1e-7;
vTol          = 1e-8;

%% RUNNING QModel (conversions from scratch)
clc

models = {PD, LS};
Names  = split('PD LS');
QM = QModelCreation(models, Names, ...
                    'Alpha',Alpha, ...
                    'Modality',Modality,...
                    'eTol',eTol,...
                    'vTol',vTol);
disp(QM)
figure(1)
spy(QM.S); title("Q-Matrix of the PD-LS Q-Model")

%% RUNNING QModel (pre-calculated conversions)

% We obtain the conversions first.
[ccD, rd] = metaCone(PD);
[ccS, rs] = metaCone(LS);

% We set up the structures necessary for the function.
CCs    = {ccD , rd;...
          ccS , rs};
models = {PD; LS};
Names  = {'PD'; 'LS'};

% Calling QModelCreation
QM     = QModelCreation(models, Names, 'preCC', CCs,...
                                     'Alpha', Alpha,...
                                     'Modality', Modality);
fprintf('Q-Model between %s and %s finished\n', Names{1},Names{2})

%% TESTING THE MODEL
% FBA

SCFA = {'EX_ac(e)';...
        'EX_but(e)'; ...
        'EX_lac_D(e)';...
        'EX_lac_L(e)';...
        'EX_succ(e)';...
        'EX_ppa(e)'};
QM_1 = changeObjective(QM, SCFA);
sol  = optimizeCbModel(QM_1);
