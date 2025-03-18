%% 
clc; 
initCobraToolbox(false)
changeCobrasolver('gurobi','ALL');

%% LOADING MODELS
% Run this section from the main repository folder

% load('./models/toymodel1_bio.mat')
load('./models/bmodel.mat')
% load('./models/iML1515.mat')
% load('./models/PD.mat'); PD = model;
% load('./models/LS.mat'); LS = model;
% load('./models/yeast-GEM.mat'); yeast8 = model; clear model

%% SETTING PARAMETERS
% metaCone can be run with several different parameters.
% Make sure you call the function according to your needs.
% This time, want to quickly check all unfiltered solutions for bmodel

Alpha         = 0.2; % fraction related to the maximum growth according to FBA
Modality      = 'fast';
Exchanges     = findRxnIDs(bmodel,bmodel.rxns(contains(bmodel.rxns,'EX_')));
biomassIDX    = 13;
keepAll       = true;
eTol          = 1e-7;
vTol          = 1e-8;

%% RUNNING MetaCone

[CC, R] = metaCone(bmodel,'Modality',Modality,...
    'Alpha',Alpha,...
    'vTol',vTol,...
    'eTol',eTol,...
    'keepAll',keepAll,...
    'biomassIndex',biomassIDX);
fprintf('\nRESULT: %i conversions found\n', R.noconv)

%% SHOWING THE RESULTS

disp([R.exchanges array2table(CC)])

% Heatmap of the conversions
figure(1)
heatmap(CC)

% We will display one conversion as (macro-) chemical equation
cci = CC(:, randi(size(CC,2)));
mets = extractBetween(R.exchanges.RxnName(cci ~= 0), 'EX_', '(e)');
signes = sign(cci);
coeffs = cci(cci ~= 0);
type = signes(cci ~= 0);
substrates = mets(type < 0);
products = mets(type > 0);
coeffs_subs = (-1)*coeffs(type < 0);
coeffs_prod = coeffs(type > 0);
csubs = cellfun(@num2str,num2cell(coeffs_subs),'UniformOutput',false);
cprod = cellfun(@num2str,num2cell(coeffs_prod),'UniformOutput',false);

leftSide = join(join([csubs substrates])', ' + ');
rightSide = join(join([cprod products])', ' + ');
equation = join([leftSide rightSide],' -> ');

disp(equation)

%% ANALYZING 

% Singular Value Decomposition
[U, S] = svd(full(R.Allsol));

figure(2)
plot(1:min(size(S)), log10(diag(S)))
xlabel('Ordered (descending) indices of singular values')
ylabel('Log_{10} singular values')

% The logarithm of the singular values, when considering all the solutions
% that metaCone can find, shows the at which number the singular values
% drop to zero.
