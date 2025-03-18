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
ax = gca;
ax.Colormap = colormap_dorei_symb;
ax.YDisplayLabels = extractBetween(R.exchanges.RxnName, 'EX_', '(e)');
xlabel("Conversion Index")
ylabel("Metabolite")
title("Heatmap showing production/consumption fashion")

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

%% ANALYZING SOLUTIONS

% Norm of Epsilon vectors
norms = full(vecnorm(R.Epsilons));

% Singular Value Decomposition
[U, S] = svd(full(R.Allsol));

% Plotting the results
figure(2); clf
plot(1:min(size(S)), log10(diag(S)))
xlabel('Ascending indices of singular values or Epsilons vectors')
ylabel('Log_{10} singular values')
hold on
yyaxis right
plot(1:numel(norms), norms)
ylabel("2-norm of Epsilon vectors")
title("Analysis of all solutions")
xlim([0 22])

% The logarithm of the singular values, when considering all the solutions
% that metaCone can find, shows the at which number the singular values
% drop to zero.

% Each solution comes with an Epsilon vector which norms should be as close
% to zero as possible.
