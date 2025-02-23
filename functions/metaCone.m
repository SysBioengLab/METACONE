function [CC, Output] = metaCone(model, varargin)
%% 
%MetaCon is a greedy LP-optimization based algorithm that aims to capture
% the metabolic potential of a COBRA model via basis calculation of a
% 'growing' Conversion Cone (CC); every conversion has a non-zero growth rate.
% The output is an MxN matrix, with M exchange reactions and N conversion
% vectors.
% The algorithm will keep track of the runtime, the rank of the matrix, the
% estimated nullity of the 'exchange' elementary matrix, the exhange
% reactions, and the elementary matrix itself.
% Optional inputs should be fed to the function as key-value pairs.
% The function will detect the exchange reactions, if not provided. 
% The function will select the first biomass rxn found, if not provided.
%
% USAGE:
%
%   [CC, output] = setCCon3_parfor(model)
%
% INPUTS:
%
%   model           : COBRA model (structure) with the following required fields:
%       * S         : m x n Stoichiometric Matrix
%       * c         : n x 1 linear objective coefficients
%       * lb        : n x 1 lower bounds on flux vector (the variables)
%       * ub        : n x 1 upper bounds on flux vector
%
% OPTIONAL INPUTS:
%
%       * Alpha     : a double between 0 and 1. 
%       * Modality  : 'full' or 'fast'
%       * biomassIndex  : biomass flux exchange index (integer)
%       * Exchanges : array (vector) with the exchanges (indices) of the model
%       * vTol      : zero flux tolerance (double)
%       * eTol      : zero epsilon tolerance (double)
%       * Nullity   : false (default) to not calculate elementary matrix
%
% OUTPUTS:
%
%       * Output    : A structure with the following fields:
%==========================================================================
%% PARSING INPUTS ===
%
p = inputParser;
p.KeepUnmatched = false;

% Default values ---
default_vTol        = 1e-8; 
default_eTol        = 1e-10; %1e-7 works for ecoli core
default_Alpha       = 0.1;
default_Nullity     = false;
default_bioIDX      = 0;
default_Exchanges   = true;
default_Modality    = 'fast';

% Check Functions ---
isvTolOk            = @(vTol) and(vTol < .5, vTol > 0);
isAlphaOk           = @(a) and(isnumeric(a), a >= 0 & a <= 1);
isNullityvalid      = @(nullity) ismember(nullity,[true, false]);
isbioIDXvalid       = @(bioIDX) class(bioIDX) == "double";
isExchsvalid        = @(exchs) or(iscellstr(exchs),isnumeric(exchs));
isModalityvalid     = @(mod) ismember(mod, {'full','fast'});
%isCobraModel as a subroutine

% Parameters ---
addRequired (p, 'model'     , @isCobraModel);
addParameter(p, 'Exchanges' , default_Exchanges, isExchsvalid);
addParameter(p, 'vTol'      , default_vTol,      isvTolOk);
addParameter(p, 'eTol'      , default_eTol,      isvTolOk);
addParameter(p, 'Alpha'     , default_Alpha,     isAlphaOk);
addParameter(p, 'Nullity'   , default_Nullity,   isNullityvalid);
addParameter(p, 'biomassIndex', default_bioIDX,  isbioIDXvalid);
addParameter(p, 'Modality'  , default_Modality,  isModalityvalid);

parse(p, model, varargin{:})

disp(p.Results)

%% INITIALIZATION ===
%

% Arguments Extraction and Initialization
S        = full(p.Results.model.S);
lb       = p.Results.model.lb;
ub       = p.Results.model.ub;
vTol     = p.Results.vTol;
eTol     = p.Results.eTol;
Alpha    = p.Results.Alpha;
Nullity  = p.Results.Nullity;
bioIDX   = p.Results.biomassIndex;
ExRxns   = p.Results.Exchanges;
Modality = p.Results.Modality;

% Proper exchanges extraction ---
switch class(ExRxns)
    case 'logical' % default
        ExRxnIDs = extractExchanges(model);
    case 'cell'
        ExRxnIDs = findRxnIDs(model, ExRxns);
    case 'double'
        ExRxnIDs = ExRxns;
end

% Finding Biomass reaction ---
if bioIDX == 0
    bioIDX = find(contains(model.rxnNames, 'biomass','IgnoreCase',true));
    bioIDX = bioIDX(1);
end

% Optional Calculation of the Elementary Matrix ---
if Nullity
    try
        [ME, EMCons] = buildElementaryMatrixCons(model, ExRxnIDs);
        nullity      = size(null(EMCons(:,ExRxnIDs)),2);
    catch
        ErrorMessage1 = "It was not possible to calculate the Elementary Matrix";
        disp(ErrorMessage1)
        nullity      = NaN;
        ME           = ErrorMessage1;
        EMCons       = NaN;
    end
else % default
    nullity = NaN;
end

% Variables Initialization ---
noexch   = length(ExRxnIDs);
% epsilon  = [];
% P_N      = 0; % Initial value of Projection matrix
FBA_init = optimizeCbModel(changeObjective(model, model.rxns(bioIDX)));
% FBA_init.f
% minBasis = []; % Solutions will be located as columns of this matrix
minBasis = FBA_init.x(ExRxnIDs); % Solutions will be located as columns of this matrix
P_N      = orth(minBasis)';
maxg = Alpha*FBA_init.f; % minimum growth per conversion
if Modality == "fast"
    % w = randn(noexch,1);
    w = randi(1e5, [1,noexch]) + rand(1,noexch); 
    % w = 1 + rand(1, noexch);
end


% disp(ExRxnIDs)
% disp(bioIDX)
% disp(nullity)
% disp(FBA_init)
disp(table(model.rxns(ExRxnIDs), minBasis))
disp(P_N)
disp(P_N*P_N')
fprintf('minimum growth: %i\n', maxg)
fprintf('Nº of exchages: %i\n', noexch)
disp(w)

CC = ExRxnIDs;
Output = bioIDX;

end

%==========================================================================
%% COMPLEMENTARY SUBROUTINES
% Subroutine to check if 'model' arg has proper fields ====================
function output = isCobraModel(model) 
if ~isstruct(model)
    output = false;
else
    if ~and(isfield(model,'S'), isnumeric(model.S))
        output = false; 
    else
        if and(isfield(model,'lb'),isfield(model,'ub'))
            if and(isnumeric(model.lb), isnumeric(model.ub))
                if and(isfield(model, 'c'), isnumeric(model.c))
                    output = true;
                end
            end
        else
            output = false;
        end
    end
end
end

% Subroutine to detect exchange reactions =================================
function ExRxnIDs = extractExchanges(model)
%Beta version 

% Option 1
ExRxns = findExcRxns(model);
ExRxnIDs = findRxnIDs(model, model.rxns(ExRxns));

% % Option 2
% S = model.S;
% ExRxnIDs = find((sum(S == 0)==(n_row-1) & sum(S == -1) == 1)); %Extract position of exchanges
% 
% % Option 3
% ExRxnIDs = find(contains(model.rxns, 'EX_'));
end

% Subroutine to calculate Elementary Matrix with the exchange reactions ===
function [ME, EMCons] = buildElementaryMatrixCons(model, ExRxnsIDs)
S = full(model.S);

%ExRxns = model.rxns(ExRxnsIDs);
Emets = findMetsFromRxns(model, model.rxns(ExRxnsIDs));
[ME, elem] = computeElementalMatrix(model, Emets);
disp("The elementary matrix was calculated for the following elements:")
disp(elem) %
EMCons = zeros(length(elem),length(model.rxns)); %
EmetsIDs = findMetIDs(model, Emets);
ME_t = ME';
%Aqui falta un chequeo de que EX_rxns:EmetsIDs = 1:1
ME_rearr = ME_t(:,(1:length(Emets))*logical(S(EmetsIDs, ExRxnIDs)));  %reordena la matriz elemental
EMCons(:,ExRxnIDs) = ME_rearr;
end