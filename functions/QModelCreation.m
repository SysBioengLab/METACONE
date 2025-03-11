function QModel = QModelCreation(models, Names, varargin)
%%
% QModelCreation generates a structure with the framework of a COBRA model,
% using the Q Matrix as the stoichiometric Matrix (S), which is composed of
% metabolic l.i. conversions of a Microorganism (MO) as columns, plus exchange
% reactions, and the rows as metabolite balances. 
% The 2 cores of the function are the metacone() algorithm and the
% QMatrixCreation() function. 
% The input cell arrays are paired element-wise, ensuring a one-to-one 
% correspondence between their elements in the output.
% The same arguments for metaCone can be passed down to this function.
% onsider that all C_ext will be calculated with the same parameters, if
% this is done.
%
% USAGE:
%
%   [QModel] = QModelCreation(models, Names, varargin)
%
% INPUTS:
%
%   models          : An array of COBRA models with the required fields:
%       * S         : m x n Stoichiometric Matrix
%       * c         : n x 1 linear objective coefficients
%       * lb        : n x 1 lower bounds on flux vector (the variables)
%       * ub        : n x 1 upper bounds on flux vector
%   Names           : A cell array with the names of the models.
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
%       * QModel    : A struct with the main fields of a COBRA Model
%       * CC        :
%       * Output    : A structure with the following fields:
%==========================================================================
%% PARSING INPUTS ===
%
p = inputParser;
p.KeepUnmatched = false;

% Default values ---
defaultCCs          = {};
default_vTol        = 1e-8; 
default_eTol        = 1e-10; %1e-7 works for ecoli core
default_Alpha       = 0.1;
default_Nullity     = false;
default_bioIDX      = 0;
% default_Exchanges   = true;
default_Modality    = 'fast';

% Check Functions ---
isModelsOk          = @(models) all(cellfun(@isCobraModel, models));
isAlphaOk           = @(a) and(isnumeric(a), a >= 0 & a <= 1);
isNullityvalid      = @(nullity) ismember(nullity,[true, false]);
isbioIDXvalid       = @(bi) and(size(bi)==size(models), isnumeric(bi));
% isExchsvalid        = @(exchs) or(iscellstr(exchs),isnumeric(exchs));
isModalityvalid     = @(mod) ismember(mod, {'full','fast'});
isvTolok            = @(vTol) and(vTol < .1, vTol > 0);
isCCsvalid          = @(ccs) and(iscell(ccs),all(size(ccs)==[numel(models) 2]));

% Parameters ---
addRequired (p, 'models'        , isModelsOk);
addRequired (p, 'Names'         , @iscell);
addOptional (p, 'preCC'         , defaultCCs      , isCCsvalid);
addParameter(p, 'vTol'          , default_vTol    , isvTolok); 
addParameter(p, 'eTol'          , default_eTol , isvTolok); 
addParameter(p, 'Alpha'         , default_Alpha   , isAlphaOk);
addParameter(p, 'Nullity'       , default_Nullity , isNullityvalid); 
addParameter(p, 'biomassIndex'  , default_bioIDX  , isbioIDXvalid);
addParameter(p, 'Modality'      , default_Modality, isModalityvalid);
% addParameter(p, 'Exchanges'     , default_Exchanges, isExchsvalid);

parse(p,models,Names,varargin{:})

fprintf('Call: \n')
disp(p.Results)

%% INITIALIZATION ===


% Arguments Extraction and Initialization
vTol     = p.Results.vTol;
eTol     = p.Results.eTol;
Alpha    = p.Results.Alpha;
Nullity  = p.Results.Nullity;
bioIDX   = p.Results.biomassIndex;
Modality = p.Results.Modality;

end
