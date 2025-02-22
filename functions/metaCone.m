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
isAlphaOk           = @(a) and(isnumeric(a), a >= 0 | a <= 1);
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
end


%% COMPLEMENTARY SUBROUTINES
% Subroutine to check if 'model' arg has proper fields
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