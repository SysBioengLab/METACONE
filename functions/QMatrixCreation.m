function [Q, Output] = QMatrixCreation(NamesAbb, CCs, Rs)
%% 
% The function calcaluates a stoichiometric matrix (Q-Matrix) based on the
% conversions in 'CCs', in turn calculated by the metaCone() function. The
% Q variable will be an M x R Matrix, with M being the cardinality of the 
% union of all exchangeable metabolites of each model, and R is the sum of
% the dimensons of each C_ext basis. 
% Given each C_ext (output of metaCone() function), this function will
% paste column-wise the conversion, while smartly ordering the rows of the
% final matrix, so that the metabolites are alphabetically sorted.
%
%
% USAGE:
%
%   [Q, Output] = QMatrixCreation(Models, CCs, Rs)
%
% INPUTS:
%
%   NamesAbb        : Cell array of COBRA models, each with the fields:
%       * S         : m x n Stoichiometric Matrix
%       * c         : n x 1 linear objective coefficients
%       * lb        : n x 1 lower bounds on flux vector (the variables)
%       * ub        : n x 1 upper bounds on flux vector
%   CCs             : Cell array of 
%
% OUTPUTS:
%
%       * Q         : The Q-Matrix.
%       * Output    : A structure with the following fields:
%==========================================================================
%% PARSING INPUTS ===
%
p = inputParser;
p.KeepUnmatched = false;


% Default values ---


% Check Functions ---
isNumericMatrix     = @(M) and(isnumeric(M), numel(M)>1);
isCCsOk             = @(C) and(iscell(C), all(cellfun(isNumericMatrix,C)));


% Parameters ---
addRequired(p, 'NamesAbb'  , @iscellstr);
addRequired(p, 'CCs'       , isCCsOk);
addRequired(p, 'Rs'        , @iscell);

parse(p, NamesAbb, CCs, Rs)


%% INITIALIZATION ===


% Arguments Extraction and Initialization
% S        = full(p.Results.model.S);
% lb       = p.Results.model.lb;
% 
% 
% % Optional Calculation of the Elementary Matrix ---
% 
% % Variables Initialization ---
% noexch              = length(ExRxnIDs);
% 
% 
% %% GREEDY LP ITERATIONS ====
% 
% 
% %% OUTPUT ===
% 
% 
% % Additional information
% Output.runtime   = runtime;

end % of MetaCone function
