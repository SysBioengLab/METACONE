function [Q, Output] = QMatrixCreation(NamesAbb, CCs, Rs)
%% (Beta) Q-Matrix Creation Function
% The function calcaluates a stoichiometric matrix (Q-Matrix) based on the
% conversions in 'CCs', in turn calculated by the metaCone() function. The
% Q variable will be an M x R Matrix, with M being the cardinality of the 
% union of all exchangeable metabolites of each model, and R is the sum of
% the dimensons of each C_ext basis. 
% Given each C_ext (output of metaCone() function), this function will
% paste column-wise the conversions, while smartly ordering the rows of the
% final matrix, so that the metabolites are alphabetically sorted.
%
%
% USAGE:
%
%   [Q, Output] = QMatrixCreation(Models, CCs, Rs)
%
% INPUTS:
%
%   NamesAbb        : Cell array of COBRA models
%   CCs             : Cell array of C_ext matrices delivered by metaCone().
%   Rs              : Structure obtained from metaCone(). The original
% models must contain a 'EX_biomass_(e)' reaction in model.rxns
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
NamesAbb = p.Results.NamesAbb;
CCs      = p.Results.CCs;
Rs       = p.Results.Rs;


% % Variables Initialization ---
nums     = numel(NamesAbb);
allExcs  = {}; 
Q        = [];
Keys     = cell(1, nums);

%% COMMON EXCHANGES SET and NORMALIZATION ====
for i = 1:numel(CCs)
    % Annotating current cc and exchanges ----
    cc = CCs{i};
    exchanges = Rs{i}.exchanges.RxnName;
    
    % Templates for finding the biomass exchanges
    a = contains(exchanges, 'EX_biomass(e)');
    b = contains(exchanges, 'EX_biomass[e]');
    
    % Extracting the biomass exchange reaction ----
    if any(a)
        
        % We tag each biommas reaction with the name abb.
        exchanges(a) = insertAfter(exchanges(a),'mass',['_' NamesAbb{i}]);
        
        % We normalize by the flux value of the biomass reaction.
        cc = cc./cc(a,:);
        
    elseif any(b)
        exchanges(b) = insertAfter(exchanges(b),'mass',['_' NamesAbb{i}]);
        cc = cc./cc(b,:);
    else
        % Termination ........
        error('Some models have no biomass exchange reaction')
    end
    % Modifying variables
    Rs{i}.exchanges.RxnName = exchanges; 
    CCs{i} = cc; 
    
    % Creating the exchanges union set.
    allExcs = union(allExcs, exchanges);
end

%% KEYS FOR EACH SPECIES ====
for i = 1:nums
    exchanges = Rs{i}.exchanges.RxnName;
    tempInd   = ismember(allExcs, exchanges);
    partexc   = allExcs(tempInd);
    tempKey   = zeros(numel(exchanges),1);
    if ~isempty(setdiff(partexc, exchanges))
        error('FATAL error while ordering exchanges')
    end

    % For each member, we get the indices of their exchanges
    for j = 1:numel(partexc)
        tempKey(j) = find(contains(exchanges, partexc(j)));
    end
    Keys{i} = tempKey;
end; clear i

%% Q-MATRIX REARRANGE ====

% The Conversions are glued and the rows rearrenged according to the keys.
for i = 1:nums
    exchanges = Rs{i}.exchanges.RxnName;
    tempInd   = ismember(allExcs, exchanges);
    Qtemp = sparse(zeros(size(allExcs,1), Rs{i}.noconv));
    conversions = CCs{i};
    Qtemp(tempInd,:) = conversions(Keys{i},:);
    Q = [Q, Qtemp];
end

% Adding the community exchanges
% "easy approach": each metabolite has a secreting echange reaction of the
% form: 'met_(e) -> '
[m, ~] = size(Q);
Q = [Q, -eye(m)];

%% OUTPUT ===

Output.allExcs   = allExcs;
Output.Keys      = Keys;

end % of main function
