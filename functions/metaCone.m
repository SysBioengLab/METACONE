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
% The 'fast'modality
% The 'full' modality performs 2J LP problems per iteration.
%
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
FBA_init = optimizeCbModel(changeObjective(model, model.rxns(bioIDX)));
minBasis = FBA_init.x(ExRxnIDs); % Solutions will be located as columns of this matrix
P_N      = orth(minBasis)'; % Initial value of Projection matrix
maxg = Alpha*FBA_init.f; % minimum growth per conversion
if Modality == "fast"
    % w = randn(noexch,1);
    w = randi(1e5, [1,noexch]) + rand(1,noexch); 
    % w = 1 + rand(1, noexch);
else % 'full'
    w = 1;
end

%% GREEDY LP ITERATIONS ====
fprintf('Starting %s modality of the algorithm: \n', Modality)
fprintf('Minimum growth per conversion set to %f: \n', maxg)
k = 0;

% Complementary Performance
nits = 29;
Epsilons = [];
allSols = [];
solspK = zeros(nits,1);
params.OutputFlag  = 0;
params.LPWarmStart = 1;

tic
while true
    % Projection Matrix update.
    P_NT = eye(noexch) - P_N'*P_N;
    
    EpsilonsK_max = zeros(noexch,1); %we will keep the epsilons
    EpsilonsK_min = zeros(noexch,1);
    
    k = k + 1;
    fprintf('Iteration: %i.\n', k);
    switch Modality
        case 'full' %------------------------------------------------------
            jx      = noexch; % size(P_NT,1)
            JMax    = zeros(size(S,2) + 1, jx);
            JMin    = zeros(size(S,2) + 1, jx);
            
            parfor i = 1:jx
                x   = P_NT(i,:); % Taking one row for each problem
                
                % We build both problems first
                LP1 = buildLPgurobi(S,x,w,lb,ub,bioIDX,ExRxnIDs,maxg,'Max'); %max
                LP2 = buildLPgurobi(S,x,w,lb,ub,bioIDX,ExRxnIDs,maxg,'Min'); %min
                
                % Solving 'max' problem first.
                solMax = gurobi(LP1, params);
                if strcmp(solMax.status,'OPTIMAL')
                    LP2.vbasis   = solMax.vbasis;
                    LP2.cbasis   = solMax.cbasis;
                    JMax(:,i)    = solMax.x;
                    EpsilonsK_max(i) = solMax.objval;
                    allSols      = [allSols sparse(solMax.x(ExRxnIDs))];
                end
                %Solvin 'min' problem with a heads-up, if possible. 
                solMin = gurobi(LP2, params);
                if strcmp(solMin.status,'OPTIMAL')
                    JMin(:,i)    = solMin.x;
                    EpsilonsK_min(i) = solMin.objval;
                    allSols      = [allSols sparse(solMin.x(ExRxnIDs))];
                end
            end
            EpsilonsK = [EpsilonsK_max' EpsilonsK_min'];
            Epsilons = [Epsilons EpsilonsK];
            Useless = abs(EpsilonsK) < eTol;
            EpsilonsK(Useless) = 0; % tolerance of Epsilon
            if all(EpsilonsK==0)
                break
            end
            J2 = [JMax, JMin];
            J2 = J2(:,~Useless);
        case 'fast' %------------------------------------------------------
    end
    [vopt, J2best] = VoptSelection(J2, ExRxnIDs, minBasis, vTol);
    minBasis = [minBasis, vopt];
    P_N = orth(minBasis)';
    if k==nits
        break
    end
end

% disp(ExRxnIDs)
% disp(bioIDX)
% disp(nullity)
% disp(FBA_init)
% disp(table(model.rxns(ExRxnIDs), minBasis))
% disp(P_NT)
% fprintf('Is it a truly agood matrix?: %f \n',P_N*P_N')
% fprintf('minimum growth: %i\n', maxg)
fprintf('Nº of exchages: %i\n', noexch)
% disp([EpsilonsK_max, EpsilonsK_min])
% disp(w)

CC = minBasis;
Output = w;

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

%% REQUIRED SUBROUTINES

function LPproblem = buildLPgurobi(S,P_NT,w,lb,ub,bioIDX,ExRxnIDs,maxg,sense)
% P_NT can be the entire matrix or just a row (x). 
% The 'w' vector is '1' for the 'full' version.

% Parameter initialization
[m,n]                   = size(S);
noexch                  = size(P_NT,1); % This should be '1' in the 'full'.
% noexch      = numel(ExRxnIDs);
% novar       = n + noexch; % No. of total variables
nores                   = m + noexch; % No. of constraints
lb(bioIDX)              = maxg; % Min growth forced

% Build the appropriate LP structure (minimal and maximal)
% Objective Function
% w_1*ep_1 + w_2*ep_2 +  ... + w_n*ep_n
LPproblem.obj           = [zeros(n,1); w.*ones(noexch,1)]; %Epsilons 
LPproblem.modelsense    = sense;

% Restrictions/Constraints
Proj                    = zeros(noexch, n);
Proj(:,ExRxnIDs)        = P_NT;
LPproblem.A             = sparse([S,     zeros(m,noexch);... % S·v = 0
                                  Proj, -eye(noexch)]);    % P·v ~= 0
LPproblem.lb            = [lb; -1e5*ones(noexch,1)]; % fluxes and epsilons
LPproblem.ub            = [ub;  1e5*ones(noexch,1)];

% Defining constraints sense
LPproblem.sense         = [repelem('=',m,1);... 
                           repelem('=',noexch,1)];

% Right-hand vector formulation
LPproblem.rhs           = zeros(nores,1);

%Extra parameters
end


% Subroutine to select among possible solutions for the k-th iteration
function [vopt, best] = VoptSelection(J2, ExRxnIDs, minBasis, vTol)
%We select which 
%fs = sum(any(J2)); %number of solutions found
fs = size(J2,2); % max n1 of posiblesolutions
vopt = J2(ExRxnIDs,:);
vopt(abs(vopt) < vTol) = 0;
if sum(any(J2))==1
    % If we only found one solution, we keep it.
    vopt = vopt(:,any(J2));
    best = J2(:,any(J2));
else
    % %we can try to implement this, to make it faster
%     if rank(vopt) == fs
%         % We also keep these solutions
%     else
%         if isempty(minBasis)
%             [~, I] = max(sum(vopt==0));
%             vopt = vopt(:,I);
%         else
%         end
%     end   
    c1 = mean(minBasis,2); %centroid
    cosSin = nan(fs,1);
    for i = 1:fs
        if ~any(J2(:,i))
            continue % we skip if there was no solution at i
        end
        vi = vopt(:,i);
        cosSin(i) = dot(vi,c1)/(norm(vi)*norm(c1)); %cosine similarity
    end
    [~, wI] = min(cosSin);
    vopt = vopt(:,wI);
    best = J2(:,wI);
end
end
