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
default_vTol        = 1e-8; 
default_eTol        = 1e-10; %1e-7 works for ecoli core
default_Alpha       = 0.1;
default_Nullity     = false;
default_bioIDX      = 0;
default_Exchanges   = true;
default_Modality    = 'fast';
default_keepAll     = false;

% Check Functions ---
isvTolOk            = @(vTol) and(vTol < .5, vTol > 0);
isAlphaOk           = @(a) and(isnumeric(a), a >= 0 & a <= 1);
isNullityvalid      = @(nullity) ismember(nullity,[true, false]);
isbioIDXvalid       = @(bioIDX) class(bioIDX) == "double";
isExchsvalid        = @(exchs) or(iscellstr(exchs),isnumeric(exchs));
isModalityvalid     = @(mod) ismember(mod, {'full','fast'});
isBool              = @(b) or(class(b)=='logical',ismember(b, [0 1]));
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
addParameter(p, 'keepAll'   , default_keepAll,   isBool);

parse(p, model, varargin{:})

fprintf('Call: \n')
disp(p.Results)

%% INITIALIZATION ===


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
keepAll  = p.Results.keepAll;

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
noexch              = length(ExRxnIDs);
FBA_init            = optimizeCbModel(changeObjective(model, model.rxns(bioIDX)));
minBasis            = FBA_init.x(ExRxnIDs); % final solutions
P_N                 = orth(minBasis)'; % Initial value of Projection matrix
maxg                = Alpha*FBA_init.f; % minimum growth per conversion
if Modality == "fast"
    % w = randn(noexch,1);
    w               = randi(1e5, [1,noexch]) + rand(1,noexch); 
    % w = 1 + rand(1, noexch);
else % 'full'
    w               = 1;
end

%% GREEDY LP ITERATIONS ====
fprintf('Starting %s modality of the algorithm: \n', Modality)
fprintf('Minimum growth per conversion set to %f: \n', maxg)
k = 0;

% Complementary Performance ---
Epsilons            = [];
allSols             = [];
params.OutputFlag   = 0;
params.LPWarmStart  = 1;

tic
while true
    % PROJECTION MATRIX UPDATE .................
    P_NT = eye(noexch) - P_N'*P_N;
    
    % Preallocation vars for termination conditions.
    EpsilonsK_max = zeros(noexch,1);
    EpsilonsK_min = zeros(noexch,1);
    k = k + 1;
    
    fprintf('Iteration: %i.\n', k);
    
    % LP OPTIMIZATION .................
    switch Modality
        case 'full' %------------------------------------------------------
            jx      = noexch; % size(P_NT,1)
            JMax    = zeros(size(S,2) + 1, jx);
            JMin    = zeros(size(S,2) + 1, jx);
            
            parfor i = 1:jx
                x   = P_NT(i,:); % taking one row for each problem
                
                % We build both problems first ---
                LP1 = buildLPgurobi(S,x,w,lb,ub,bioIDX,ExRxnIDs,maxg,'Max'); %max
                LP2 = buildLPgurobi(S,x,w,lb,ub,bioIDX,ExRxnIDs,maxg,'Min'); %min
                
                % Solving 'max' problem first ---
                solMax = gurobi(LP1, params);
                if strcmp(solMax.status,'OPTIMAL')
                    LP2.vbasis       = solMax.vbasis;
                    LP2.cbasis       = solMax.cbasis;
                    JMax(:,i)        = solMax.x;
                    EpsilonsK_max(i) = solMax.objval;
                    if keepAll %only saved if flagged
                        allSols      = [allSols sparse(solMax.x(ExRxnIDs))];
                    end
                end
                %Solvin 'min' problem with a heads-up, if possible. ---
                solMin = gurobi(LP2, params);
                if strcmp(solMin.status,'OPTIMAL')
                    JMin(:,i)        = solMin.x;
                    EpsilonsK_min(i) = solMin.objval;
                    if keepAll
                        allSols      = [allSols sparse(solMin.x(ExRxnIDs))];
                    end
                end
            end
            
            % We keep all the epsilons ---
            EpsilonsK          = [EpsilonsK_max' EpsilonsK_min'];
            Epsilons           = [Epsilons EpsilonsK];
            
            % which rows gave null epsilons?
            Useless            = abs(EpsilonsK) < eTol;
            EpsilonsK(Useless) = 0; % tolerance of Epsilon
            
            % TERMINATION .................
            if all(EpsilonsK==0)
                break
            elseif k == noexch
                break
            end
            J2 = [JMax, JMin];
            J2 = J2(:,~Useless);
        case 'fast' %------------------------------------------------------
            % Aux. var. for preaallocation. ---
            J2      = zeros(size(S,2) + noexch,2);
            
            % Building both problems with the complete P_NT ---
            LP1     = buildLPgurobi(S,P_NT,w',lb,ub,bioIDX,ExRxnIDs,maxg,'Max'); %max
            LP2     = buildLPgurobi(S,P_NT,w',lb,ub,bioIDX,ExRxnIDs,maxg,'Min'); %min
            
            % Annotating epsilon indices.---
            eps_ids = (size(S,2)+1):size(LP1.A,2);
            
            % Solving the 'max' problem ---
            solMax  = gurobi(LP1, params);
            if strcmp(solMax.status,'OPTIMAL')
                LP2.vbasis         = solMax.vbasis;
                LP2.cbasis         = solMax.cbasis;
                J2(:,1)            = solMax.x;
                EpsilonsK_max(:,1) = solMax.x(eps_ids);
                if keepAll
                    allSols        = [allSols sparse(solMax.x(ExRxnIDs))];
                    Epsilons       = [Epsilons sparse(solMax.x(eps_ids))];
                end
            end
            % Solving the 'min' problem with warm-up ---
            solMin = gurobi(LP2, params);
            if strcmp(solMin.status,'OPTIMAL')
                J2(:,2)            = solMin.x;
                EpsilonsK_min(:,1) = solMin.x(eps_ids);
                if keepAll
                    allSols        = [allSols sparse(solMin.x(ExRxnIDs))];
                    Epsilons       = [Epsilons sparse(solMax.x(eps_ids))];
                end
                
            end
            
            % We retain all the epsilons ---
            EpsilonsK              = [EpsilonsK_max EpsilonsK_min];
            EpsilonsK(abs(EpsilonsK) < eTol) = 0; % tolerance of Epsilon
            
            % which LP gave null epsilsons? ---
            Useless                = ~any(EpsilonsK); 
            
            % TERMINATION .................
            if all(EpsilonsK==0)
                disp("All Epsilons are too small")
                break
            elseif k == noexch
                break
            end
            
            % We filter out the nullepsilons LPs. ---
            J2 = J2(:,~Useless);
            %--------------------------------------------------------------
    end
    % Selecting and saving the best solution of the iteration. ---
    [vopt, J2best] = VoptSelection(J2, ExRxnIDs, minBasis, vTol);
    minBasis       = [minBasis, vopt];
    
    % PROJECTION MATRIX UPDATE .................
    P_N            = orth(minBasis)';
end
runtime            = toc;

%% OUTPUT ===


% Additional information
Output.runtime   = runtime;

end % of MetaCone function
