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
% -- May 2025.
%
% USAGE:
%
%   [QModel] = QModelCreation(models, Names, varargin)
%
% INPUTS:
%
%   models       : A (q,1) array of COBRA models with the required fields:
%       * S      : m x n Stoichiometric Matrix
%       * c      : n x 1 linear objective coefficients
%       * lb     : n x 1 lower bounds on flux vector (the variables)
%       * ub     : n x 1 upper bounds on flux vector
%   Names        : A (q,1) cell array with the names of the models.
%
% OPTIONAL INPUTS:
%
%   preCC        : A (q,2) struct with the metaCone results for each model
%       * preCC(:,1) column should be the C_ext (matrices) ordered as Names
%       * preCC(:,2) column should be the Output (structs) ordered as Names
%   Alpha        : a double between 0 and 1. 
%   Modality     : 'full' or 'fast' character string.
%   biomassIndex : A vector of biomass exchanges indices (integers)
%   Exchanges    : array (vector) with the exchanges (indices) of the model
%   vTol         : zero flux tolerance (double)
%   eTol         : zero epsilon tolerance (double)
%   Nullity      : false (default) to not calculate elementary matrix
%
% OUTPUTS:
%
%   QModel       : A struct with the main fields of a COBRA Model
%   CC           :
%   Output       : A structure with the following fields:
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


% Arguments Extraction and Initialization ---
models   = p.Results.models;
q        = numel(models);
Names    = p.Results.Names;
vTol     = p.Results.vTol;
eTol     = p.Results.eTol;
Nullity  = p.Results.Nullity;
CCRs     = p.Results.preCC;
Alpha    = p.Results.Alpha;
bioIDX   = p.Results.biomassIndex;
Modality = p.Results.Modality;

% Finding Biomass reaction ---
if numel(bioIDX) == 1
    foundBios = zeros(q,1); % indices of biomass rxns
    
    % We create the exception in case biomass rxn is not present.
    e1Mess      = 'Some models without biomass rxn. Check cause of MException.last';
    e1Ide       = 'QModelCreation:incompleteArguments';
    e1Exception = MException(e1Ide, e1Mess);
    for i = 1:q
        model        = models{i};
        temBios      = find(contains(model.rxnNames, 'biomass','IgnoreCase',true));
        foundBios(i) = temBios(1);
        if foundBios(i) == 0
            disp('Oh Shit!')
            mid_i       = 'QModelCreation:ModelWithoutBiomassRxn';
            msg_i       = sprintf('%s model without biomass rxn',Names{i});
            noBio       = MException(mid_i, msg_i);
            e1Exception = addCause(e1Exception, noBio);
        end
    end
    if any(foundBios==0)
        throw(e1Exception);
    end
else
    foundBios = bioIDX;
end

%% BASES CALCULATION ===

CCs = cell(q,1);
Res = cell(q,1);

disp(bioIDX)
disp(foundBios)

if isempty(CCRs)
    fprintf('No bases detected as argument.\n')
    for i = 1:q
        fprintf('Initializing metaCone() routine of %s\n',Names{i})
        [CCs{i}, Res{i}] = metaCone(models{i},...
            'Alpha',Alpha',...
            'biomassIndex',foundBios(i),...
            'vTol',vTol,...
            'eTol',eTol,...
            'Nullity',Nullity,...
            'Modality',Modality);
        fprintf('C_ext basis of %s completed\n',Names{i})
    end
else
    CCs = CCRs(:,1);
    Res = CCRs(:,2);
end

%% Q-MATRIX CALCULATION ===

try
    disp("Started Q-Matrix Calculation")
    [Q, Output] = QMatrixCreation(Names, CCs, Res);
catch ME
    QModel = {CCs; Res};
    disp(ME)
    return
end
allExcs = Output.allExcs;

%% Q-MODEL CONSTRUCTION ===
% The model has the structure of a COBRA model.

% We create the template of the model
QModel      = createModel();
QModel.S    = Q;
QModel.c    = zeros(size(Q,2),1);
QModel.b    = zeros(size(Q,1),1);
QModel.lb   = repmat(-100, size(Q,2),1);
QModel.ub   = repmat(100, size(Q,2),1);

%We set the names of the columns of the Q Matrix
QModel.rxns = {};
nparts      = q + size(Q,1);
Ctemp       = [];
for i = 1:q
    numccs          = size(CCs{i},2);
    QModel.rxns     = [QModel.rxns;   strcat(join([repmat({'a_'},numccs,1) split(num2str(1:numccs))],''),Names{i})];
    CtempBlock      = zeros(numccs,nparts);
    CtempBlock(:,i) = ones(numccs,1);
    Ctemp           = [Ctemp; CtempBlock];
end
QModel.rxns = [QModel.rxns; allExcs];
% numel(QModel.rxns) == size(Q,2)

%Participation Matrix
Ctemp               = [Ctemp; zeros(size(Q,1),q) eye(size(Q,1))];
if ~and(size(Ctemp,1)==size(Q,2), size(Ctemp,2)==nparts)
    error('ERROR! Cp matrix was not properly built')
end
QModel.Cp = Ctemp;


% Here, one must define which metabolites can only "enter", which can only
% go out, and which can do both. 
tempExLowerBounds = zeros(size(Q,1),q);
tempExUpperBounds = tempExLowerBounds;
for i = 1:length(models)
    model = models{i};
    exchanges = logical(findExcRxns(model)); %exchanges indexes

    %Forward
    forward = and(model.lb>=0, model.ub > 0); %forward reactions
    fExrxns = model.rxns(logical(exchanges.*forward)); %forward exch rxns
    tempExLowerBounds(:,i) = ismember(allExcs, fExrxns);

    %Reverse
    reverse = and(model.lb< 0, model.ub<= 0); %reverse reactoins
    rExrxns = model.rxns(logical(exchanges.*reverse)); %reverse exch rxns
    tempExUpperBounds(:,i) = ismember(allExcs, rExrxns);
end

forIrrev = sum(tempExLowerBounds, 2) == numel(models);
QModel.lb(end-(size(forIrrev,1)-1):end) = (forIrrev==0)*(-1e3);
revIrrev = sum(tempExUpperBounds, 2) == numel(models);
QModel.ub(end-(size(forIrrev,1)-1):end) = (revIrrev==0)*(1e3);
%QModel.lb(1:(end-size(Q,1))) = -100;
QModel.lb(1:(end-size(Q,1))) = -0;

QModel.mets = extractAfter(allExcs, 'EX_');
QModel.rev = and(QModel.lb < 0, QModel.ub > 0);

end
