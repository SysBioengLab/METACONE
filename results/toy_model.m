%% Toy model construction
% A small model with exchange reactions, a biomass proxy reaction and
% formulas for the metabolites was needed initially to test early versions
% of the MetaCone algorithm

clc; clear
initCobraToolbox(false)
changeCobraSolver('gurobi', 'ALL')

%% COBRA Model Creation

reactionIdentifiers = split('LDH_D L2M_Da ACOAHi GLXS MDH OAADC ME1 L_LACta O2t ACt2r CO2t H2Ot GLXte Ht')
exchange = split('EX_lac_L[e] EX_o2[e] EX_ac[e] EX_co2[e] EX_h2o[e] EX_glx[e] EX_h[e]');
reactionFormulas = {...
    'lac_D[c] + nad[c] <=> h[c] + nadh[c] + pyr[c]'
    'lac_D[c] + o2[c] <=> ac[c] + co2[c] + h2o[c]'
    'accoa[c] + h2o[c] -> ac[c] + coa[c] + h[c]' %rev?
    'coa[c] + h[c] + mal_L[c] <=> accoa[c] + glx[c] + h2o[c]'
    'mal_L[c] + nad[c] <=> h[c] + nadh[c] + oaa[c]'
    'h[c] + oaa[c] -> co2[c] + pyr[c]' %rev?
    'mal_L[c] + nad[c] -> co2[c] + nadh[c] + pyr[c]' %rev?
    'lac_D[e] <=> lac_D[c]' 
    'o2[e] <=> o2[c]'
    'ac[e] + h[e] <=> ac[c] + h[c]' 
    'co2[e] <=> co2[c]'
    'h2o[e] <=> h2o[c]'
    'glx[c] -> glx[e]' %rev?
    'h[e] <=> h[c]' 
    'lac_D[e] <=>'
    'o2[e] <=>'
    'ac[e] <=>'
    'co2[e] <=>'
    'h2o[e] <=>'
    'glx[e] <=>'
    'h[e] <=>'}
reactionNames = {...
    'D-Lactate Dehydrogenase'
    'D-Lactate 2-monooxygenase'
    'Acetyl Coenzyme A Hydrolase'
    'Glyoxylate synthase, reversible'
    'Malate Dehydrogenase'
    'Oxaloacetate Decarboxylase'
    'Malic enzyme (NAD)'
    'Lactate Transport by FictionalProtein'
    'O2 Transport (Diffusion)'
    'Acetate Reversible Transport via Proton Symport'
    'CO2 Transporter via Diffusion'
    'H2O Transport via Diffusion'
    'Transport of Glyoxylate, Extracellular'
    'Transport of Protons, fictional'
    'Exchange of L-Lactate'
    'Exchange of Oxygen'
    'Exchange of Acetate'
    'Exchange of Carbon Dioxide'
    'Exchange of Water'
    'Exchange of Glyoxylate'
    'Exchange of Proton'}
reactionIdentifiers = [reactionIdentifiers; exchange];
% rev = ~or(toymodel1.lb==0, toymodel1.ub==0);
rev = ones(length(reactionNames), 1); 
rev([3 7 13]) = 0; % TESTING AS IRREVERSIBLE

toymodel1 = createModel(reactionIdentifiers, reactionNames, reactionFormulas, 'revFlagList', rev);

% WE ADD CHARGES AND FORMULAS
MFC = {... 
    'pyr'       'C3H3O3'            -1     %
    'lac_D'       'C3H5O3'            -1     %  
    'nad'       'C21H26N7O14P2'     -1     %
    'nadh'      'C21H27N7O14P2'     -2     %
    'h'         'H1'                +1     %
    'o2'        'O2'                0     %
    'ac'        'C2H3O2'            -1     %
    'accoa'     'C23H34N7O17P3S1'   -4     %
    'co2'       'C1O2'              0     %
    'h2o'       'H2O1'              0     %
    'coa'       'C21H32N7O16P3S1'   -4     %
    'glx'       'C2H1O3'            -1     %
    'mal_L'     'C4H4O5'            -2     %
    'oaa'       'C4H2O5'            -2     % 
    };

mets = extractBefore(toymodel1.mets, '[')
toymodel1.metCharge = cell(length(mets),1);
toymodel1.metFormulas = cell(length(mets),1)
cargas = MFC(:,3); formulas = MFC(:,2);
for i = 1:length(mets)
    met_i = mets(i);
    ix = strcmp(met_i, MFC(:,1));
    toymodel1.metCharge(i) = cargas(ix);
    toymodel1.metFormulas(i)= formulas(ix);
end
clear mets met_i i ix MFC rev reactionFormulas reactionIdentifiers reactionNames exchange cargas formulas 

toymodel1_bio = addReaction(toymodel1, 'biomass', 'reactionFormula', 'pyr[c] + nadh[c] + accoa[c] -> biomass[c] + nad[c] + coa[c] + co2[c] + h2o[c]');
toymodel1_bio = addReaction(toymodel1_bio, 'EX_biomass[c]', 'reactionFormula', '1 biomass[c] ->');
toymodel1_bio = changeObjective(toymodel1_bio, 'EX_biomass[c]');
FBA_tm = optimizeCbModel(toymodel1_bio, 'max');
FBA_tm.f
FBA_tm.x

[minFlux, maxFlux] = fluxVariability(toymodel1_bio, 0)
bar(maxFlux)
hold on
bar(minFlux)
save('./toymodel1_bio', 'toymodel1_bio')
