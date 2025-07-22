%% E. coli core modification model
% The model needs to include a biomass meatbolite for the MetaCone
% algorithm to work on it.

initCobraToolbox
changeCobraSolver('gurobi', 'ALL')
load ecoli_core_model.mat

bmodel = addMetabolite(model, 'biomass[c]', 'Biomass metabolite');
bmodel = addMetabolite(bmodel, 'biomass[e]', 'Biomass metabolite Extracellular');
bmodel = addReaction(bmodel, 'balance1', 'reactionName', 'Creation of extracellular biommas', 'reactionFormula', 'biomass[c] -> biomass[e]');
bmodel = addReaction(bmodel, 'EX_bio[e]', 'reactionName', 'Export of biomass', 'reactionFormula', 'biomass[e] -> ');
bmodel.S(73, 13) = 1; %The biomass rxns is modified
printRxnFormula(bmodel, bmodel.rxns(13));

save('./bmodel.mat', 'bmodel');
