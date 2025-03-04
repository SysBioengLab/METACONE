function M = addBiomassMet(model, metID, metName, bioIDX)
%% 
% Function that modifies a COBRA model and adds a virtual metabolite to the
% biomass reaction as EX_biomass_e orEX_biomass_(e)

% Add biomass metabolite
M = addMetabolite(model, metID, metName);

% Find growth reaction
if bioIDX == 0
    cans = model.rxns(contains(model.rxns, 'biomass', 'IgnoreCase', true));
    bioRxn = cans(1);
    bioIDX = findRxnIDs(M, bioRxn);
end

% Modifying S matrix
biomIDX = findMetIDs(M, metID);
M.S(biomIDX, bioIDX) = 1;

% Adding the exchange reaction
M = addReaction(M, 'reactionFormula', [metID ' -> ']);

end