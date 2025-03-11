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