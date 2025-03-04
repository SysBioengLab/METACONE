function output = isCobraModel(model) 
% Complementary Subroutine for metaCone and QMatrix creation to check if 
%'model' arg has proper fields ====================

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
