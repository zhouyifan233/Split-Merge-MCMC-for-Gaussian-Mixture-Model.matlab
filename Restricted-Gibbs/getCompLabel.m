function CompLabel = getCompLabel(Models, idx)
    CompLabel = [];
    for i = 1:size(Models, 2)
        if Models{i}.SearchIdx(idx)
            CompLabel = [CompLabel, i];
        end
    end
    
    if isempty(CompLabel)
       error('The sample does not have a class...');
    elseif numel(CompLabel) > 1
        error('The sample has more than one class...');
    end
    
end
