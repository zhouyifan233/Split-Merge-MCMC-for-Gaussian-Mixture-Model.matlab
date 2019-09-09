function NumsOfSample = getNumsOfSample(Models)
    NumsOfSample = [];
    for i = 1:size(Models, 2)
        NumsOfSample = [NumsOfSample, Models{i}.getNumSample()];
    end
end