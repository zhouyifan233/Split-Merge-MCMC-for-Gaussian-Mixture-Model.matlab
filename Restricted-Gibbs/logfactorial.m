function [output] = logfactorial(inputarray1, inputarray2)
%LOGFACTORIAL Summary of this function goes here
%   Detailed explanation goes here
output = 0;

for i = 1:size(inputarray1, 2)
    output = output + log(inputarray1(i));
end
for i = 1:size(inputarray2, 2)
    output = output - log(inputarray2(i));
end

end

