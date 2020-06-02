function [centile_per_element]= BasicFunc_compute_centile_per_array(in_array)

%% Summary:
% function to compute percentile for each array element
%
% Input:  in_array: stat array
% Output: centile_per_element: centile values per element 
%                              (range: 1 to 100)
%
%%

% Compute centile
n_elem = length(in_array);
centile_per_element = zeros(n_elem,1);

for i=1:n_elem
    nless = sum(in_array < in_array(i));
    nequal = sum(in_array == in_array(i));
    centile_per_element(i,1) = 100 * (nless + 0.5*nequal) / n_elem;
end

end