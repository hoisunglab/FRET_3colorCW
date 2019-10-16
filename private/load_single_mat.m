function output=load_single_mat(input_mat)
% Directly returns variable when a MAT file contains only one variable
variables=load(input_mat);
names=fieldnames(variables);
if isempty(names) || numel(names)>1
    assert(false,'There are more than one variable');
end

output=variables.(names{1});

end