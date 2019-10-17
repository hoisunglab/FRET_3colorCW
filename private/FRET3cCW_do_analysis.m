function output_params=FRET3cCW_do_analysis(init_params,LU_bounds,input_files,analysis_mode)
% input_data=cell(3,2);
% input_data{1}=DA12 photon trajectory
% input_data{2}=DA1 photon trajectory
% input_data{3}=DA2 photon trajectory
% input_data{4}=DA12 photon index
% input_data{5}=DA1 photon index
% input_data{6}=DA2 photon index

% analysis_mode{1}=folding or binding
% analysis_mode{2}=gpu on or off

input_data=input_bundle(input_files,analysis_mode{1});

if strcmpi(analysis_mode{1},'Folding')
    resparams=mlhrateeiglingen3cFolding(init_params,LU_bounds,input_data,analysis_mode{2});
elseif strcmpi(analysis_mode{1},'Binding')
    resparams=mlhrateeiglingen3cBinding(init_params,LU_bounds,input_data,analysis_mode{2});
else
    errordlg('Unexpected analysis mode');
    assert(false,'Unexpected analysis mode');
end

output_params=resparams;