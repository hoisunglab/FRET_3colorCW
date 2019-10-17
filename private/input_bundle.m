function input_file_bundle=input_bundle(input_files,analysis_mode)

if strcmpi(analysis_mode,'Folding')
    trj_index=[1 2 3];
    index_index=[4 5 6];
    
    burstbint3rall=[];
    
    for iFile=1:3 % the first photon trajectories
        if isempty(input_files{trj_index(iFile)})
            errordlg('An input file is empty');
            assert(false,'An input file is empty');
        end
        burstbint3rall=[burstbint3rall;load_single_mat(input_files{trj_index(iFile)})];
    end
    
    cumindex=[0];
    stateid=[];
    states_order=[1 2 3]; %1: DA1A2, 2:DA1, 3:DA2
    for iFile=1:3
        if isempty(input_files{index_index(iFile)})
            errordlg('An input file is empty');
            assert(false,'An input file is empty');
        end
        cumindex_temp=load_single_mat(input_files{index_index(iFile)});
        cumindex_temp=cumindex_temp(2:end)+cumindex(end);
        cumindex=[cumindex;cumindex_temp];
        stateid=[stateid;zeros(numel(cumindex_temp),1)+states_order(iFile)];
        %indexone{iFile}=(1:(numel(cumindex{iFile})-1))';
    end
    indexone=(1:(numel(cumindex)-1))';
    input_file_bundle={burstbint3rall, cumindex, indexone, stateid};
elseif strcmpi(analysis_mode,'Binding')
    burstbint3rall=cell(1,2);
    cumindex=cell(1,2);
    indexone=cell(1,2);
    
    trj_index=[1 3];
    index_index=[4 6];
    for iFile=1:2 % photon trajectories
        if isempty(input_files{trj_index(iFile)})
            errordlg('An input file is empty');
            assert(false,'An input file is empty');
        end
        burstbint3rall{iFile}=load_single_mat(input_files{trj_index(iFile)});
    end
    for iFile=1:2
        if isempty(input_files{index_index(iFile)})
            errordlg('An input file is empty');
            assert(false,'An input file is empty');
        end
        cumindex{iFile}=load_single_mat(input_files{index_index(iFile)});
        indexone{iFile}=(1:(numel(cumindex{iFile})-1))';
    end
    input_file_bundle={burstbint3rall, cumindex, indexone};
else
    errordlg('Unexpected analysis mode');
    assert(false,'Unexpected analysis mode');
    input_file_bundle=[];
end