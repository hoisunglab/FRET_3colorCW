function [outputArg1] = mlhrateeiglingen3cBinding(init_params,LU_bounds,input_data,gpu_mode)
%mlhrateeiglingen3cBinding Excute maximum-likelihood analysis depending on
%its gpu availability
%   Detailed explanation goes here

chans=[1 3 4]; %pre-defined.
isfastbinding=1; %this is for binding analysis.
stateidML=[1 3]; %in this analysis we assume that state1 is a bound state
ubid=2; %state2 is a unbound state
isA2bleach=1;    % True: If A2 bleaches or blink in the bound state, transition kinetics is included in the model.
fixparam={[0.05 0.0387]', [0 0]', [0 0]', 0.2};

burstbint3rall=input_data{1};
cumindex=input_data{2};
indexone=input_data{3};

%param_index=[1 5 2 6 3 11 12 9 10]'; %rearrange the params
param_index=[1 5 2 6 3 4 8 11 12 9 10 13 15 14 16]';
initparams_new=init_params(param_index);
LUbounds=LU_bounds(param_index,:);
if gpu_mode==0
    %resparams=mlhrateeiglingen3c1ExC(initparams_new,LUbounds,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans);
    resparams=mlhrateeiglingen3c1ExABIC(initparams_new,LUbounds,fixparam,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans);
elseif gpu_mode==1
    %Not yet implemented
    %resparams=mlhrateeiglingen3c1ExC(initparams_new,LUbounds,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans);
    resparams=mlhrateeiglingen3c1ExABIC(initparams_new,LUbounds,fixparam,burstbint3rall,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateidML,chans);
else
    errordlg('Unexpected gpu mode');
    assert(false,'Unexpected gpu mode');
end
init_params(param_index)=resparams;

outputArg1=init_params;
end

