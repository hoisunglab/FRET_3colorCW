function resparams = mlhrateeiglingen3cFolding(init_params,LU_bounds,input_data,gpu_mode)
%mlhrateeiglingen3cFolding Excute maximum-likelihood analysis depending on its
%gpu availability
%   Detailed explanation goes here

chans=[3 4 1];

burstbint3r=input_data{1};
cumindex=input_data{2};
indexone=input_data{3};
stateid=input_data{4};
exid=zeros(size(burstbint3r,1),1)+1;
concParams=[0.2448    0.0594    0.0200];

%initparams{1}=[eff1s1 ... eff1sn ... eff2s1 ... eff2sn] for DA1A2
initparams{1}=[init_params(1) init_params(5) init_params(2) init_params(6)]'; %for a3D, [E1f E1u E2f E2u];
LUbounds{1}=[LU_bounds(1,1) LU_bounds(1,2); LU_bounds(5,1) LU_bounds(5,2); LU_bounds(2,1) LU_bounds(2,2); LU_bounds(6,1) LU_bounds(6,2)];
fixed{1}=[0 0 0 0]';

%initparams{2}=[eff1s1 ... eff1sn leakage] for DA1
initparams{2}=[init_params(3) init_params(7) concParams(1)]';
%initparams{2}=[0.95 0.34 eff12D]';
LUbounds{2}=[LU_bounds(3,1) LU_bounds(3,2); LU_bounds(7,1) LU_bounds(7,2); 0 0.4];
fixed{2}=[0 0 1]';
assert(numel(initparams{2})==numel(fixed{2}));assert(numel(initparams{2})*2==numel(LUbounds{2}));

%initparams{3}=[eff2s1 ... eff2sn] for DA2
initparams{3}=[init_params(4) init_params(8)]';
LUbounds{3}=[LU_bounds(4,1) LU_bounds(4,2); LU_bounds(8,1) LU_bounds(8,2)];
fixed{3}=[0 0]';
assert(numel(initparams{3})==numel(fixed{3}));assert(numel(initparams{3})*2==numel(LUbounds{3}));

%initparams{4}=[eff12s1 ... eff12sn] for A1A2 (acceptor
%excitation)
initparams{4}=[0.8377 0.6554]';
LUbounds{4}=[0.75 0.9; 0.55 0.7];
fixed{4}=[0 0]'+1;
assert(numel(initparams{4})==numel(fixed{4}));assert(numel(initparams{4})*2==numel(LUbounds{4}));

%initparams{5}=[eff1D eff2D] for D-only
initparams{5}=[concParams(2) concParams(3)]';
LUbounds{5}=[0 0.1; 0 0.1;];
fixed{5}=[1 1]';
assert(numel(initparams{5})==numel(fixed{5}));assert(numel(initparams{5})*2==numel(LUbounds{5}));

%initparams{6}=[effA2 gamma] for A2 only (A1 excitation).
%gamma=A2onlyCountrate/A1A2countrate
%initparams{6}=[concParams{iConc}(1) concParams{iConc}(2)]';
initparams{6}=[0.9 0.99]';
LUbounds{6}=[0.8 1;0 1];
fixed{6}=[1 1]'; %gamma should be fixed, always, setting effA2=1 make program crash sometimes
assert(numel(initparams{6})==numel(fixed{6}));assert(numel(initparams{6})*2==numel(LUbounds{6}));

%initparams{7}=[k(1) ~ k(n-1)]  k(1) is the sum of foward and backward rates between state 1 and 2;
initparams{7}=[init_params(9)]';
LUbounds{7}=[LU_bounds(9,1) LU_bounds(9,2);];
fixed{7}=[0]';
assert(numel(initparams{7})==numel(fixed{7}));assert(numel(initparams{7})*2==numel(LUbounds{7}));

%initparams{8}=[f(1) ~ f(n-1)] p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
%Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
initparams{8}=[init_params(10)]';
LUbounds{8}=[LU_bounds(10,1) LU_bounds(10,2);];
fixed{8}=[0]';
assert(numel(initparams{8})==numel(fixed{8}));assert(numel(initparams{8})*2==numel(LUbounds{8}));

initparams{9}=[900 0.99999 900 0.99999]'; %for 2-state
LUbounds{9}=[0 1000;0.96 1; 0.5 1000;0.98 1];
fixed{9}=[1 1 1 1]';

%Donly background by A1 excitation
initparams{10}=[0.4 0.1]'; %[EA1 EA2]
LUbounds{10}=[0 1; 0 1];
fixed{10}=[1 1]';

%resparams=mlhrateeiglingen3cFolding(init_params,LU_bounds,input_data,analysis_mode{2});
if gpu_mode==0
    resparams=mlhrateeiglingen3cAlexABIC(initparams,LUbounds,burstbint3r,cumindex,indexone, fixed, stateid,exid,chans);
elseif gpu_mode==1
    denatConc=indexone*0+1;
    maxThread=13;
    GPU_load=0.1;
    resparams=mlhrateeiglingen3c1Ex_G(initparams,LUbounds,burstbint3r,cumindex,indexone, fixed, stateid,exid,chans,denatConc,maxThread,GPU_load);
else
    %
end
param_index=[1 3 5 8 2 4 6 9 16 17];
init_params(1:10)=resparams(param_index);
resparams=init_params;
end

