assert(false,'This script is not for the direct excution.');

% This script includes example of 3-color CW folding analysis.
% The evaluation of maximum-likelihood  of 3-color CW is identical to that
% of 3-color ALEX. The analysis program can also be used for the ALEX data.

%% Input data

% Example data.
p = pwd; %put your installation directory

input_files{1}=[p '\Folding\trj_DA12.mat'];
input_files{2}=[p '\Folding\trj_DA1.mat'];
input_files{3}=[p '\Folding\trj_DA2.mat'];
input_files{4}=[p '\Folding\index_DA12.mat'];
input_files{5}=[p '\Folding\index_DA1.mat'];
input_files{6}=[p '\Folding\index_DA2.mat'];

input_data=input_bundle(input_files,'Folding');

burstbint3r=input_data{1};
cumindex=input_data{2};
indexone=input_data{3};
stateid=input_data{4};

exid=zeros(size(burstbint3r,1),1)+1; %this array indicates the excitation laser. 1: donor excitation, 2: acceptor1 excitation.
denatConc=indexone*0+1; %reserved parameter for denaturant conc. dependent reactions. Currently not in use.
%% pre-defined parameters 
chans=[3 4 1]; %[accpetor1 donor acceptor2]
concParams=[0.2448    0.0594    0.0200]; %[acceptor1->acceptor2 leakage (apparent), donor->acceptor1 leakage (apparent), donor->acceptor2 leakage (apparent)]

%% Parameter initialization: 2-state analysis

% initparams: parameter initialization
% LUbounds: lower and upper bounds
% fixed: parameter fixation. 1: fixed, 0: non-fixed.

%initparams{1}=[eff1s1 ... eff1sn ... eff2s1 ... eff2sn] for DA1A2
initparams{1}=[0.1549 0.2638 0.7107 0.5337]'; %for 2-state folding, [E1f E1u E2f E2u];
LUbounds{1}=[0.14 0.17; 0.24 0.28; 0.67 0.73; 0.51 0.56];
fixed{1}=[0 0 0 0]';
assert(numel(initparams{1})==numel(fixed{1}));assert(numel(initparams{1})*2==numel(LUbounds{1}));

%initparams{2}=[eff1s1 ... eff1sn leakage] for DA1
initparams{2}=[0.8938 0.7751 concParams(1)]';
LUbounds{2}=[0.85 0.94; 0.76 0.81; 0 0.4];
fixed{2}=[0 0 1]';
assert(numel(initparams{2})==numel(fixed{2}));assert(numel(initparams{2})*2==numel(LUbounds{2}));

%initparams{3}=[eff2s1 ... eff2sn] for DA2
initparams{3}=[0.2999 0.2147]';
LUbounds{3}=[0.26 0.32; 0.17 0.24];
fixed{3}=[0 0]';
assert(numel(initparams{3})==numel(fixed{3}));assert(numel(initparams{3})*2==numel(LUbounds{3}));

%initparams{4}=[eff12s1 ... eff12sn] for A1A2 (acceptor
%excitation)
initparams{4}=[0.8377 0.6554]';
LUbounds{4}=[0.75 0.9; 0.55 0.7];
fixed{4}=[0 0]'+1; %we don't fit these parameters because they are for ALEX
assert(numel(initparams{4})==numel(fixed{4}));assert(numel(initparams{4})*2==numel(LUbounds{4}));

%initparams{5}=[eff1D eff2D] for D-only
initparams{5}=[concParams(2) concParams(3)]';
LUbounds{5}=[0 0.1; 0 0.1;];
fixed{5}=[1 1]'; % these are pre-determined parameters
assert(numel(initparams{5})==numel(fixed{5}));assert(numel(initparams{5})*2==numel(LUbounds{5}));

%initparams{6}=[effA2 gamma] for A2 only (A1 excitation).
%gamma=A2onlyCountrate/A1A2countrate.
%These parameters are reserved for future usage. Currently not in use.
initparams{6}=[0.9 0.99]';
LUbounds{6}=[0.8 1;0 1];
fixed{6}=[1 1]'; %gamma should be fixed, always, However avoid effA2=1 which makes program crash sometimes
assert(numel(initparams{6})==numel(fixed{6}));assert(numel(initparams{6})*2==numel(LUbounds{6}));

%initparams{7}=[k(1) ~ k(n-1)]  k(1) is the sum of foward and backward rates between state 1 and 2;
initparams{7}=[1.1]';
LUbounds{7}=[0.9 1.5;];
fixed{7}=[0]';
assert(numel(initparams{7})==numel(fixed{7}));assert(numel(initparams{7})*2==numel(LUbounds{7}));

%initparams{8}=[f(1) ~ f(n-1)] p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
%Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
initparams{8}=[0.6479]';
LUbounds{8}=[0.6 0.7;];
fixed{8}=[0]';
assert(numel(initparams{8})==numel(fixed{8}));assert(numel(initparams{8})*2==numel(LUbounds{8}));

%initparams{9}=[kb1 pb1 kb2 pb2] blinking description for A1 and A2.
initparams{9}=[5.84 0.9884 1.07 0.9967]'; %
LUbounds{9}=[0 20;0.96 1; 0.5 1.5;0.98 1];

initparams{9}=[900 0.99999 900 0.99999]'; %for 2-state without blinking. Simple put high k value and p value close to 1 to avoid dark state.
LUbounds{9}=[0 1000;0.96 1; 0.5 1000;0.98 1];
fixed{9}=[1 1 1 1]'; %fix for a model without blinking.

%Donly background by A1 excitation
%Reserved parameters future usage in ALEX analysis. Currently not in use.. 
initparams{10}=[0.4 0.1]'; %[EA1 EA2]
LUbounds{10}=[0 1; 0 1];
fixed{10}=[1 1]'; %don't fit these parameters

%% Maximum-likelihood fitting

%without GPU
resparams=mlhrateeiglingen3cAlexABIC(initparams,LUbounds,burstbint3r,cumindex,indexone, fixed, stateid,exid,chans,denatConc);

%with GPU
%maxThread=13; %total number of CPU threads + 1 (1 additional thread for GPU management)
%GPU_load=0.1; %determine what portion of photon trajectories will be analyzed on GPU
%resparams=mlhrateeiglingen3c1Ex_G(initparams,LUbounds,burstbint3r,cumindex,indexone, fixed, stateid,exid,chans,denatConc,maxThread,GPU_load);

%% error estimation

% Covariance matrices are used to estimate errors.
drvintv=0.005; %relative perturbation for covariance matrices
drvpnt=5;  % number of perturbations on each side. Total 2*drvpnt perturbations will be evaluated around the determined parameters.

%without GPU
[errorparams, logmlh, bic]=mlhrateeiglingen3cAlexABICerror(resparams, LUbounds,burstbint3r,cumindex,indexone, fixed, stateid,exid,chans,[],drvintv,drvpnt);

%with GPU
%[errorparams, logmlh, bic]=mlhrateeiglingen3c1Ex_Gerror(resparams, LUbounds,burstbint3r,cumindex,indexone, fixed, stateid,exid,chans,[],drvintv,drvpnt,maxThread,GPU_load);

%% n-linear models

% For n-linear models simply put additional parameters in to the following 
% parameter sets.

% initparams{1}=[eff1s1 ... eff1sn ... eff2s1 ... eff2sn] for DA1A2
% initparams{2}=[eff1s1 ... eff1sn leakage] for DA1
% initparams{3}=[eff2s1 ... eff2sn] for DA2
% initparams{4}=[eff12s1 ... eff12sn] for A1A2 (acceptor excitation). Even
% though initparams{4} is a dummy in CW analysis, you need to match the
% number of parameters.
% initparams{7}=[k(1) ~ k(n-1)]  k(1) is the sum of foward and backward rates between state 1 and 2;
% initparams{8}=[f(1) ~ f(n-1)] p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
% Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
